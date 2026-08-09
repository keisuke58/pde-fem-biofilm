/* biofilm_py_eval.c — C shim connecting the Fortran USERMAT hook to the Python
 * material server (coupling/material_server.py).
 *
 * Implements the symbol declared by `usermat_py_hook.f`:
 *
 *     int biofilm_py_eval(const double *F9, const double *Fv9,
 *                         const double *params7, double *stress6,
 *                         double *Fvnew9, double *dsde36);
 *
 * Mechanism: a persistent local TCP connection to material_server.py, one
 * newline-delimited JSON frame per Gauss-point evaluation (see protocol.py).
 * The socket is opened lazily on first use and reused for the whole run — a
 * fresh connection per evaluation would dominate the cost.
 *
 * params7 = {alpha, C10, C01, D1, eta, mtype, dt}
 * Voigt order: Abaqus 11,22,33,12,13,23.
 * Returns 0 on success, nonzero on any failure — the Fortran side then falls
 * back to the verified inline core, so a dead server degrades to the reference
 * law rather than aborting the solve.
 *
 * Build (example):
 *     cc -c -fPIC biofilm_py_eval.c -o biofilm_py_eval.o
 * then link the .o together with the USERMAT objects.
 *
 * Host/port overridable via env: BIOFILM_PY_HOST, BIOFILM_PY_PORT.
 */
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <unistd.h>

#define RECV_CAP 65536

static int g_fd = -1;              /* persistent connection, -1 = not connected */

static int py_connect(void)
{
    const char *host = getenv("BIOFILM_PY_HOST");
    const char *port = getenv("BIOFILM_PY_PORT");
    struct sockaddr_in addr;
    int fd, one = 1;

    if (!host) host = "127.0.0.1";

    fd = socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) return -1;

    memset(&addr, 0, sizeof addr);
    addr.sin_family = AF_INET;
    addr.sin_port = htons((unsigned short)(port ? atoi(port) : 8765));
    if (inet_pton(AF_INET, host, &addr.sin_addr) != 1) { close(fd); return -1; }
    if (connect(fd, (struct sockaddr *)&addr, sizeof addr) != 0) { close(fd); return -1; }

    /* per-call latency matters far more than throughput here */
    setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &one, sizeof one);
    return fd;
}

static int send_all(int fd, const char *buf, size_t n)
{
    while (n > 0) {
        ssize_t w = write(fd, buf, n);
        if (w <= 0) return -1;
        buf += w; n -= (size_t)w;
    }
    return 0;
}

/* Read one newline-terminated frame. */
static int recv_line(int fd, char *buf, size_t cap)
{
    size_t used = 0;
    while (used + 1 < cap) {
        ssize_t r = read(fd, buf + used, 1);
        if (r <= 0) return -1;
        if (buf[used] == '\n') { buf[used] = '\0'; return 0; }
        used += (size_t)r;
    }
    return -1;                                   /* frame longer than the buffer */
}

/* Minimal JSON array extractor: finds "key":[ ... ] and parses `count` doubles.
 * The response schema is fixed and machine-generated (protocol.py), so a full
 * JSON parser would be overkill; anything unexpected is reported as failure. */
static int parse_array(const char *json, const char *key, double *out, int count)
{
    const char *p = strstr(json, key);
    int i;
    if (!p) return -1;
    p = strchr(p, '[');
    if (!p) return -1;
    p++;
    for (i = 0; i < count; i++) {
        char *end;
        out[i] = strtod(p, &end);
        if (end == p) return -1;
        p = end;
        while (*p == ' ' || *p == ',') p++;
    }
    return 0;
}

static int parse_scalar(const char *json, const char *key, double *out)
{
    const char *p = strstr(json, key);
    char *end;
    if (!p) return -1;
    p = strchr(p, ':');
    if (!p) return -1;
    *out = strtod(p + 1, &end);
    return (end == p + 1) ? -1 : 0;
}

int biofilm_py_eval(const double *F9, const double *Fv9, const double *params7,
                    double *stress6, double *Fvnew9, double *dsde36)
{
    char req[2048], resp[RECV_CAP];
    int n, i, attempt;
    double detFe;

    n = snprintf(req, sizeof req,
        "{\"F\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g],"
        "\"Fv\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g],"
        "\"alpha\":%.17g,\"C10\":%.17g,\"C01\":%.17g,\"D1\":%.17g,"
        "\"eta\":%.17g,\"mtype\":%.17g,\"dt\":%.17g}\n",
        F9[0], F9[1], F9[2], F9[3], F9[4], F9[5], F9[6], F9[7], F9[8],
        Fv9[0], Fv9[1], Fv9[2], Fv9[3], Fv9[4], Fv9[5], Fv9[6], Fv9[7], Fv9[8],
        params7[0], params7[1], params7[2], params7[3],
        params7[4], params7[5], params7[6]);
    if (n <= 0 || (size_t)n >= sizeof req) return 1;

    /* One reconnect retry: the server may have been restarted mid-run. */
    for (attempt = 0; attempt < 2; attempt++) {
        if (g_fd < 0) g_fd = py_connect();
        if (g_fd < 0) return 2;

        if (send_all(g_fd, req, (size_t)n) == 0 &&
            recv_line(g_fd, resp, sizeof resp) == 0)
            break;

        close(g_fd);
        g_fd = -1;
        if (attempt == 1) return 3;
    }

    if (strstr(resp, "\"error\"")) return 4;
    if (parse_array(resp, "\"stress\"", stress6, 6) != 0) return 5;
    if (parse_array(resp, "\"Fv_new\"", Fvnew9, 9) != 0) return 6;
    if (parse_array(resp, "\"dsdePl\"", dsde36, 36) != 0) return 7;
    if (parse_scalar(resp, "\"detFe\"", &detFe) != 0) return 8;

    for (i = 0; i < 6; i++)  if (stress6[i] != stress6[i]) return 9;   /* NaN guard */
    for (i = 0; i < 9; i++)  if (Fvnew9[i] != Fvnew9[i])  return 9;
    return 0;
}

/* Optional: close the connection at the end of a run. */
void biofilm_py_close(void)
{
    if (g_fd >= 0) { close(g_fd); g_fd = -1; }
}
