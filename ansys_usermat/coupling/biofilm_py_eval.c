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
 *     cc -c -fPIC biofilm_py_eval.c -o biofilm_py_eval.o        (Linux)
 *     gcc -c biofilm_py_eval.c -o biofilm_py_eval.o -lws2_32     (MinGW on
 *         Windows -- this particular MinGW-w64 build (portable WinLibs,
 *         see dev-env.ps1) does not ship arpa/inet.h etc., so it takes the
 *         same Winsock2 path as MSVC below, not the POSIX one; add
 *         -lws2_32 explicitly since MinGW does not honor MSVC's #pragma
 *         comment auto-link)
 *     cl /c biofilm_py_eval.c                                    (MSVC,
 *         from a vcvars64 shell; links ws2_32.lib automatically via the
 *         #pragma comment below)
 * then link the .o/.obj together with the USERMAT objects.
 *
 * Host/port overridable via env: BIOFILM_PY_HOST, BIOFILM_PY_PORT.
 *
 * Windows port, 2026-09-03: guarded on _WIN32 (true for both MSVC and
 * MinGW) rather than _MSC_VER, because this MinGW-w64 build lacks the
 * POSIX socket compatibility headers other MinGW distributions ship --
 * confirmed empirically (`arpa/inet.h: No such file or directory`), not
 * assumed. send()/recv() are used uniformly instead of read()/write() --
 * both platforms' socket layers support them identically, so that one
 * substitution is enough to avoid a second implementation of
 * send_all()/recv_line().
 */
#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "ws2_32.lib")
typedef SOCKET sock_t;
#define SOCK_INVALID INVALID_SOCKET
#define SOCK_CLOSE(s) closesocket(s)
#else
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <sys/socket.h>
#include <unistd.h>
typedef int sock_t;
#define SOCK_INVALID (-1)
#define SOCK_CLOSE(s) close(s)
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RECV_CAP 65536

static sock_t g_fd = SOCK_INVALID;   /* persistent connection */

#ifdef _WIN32
static int g_wsa_ready = 0;
static void ensure_wsa(void)
{
    if (!g_wsa_ready) {
        WSADATA wd;
        WSAStartup(MAKEWORD(2, 2), &wd);
        g_wsa_ready = 1;
    }
}
#endif

static sock_t py_connect(void)
{
    const char *host = getenv("BIOFILM_PY_HOST");
    const char *port = getenv("BIOFILM_PY_PORT");
    struct sockaddr_in addr;
    sock_t fd;
    int one = 1;

#ifdef _WIN32
    ensure_wsa();
#endif
    if (!host) host = "127.0.0.1";

    fd = socket(AF_INET, SOCK_STREAM, 0);
    if (fd == SOCK_INVALID) return SOCK_INVALID;

    memset(&addr, 0, sizeof addr);
    addr.sin_family = AF_INET;
    addr.sin_port = htons((unsigned short)(port ? atoi(port) : 8765));
    if (inet_pton(AF_INET, host, &addr.sin_addr) != 1) { SOCK_CLOSE(fd); return SOCK_INVALID; }
    if (connect(fd, (struct sockaddr *)&addr, sizeof addr) != 0) { SOCK_CLOSE(fd); return SOCK_INVALID; }

    /* per-call latency matters far more than throughput here */
    setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, (const char *)&one, sizeof one);
    return fd;
}

static int send_all(sock_t fd, const char *buf, size_t n)
{
    while (n > 0) {
        int w = send(fd, buf, (int)n, 0);
        if (w <= 0) return -1;
        buf += w; n -= (size_t)w;
    }
    return 0;
}

/* Read one newline-terminated frame. */
static int recv_line(sock_t fd, char *buf, size_t cap)
{
    size_t used = 0;
    while (used + 1 < cap) {
        int r = recv(fd, buf + used, 1, 0);
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
        if (g_fd == SOCK_INVALID) g_fd = py_connect();
        if (g_fd == SOCK_INVALID) return 2;

        if (send_all(g_fd, req, (size_t)n) == 0 &&
            recv_line(g_fd, resp, sizeof resp) == 0)
            break;

        SOCK_CLOSE(g_fd);
        g_fd = SOCK_INVALID;
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
    if (g_fd != SOCK_INVALID) { SOCK_CLOSE(g_fd); g_fd = SOCK_INVALID; }
}
