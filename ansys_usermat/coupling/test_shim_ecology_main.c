/* test_shim_ecology_main.c — tiny driver exercising biofilm_ecology_eval()
 * against a running material_server.py, mirroring test_shim_main.c.
 *
 * Reads from stdin (whitespace-separated): g(12) theta(20) dt_h
 * Prints g_new(12), one value per line (%.17g). Exit code = the shim's
 * return code (0 = success).
 *
 *     cc test_shim_ecology_main.c biofilm_py_eval.c -o test_shim_ecology
 */
#include <stdio.h>

int biofilm_ecology_eval(const double *g12, const double *theta20, double dt_h,
                         double *g_new12);
void biofilm_py_close(void);

int main(void)
{
    double g12[12], theta20[20], dt_h, g_new12[12];
    int i, rc;

    for (i = 0; i < 12; i++) if (scanf("%lf", &g12[i]) != 1) return 100;
    for (i = 0; i < 20; i++) if (scanf("%lf", &theta20[i]) != 1) return 100;
    if (scanf("%lf", &dt_h) != 1) return 100;

    rc = biofilm_ecology_eval(g12, theta20, dt_h, g_new12);
    if (rc != 0) { fprintf(stderr, "biofilm_ecology_eval rc=%d\n", rc); biofilm_py_close(); return rc; }

    for (i = 0; i < 12; i++) printf("%.17g\n", g_new12[i]);
    biofilm_py_close();
    return 0;
}
