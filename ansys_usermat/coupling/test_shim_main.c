/* test_shim_main.c — tiny driver used by tests/test_coupling_shim.py to exercise
 * biofilm_py_eval() against a running material_server.py.
 *
 * Reads from stdin (whitespace-separated):
 *     F(9)  Fv(9)  alpha C10 C01 D1 eta mtype dt
 * Prints stress(6), Fv_new(9), then dsdePl(36), one value per line (%.17g).
 * Exit code = the shim's return code (0 = success).
 *
 *     cc test_shim_main.c biofilm_py_eval.c -o test_shim
 */
#include <stdio.h>

int biofilm_py_eval(const double *F9, const double *Fv9, const double *params7,
                    double *stress6, double *Fvnew9, double *dsde36);
void biofilm_py_close(void);

int main(void)
{
    double F9[9], Fv9[9], params7[7], stress6[6], Fvnew9[9], dsde36[36];
    int i, rc;

    for (i = 0; i < 9; i++) if (scanf("%lf", &F9[i]) != 1) return 100;
    for (i = 0; i < 9; i++) if (scanf("%lf", &Fv9[i]) != 1) return 100;
    for (i = 0; i < 7; i++) if (scanf("%lf", &params7[i]) != 1) return 100;

    rc = biofilm_py_eval(F9, Fv9, params7, stress6, Fvnew9, dsde36);
    if (rc != 0) { biofilm_py_close(); return rc; }

    for (i = 0; i < 6; i++)  printf("%.17g\n", stress6[i]);
    for (i = 0; i < 9; i++)  printf("%.17g\n", Fvnew9[i]);
    for (i = 0; i < 36; i++) printf("%.17g\n", dsde36[i]);
    biofilm_py_close();
    return 0;
}
