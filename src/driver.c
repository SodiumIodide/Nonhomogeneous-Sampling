#ifndef _CONFIG_H
#include <config.h>
#endif

#ifndef _RUNNING_STATISTICS_H
#include "running_statistics.h"
#endif

#ifndef _CONSTANTS_H
#include "constants.h"
#endif

#ifdef COX_FULL_CONVERGE
    #ifndef _COX_FULL_CONVERGE_H
        #include "cox_full_converge.h"
        #define SOLVE() do {\
            cox_full_converge(rng, &phi_mat_stat_0, &phi_mat_stat_1);
        }
    #endif
#elif defined S2
    #ifndef _S2_H
        #include "s2.h"
    #endif
    #define SOLVE() do {\
        s2(rng, &phi_mat_stat_0, &phi_mat_stat_1, 0.0, 0.0);\
    } while (0)
#elif defined PUREABS
    #ifndef _PURE_ABS_H
        #include "pure_abs.h"
    #endif
    #define SOLVE() do {\
        pure_abs(rng, &phi_mat_stat_0, &phi_mat_stat_1, 0.0, 0.0);\
    } while (0)
#elif defined GEOMETRYTIME
    #ifndef _GEOMETRY_TIME_H
        #include "geometry_time.h"
    #endif
    #define SOLVE() do {\
        geometry_time(rng);\
    } while (0)
#elif defined ANALYTICAL
    #ifndef _ANALYTICAL_H
        #include "analytical.h"
    #endif
    #define SOLVE() do {\
        analytical(rng, &phi_mat_stat_0, &phi_mat_stat_1, 0.0, 0.0);\
    } while (0)
#elif defined VOLFRAC
    #ifndef _VOLUME_FRACTION_H
        #include "volume_fraction.h"
    #endif
    #define SOLVE() do {\
        volume_fraction(rng, &phi_mat_stat_0, &phi_mat_stat_1, 0.0, 0.0);\
    } while (0)
#endif

#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    // Set and seed generator
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, SEED);

#ifndef GEOMETRY_TIME
    // Cell iterator
    long int c;

    // Solution variables - statistics flux variables
    runningstat *phi_mat_stat_0 = malloc(sizeof *phi_mat_stat_0 * NUM_CELLS);
    runningstat *phi_mat_stat_1 = malloc(sizeof *phi_mat_stat_1 * NUM_CELLS);

    for (c = 0; c < NUM_CELLS; c++) {
        *(phi_mat_stat_0 + c) = init_runningstat();
        *(phi_mat_stat_1 + c) = init_runningstat();
    }
#endif

    SOLVE();

#ifndef GEOMETRY_TIME
    printf("\nCalculation done\n");

    // File handle
    FILE *fp;
    // File variables
    double delta_x = END_DIST / (double)NUM_CELLS;
    double distance = 0.0;
    // Write file
    OPEN_FILE();
    fprintf(fp, "distance,flux0,varflux0,flux1,varflux1\n");
    for (c = 0; c < NUM_CELLS; c++) {
        distance += delta_x;
        fprintf(fp, "%f,%f,%f,%f,%f\n", distance - delta_x / 2.0, mean(&*(phi_mat_stat_0 + c)), variance(&*(phi_mat_stat_0 + c)), mean(&*(phi_mat_stat_1 + c)), variance(&*(phi_mat_stat_1 + c)));
    }
    printf("File written\n");

    // Cleanup
    gsl_rng_free(rng);
    fclose(fp);
    free(phi_mat_stat_0);
    free(phi_mat_stat_1);
#endif

    return 0;
}
