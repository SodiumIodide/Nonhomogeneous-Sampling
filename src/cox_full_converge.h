#ifndef _COX_FULL_CONVERGE_H
#define _COX_FULL_CONVERGE_H

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>

#ifndef _RUNNING_STATISTICS_H
#include "running_statistics.h"
#endif

#ifndef _MESHMAP_H
#include "meshmap.h"
#endif

#ifndef _CONSTANTS_H
#include "constants.h"
#endif

#ifndef _LINSPACE_H
#include "linspace.h"
#endif

#ifndef _CONFIG_H
#include "config.h"
#endif

#ifdef S2
    #ifndef _S2_H
        #include "s2.h"
    #endif
    #define INNER_SOLVE() do {\
        s2(rng, &inner_phi_mat_stat_0, &inner_phi_mat_stat_1, chord_0, chord_1);\
    } while (0)
#elif defined PUREABS
    #ifndef _PURE_ABS_H
        #include "pure_abs.h"
    #endif
    #define INNER_SOLVE() do {\
        pure_abs(rng, &inner_phi_mat_stat_0, &inner_phi_mat_stat_1, chord_0, chord_1);\
    } while (0)
#elif defined ANALYTICAL
    #ifndef _ANALYTICAL_H
        #include "analytical.h"
    #endif
    #define INNER_SOLVE() do {\
        analytical(rng, &inner_phi_mat_stat_0, &inner_phi_mat_stat_1, chord_0, chord_1);\
    } while (0)
#elif defined VOLFRAC
    #ifndef _VOLUME_FRACTION_H
        #include "volume_fraction.h"
    #endif
    #define INNER_SOLVE() do {\
        volume_fraction(rng, &inner_phi_mat_stat_0, &inner_phi_mat_stat_1, chord_0, chord_1);\
    } while (0)
#endif

#ifdef COX_UNIFORM
    #define SAMPLE() do {\
        cox_sample_uniform(rng, &chord_0, &chord_1);\
    } while (0)
#elif defined COX_EXPONENTIAL
    #define SAMPLE() do {\
        cox_sample_exponential(rng, &chord_0, &chord_1);\
    } while (0)
#else
    #define SAMPLE() do {\
        cox_sample_gaussian(rng, &chord_0, &chord_1);\
    } while (0)
#endif

#define PI 4.0 * atan(1.0)

#define CLEANUP() do {\
    free(inner_phi_mat_stat_0);\
    free(inner_phi_mat_stat_1);\
} while (0)

void cox_sample_uniform(const gsl_rng* rng, double* chord_0, double* chord_1) {
    double delta_0 = COX_END_VALUE[0] - COX_START_VALUE[0];
    double delta_1 = COX_END_VALUE[1] - COX_START_VALUE[1];

    *chord_0 = COX_START_VALUE[0] + gsl_rng_uniform_pos(rng) * delta_0;
    *chord_1 = COX_START_VALUE[1] + gsl_rng_uniform_pos(rng) * delta_1;
}

void cox_sample_exponential(const gsl_rng* rng, double* chord_0, double* chord_1) {
    double mean_0 = (COX_END_VALUE[0] + COX_START_VALUE[0]) / 2.0;
    double mean_1 = (COX_END_VALUE[1] + COX_START_VALUE[1]) / 2.0;

    *chord_0 = -mean_0 * log(gsl_rng_uniform_pos(rng));
    *chord_1 = -mean_1 * log(gsl_rng_uniform_pos(rng));
}

void cox_sample_gaussian(const gsl_rng* rng, double* chord_0, double* chord_1) {
    double mean_0 = (COX_END_VALUE[0] + COX_START_VALUE[0]) / 2.0;
    double mean_1 = (COX_END_VALUE[1] + COX_START_VALUE[1]) / 2.0;

    double uniform[2], gaussian[2];
    double radius, angle;

    do {
        uniform[0] = gsl_rng_uniform_pos(rng);
        uniform[1] = gsl_rng_uniform_pos(rng);
        radius = sqrt(-2.0 * log(uniform[0]));
        angle = 2.0 * PI * uniform[1];
        gaussian[0] = radius * sin(angle);
        gaussian[1] = radius * cos(angle);

        *chord_0 = mean_0 + sqrt(GAUSSIAN_VARIANCE[0]) * gaussian[0];
        *chord_1 = mean_1 + sqrt(GAUSSIAN_VARIANCE[1]) * gaussian[1];
    } while ((*chord_0 <= 0.0) || (*chord_1 <= 0.0));
}

void cox_full_converge(const gsl_rng* rng, runningstat** overall_flux_0, runningstat** overall_flux_1) {
    // Iterator
    long int c;
    long int iterations_counter = 0;

    // Material chord lengths
    double chord_0, chord_1;

    // Structured placeholders
    runningstat *inner_phi_mat_stat_0 = malloc(sizeof *inner_phi_mat_stat_0 * NUM_CELLS);
    runningstat *inner_phi_mat_stat_1 = malloc(sizeof *inner_phi_mat_stat_1 * NUM_CELLS);
    for (c = 0; c < NUM_CELLS; c++) {
        *(inner_phi_mat_stat_0 + c) = init_runningstat();
        *(inner_phi_mat_stat_1 + c) = init_runningstat();
    }

    int cont_calc_outer = 1;

    while(cont_calc_outer) {
        // Reset running statistics
        for (c = 0; c < NUM_CELLS; c++) {
            clear(inner_phi_mat_stat_0);
            clear(inner_phi_mat_stat_1);
        }

        // Sample the chord lengths
        SAMPLE();

        INNER_SOLVE();

        for (c = 0; c < NUM_CELLS; c++) {
            push(&*((*overall_flux_0) + c), mean(&*(inner_phi_mat_stat_0 + c)));
            push(&*((*overall_flux_1) + c), mean(&*(inner_phi_mat_stat_1 + c)));
        }

        iterations_counter++;

        if (iterations_counter % NUM_SAY == 0) {
            printf("Overall Iteration Number: %ld / %ld\n", iterations_counter, NUM_COX_REALIZATIONS_OUTER);
        }
    }

    CLEANUP();
}

#endif
