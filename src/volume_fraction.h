#ifndef _VOLUME_FRACTION_H
#define _VOLUME_FRACTION_H

#include <stdio.h>

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

#ifndef _CONFIG_H
#include "config.h"
#endif

#define CLEANUP() do {\
    free(x_delta);\
    free(x_arr);\
    free(materials);\
    free(unity);\
    free(m_vfrac_0);\
    free(m_vfrac_1);\
} while (0)

void volume_fraction(const gsl_rng* rng, runningstat** s_vfrac_0, runningstat** s_vfrac_1) {
    // Iterators
    long int i, c;
    // Geometry variables
    // Initially allocate to 1-length
    double *x_delta = malloc(sizeof *x_delta * NUM_DIVS);
    double *x_arr = malloc(sizeof *x_arr * NUM_DIVS);
    int *materials = malloc(sizeof *materials * NUM_DIVS);
    long num_r_cells = 0;
    // Generic return
    int success = 0;
    double delta_x = END_DIST / (double)NUM_CELLS;
    // Solution variables
    double *unity = malloc(sizeof *unity * NUM_DIVS);
    // Structured results from map onto array
    double *m_vfrac_0 = malloc(sizeof *m_vfrac_0 * NUM_CELLS);
    double *m_vfrac_1 = malloc(sizeof *m_vfrac_1 * NUM_CELLS);
    // Buffer variable for reallocation
    double *buf_unity = NULL;

    for (i = 0; i < NUM_REALIZATIONS; i++) {
        // Fill Markovian geometry
        GEOMETRY_GEN();
        if (!success) {
            CLEANUP();
            printf("Failed to generate geometry\n");
            exit(EXIT_FAILURE);
        }

        buf_unity = realloc(unity, sizeof *unity * num_r_cells);
        if (buf_unity == NULL) {
            CLEANUP();
            printf("Error in reallocation: flux/unity!\n");
            exit(EXIT_FAILURE);
        }
        unity = buf_unity;
        buf_unity = NULL;

        // Initialize unity
        for (c = 0; c < num_r_cells; c++) {
            *(unity + c) = 1.0;
        }

        // Map values onto parallel arrays
        material_calc(&unity, &x_delta, &materials, &num_r_cells, delta_x, NUM_CELLS, &m_vfrac_0, &m_vfrac_1);

        // Statistically track volume fraction values
        for (c = 0; c < NUM_CELLS; c++) {
            push(&*((*s_vfrac_0) + c), *(m_vfrac_0 + c));
            push(&*((*s_vfrac_1) + c), *(m_vfrac_1 + c));
        }

        if ((i + 1) % NUM_SAY == 0) {
            printf("\rRealization Number: %ld / %ld", i + 1, NUM_REALIZATIONS);
            fflush(stdout);
        }
    }
    CLEANUP();
}

#endif
