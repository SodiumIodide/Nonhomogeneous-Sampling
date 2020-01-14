#ifndef _PURE_ABS_H
#define _PURE_ABS_H

#include <stdio.h>
#include <stdlib.h>

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
    free(flux);\
    free(s_flux_0);\
    free(s_flux_1);\
    free(m_flux_0);\
    free(m_flux_1);\
} while (0)

void pure_abs() {
    // Seed built-in generator
    srand(SEED);
    // Material absorption cross-sections
    double sigma_a[2] = {
        (1.0 - SCAT_COEFF[0]) * SIGMA_T[0],
        (1.0 - SCAT_COEFF[1]) * SIGMA_T[1]
    };  // cm^-1
    // Geometry variables
    // Iterator
    long int i, j;
    // Initially allocate to 1-length
    double *x_delta = malloc(sizeof *x_delta * NUM_DIVS);
    double *x_arr = malloc(sizeof *x_arr * NUM_DIVS);
    int *materials = malloc(sizeof *materials * NUM_DIVS);
    long num_r_cells = 0;
    double delta_x = END_DIST / (double)NUM_CELLS;
    // Solution variables
    double *flux = malloc(sizeof *flux * NUM_DIVS);
    // Statistics variables - structured
    runningstat *s_flux_0 = malloc(sizeof *s_flux_0 * NUM_CELLS);
    runningstat *s_flux_1 = malloc(sizeof *s_flux_1 * NUM_CELLS);
    for (i = 0; i < NUM_CELLS; i++) {
        s_flux_0[i] = init_runningstat();
        s_flux_1[i] = init_runningstat();
    }
    // Structured results from map onto array
    double *m_flux_0 = malloc(sizeof *m_flux_0 * NUM_CELLS);
    double *m_flux_1 = malloc(sizeof *m_flux_1 * NUM_CELLS);
    // Buffer variable for reallocation
    double *buf_flux = NULL;
    // Generic return
    int success = 0;

    // File handle
    FILE *fp;
    // Distance counter
    double distance = 0.0;

    for (i = 0; i < NUM_REALIZATIONS; i++) {
        // Fill Markovian geometry
        GEOMETRY_GEN();
        if (!success) {
            CLEANUP();
            printf("Failed to generate geometry\n");
            exit(EXIT_FAILURE);
        }

        buf_flux = realloc(flux, sizeof *flux * num_r_cells);
        if (buf_flux == NULL) {
            CLEANUP();
            printf("Error in reallocation: flux!\n");
            exit(EXIT_FAILURE);
        }
        flux = buf_flux;

        // Initialize flux
        *(flux) = FLUX_INIT;

        // Exponential solution to pure absorber problem
        for (j = 1; j < num_r_cells; j++) {
            *(flux + j) = *(flux + j - 1) * exp(-sigma_a[*(materials + j)] * *(x_delta + j - 1));
        }

        // Map flux values onto parallel arrays
        material_calc(&flux, &x_delta, &materials, &num_r_cells, delta_x, NUM_CELLS, &m_flux_0, &m_flux_1);

        // Statistically track flux values
        for (j = 0; j < NUM_CELLS; j++) {
            if (*(m_flux_0 + j) != 0.0) {
                push(&*(s_flux_0 + j), *(m_flux_0 + j));
            }
            if (*(m_flux_1 + j) != 0.0) {
                push(&*(s_flux_1 + j), *(m_flux_1 + j));
            }
        }

        if ((i + 1) % NUM_SAY == 0) {
            printf("\rRealization Number: %ld / %ld", i + 1, NUM_REALIZATIONS);
            fflush(stdout);
        }
    }
    printf("\n");

    // Save flux data
    OPEN_FILE();
    fprintf(fp, "distance,flux0,varflux0,flux1,varflux1\n");
    for (i = 0; i < NUM_CELLS; i++) {
        distance += delta_x;
        fprintf(fp, "%f,%f,%f,%f,%f\n", distance - delta_x / 2.0, mean(&*(s_flux_0 + i)), variance(&*(s_flux_0 + i)), mean(&*(s_flux_1 + i)), variance(&*(s_flux_1 + i)));
    }

    printf("Calculation done\n");

    // Cleanup
    fclose(fp);
    CLEANUP();
}

#endif
