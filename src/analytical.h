#ifndef _ANALYTICAL_H
#define _ANALYTICAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

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

#define CLEANUP() do {\
    free(x_delta);\
    free(x_arr);\
    free(materials);\
    free(phi_mat_0);\
    free(phi_mat_1);\
    free(phi_mat_stat_0);\
    free(phi_mat_stat_1);\
    free(cell_vector);\
} while (0)

#define CLEANUP_LOOP() do{\
    free(alpha_beta_vec);\
    free(mult_mat);\
    free(phi_morph);\
    free(pivot);\
} while (0)

int dgesv_(int *n, int *nrhs, double *a,
           int *lda, int *ipiv, double *b, int *ldb, int *info);

void analytical() {
    // Seed built-in generator
    srand(SEED);
    // Iterators
    long int c, i;
    // Material absorption cross-sections
    double sigma_a[2] = {
        (1.0 - SCAT_COEFF[0]) * SIGMA_T[0],
        (1.0 - SCAT_COEFF[1]) * SIGMA_T[1]
    };  // cm^-1
    // Material length scale (not used if sigma_a is zero, should fill with NaN values regardless)
    double length_scale[2] = {
        sqrt(1.0 / (sigma_a[0] * SIGMA_T[0])),
        sqrt(1.0 / (sigma_a[1] * SIGMA_T[1]))
    };  // cm
    // Geometry variables
    // Initially allocate to 1-length
    double *x_delta = malloc(sizeof *x_delta * NUM_DIVS);
    double *x_arr = malloc(sizeof *x_arr * NUM_DIVS);
    int *materials = malloc(sizeof *materials * NUM_DIVS);
    long num_r_cells = 0;
    double delta_x = END_DIST / (double)NUM_CELLS;
    // Solution variables
    // Statistics flux variables - structured
    runningstat *phi_mat_stat_0 = malloc(sizeof *phi_mat_stat_0 * NUM_CELLS);
    runningstat *phi_mat_stat_1 = malloc(sizeof *phi_mat_stat_1 * NUM_CELLS);
    // Structured flux results from map onto array
    double *phi_mat_0 = malloc(sizeof *phi_mat_0 * NUM_CELLS);
    double *phi_mat_1 = malloc(sizeof *phi_mat_1 * NUM_CELLS);
    // Generic return
    int success = 0;

    // File handle
    FILE *fp;
    // Distance counter
    double distance = 0.0;

    // X values of flux for plotting
    double *cell_vector = malloc(sizeof *cell_vector * NUM_CELLS);
    linspace(cell_vector, 0.0, END_DIST, NUM_CELLS);

    // Boundary conditions: 1/cm^3-s-MeV-strad
    // Left boundary
    double psi_bound_l = INCIDENT_ANGULAR_FLUX;
    // Right boundary
    double psi_bound_r = REVERSE_ANGULAR_FLUX;

    // Initialize
    for (c = 0; c < NUM_CELLS; c++) {
        *(phi_mat_0 + c) = 0.0;
        *(phi_mat_1 + c) = 0.0;
        *(phi_mat_stat_0 + c) = init_runningstat();
        *(phi_mat_stat_1 + c) = init_runningstat();
    }

    // Calculation: iterations
    int cont_calc_outer = 1;
    // Start counter at 0
    long int iterations_outer = 0;
    // Matrix filling assists
    int material_prev, material_curr;
    double x_loc;
    // Lapack variables
    int info;
    int rhs = 1;

    // Outer loop over overall mixed system
    while (cont_calc_outer) {
        // Fill Markovian geometry
        GEOMETRY_GEN();
        if (!success) {
            CLEANUP();
            printf("Failed to generate geometry\n");
            exit(EXIT_FAILURE);
        }

        // Calculational arrays
        int arr_size = (num_r_cells / NUM_DIVS) * 2;
        double *alpha_beta_vec = calloc(arr_size, sizeof *alpha_beta_vec);
        double *mult_mat = calloc(arr_size * arr_size, sizeof *mult_mat);
        double *phi_morph = calloc(num_r_cells, sizeof *phi_morph);
        int *pivot = calloc(arr_size, sizeof *pivot);

        // Initialize boundary terms
        *(alpha_beta_vec + 0) = 2.0 * psi_bound_l;
        *(alpha_beta_vec + arr_size - 1) = 2.0 * psi_bound_r;

        // First row
        if (sigma_a[*(materials + 0)] != 0.0) {
            *(mult_mat + 0 * arr_size + 0) = 1.0 - 1.0 / (length_scale[*(materials + 0)] * SIGMA_T[*(materials + 0)]);
            *(mult_mat + 1 * arr_size + 0) = 1.0 + 1.0 / (length_scale[*(materials + 0)] * SIGMA_T[*(materials + 0)]);
        } else {
            *(mult_mat + 0 * arr_size + 0) = -1.0 / SIGMA_T[*(materials + 0)];
            *(mult_mat + 1 * arr_size + 0) = 1.0;
        }

        // Store in column-major format for LAPACK
        i = NUM_DIVS;
        for (c = 1; c < arr_size - 1; c += 2) {
            // Assign filler helpers
            material_prev = *(materials + i - 1);
            material_curr = *(materials + i);
            x_loc = *(x_arr + i);  // Left-cell values
            // First row
            // Previous step
            if (sigma_a[material_prev] != 0.0) {
                *(mult_mat + (c - 1) * arr_size + c) = -exp(x_loc / length_scale[material_prev]);
                *(mult_mat + (c + 0) * arr_size + c) = -exp(-x_loc / length_scale[material_prev]);
            } else {
                *(mult_mat + (c - 1) * arr_size + c) = -x_loc;
                *(mult_mat + (c + 0) * arr_size + c) = -1.0;
            }
            // Current step
            if (sigma_a[material_curr] != 0.0) {
                *(mult_mat + (c + 1) * arr_size + c) = exp(x_loc / length_scale[material_curr]);
                *(mult_mat + (c + 2) * arr_size + c) = exp(-x_loc / length_scale[material_curr]);
            } else {
                *(mult_mat + (c + 1) * arr_size + c) = x_loc;
                *(mult_mat + (c + 2) * arr_size + c) = 1.0;
            }
            // Second row
            // Previous step
            if (sigma_a[material_prev] != 0.0) {
                *(mult_mat + (c - 1) * arr_size + c + 1) = -1.0 / (SIGMA_T[material_prev] * length_scale[material_prev]) * exp(x_loc / length_scale[material_prev]);
                *(mult_mat + (c + 0) * arr_size + c + 1) = 1.0 / (SIGMA_T[material_prev] * length_scale[material_prev]) * exp(-x_loc / length_scale[material_prev]);
            } else {
                *(mult_mat + (c - 1) * arr_size + c + 1) = -1.0 / SIGMA_T[material_prev];
                *(mult_mat + (c + 0) * arr_size + c + 1) = 0.0;
            }
            // Current step
            if (sigma_a[material_curr] != 0.0) {
                *(mult_mat + (c + 1) * arr_size + c + 1) = 1.0 / (SIGMA_T[material_curr] * length_scale[material_curr]) * exp(x_loc / length_scale[material_curr]);
                *(mult_mat + (c + 2) * arr_size + c + 1) = -1.0 / (SIGMA_T[material_curr] * length_scale[material_curr]) * exp(-x_loc / length_scale[material_curr]);
            } else {
                *(mult_mat + (c + 1) * arr_size + c + 1) = 1.0 / SIGMA_T[material_curr];
                *(mult_mat + (c + 2) * arr_size + c + 1) = 0.0;
            }
            // Track interface count
            i += NUM_DIVS;
        }

        // Last row
        if (sigma_a[*(materials + num_r_cells - 1)] != 0.0) {
            *(mult_mat + (arr_size - 1) * arr_size - 1) = (1.0 + 1.0 / (length_scale[*(materials + num_r_cells - 1)] * SIGMA_T[*(materials + num_r_cells - 1)])) * exp((*(x_arr + num_r_cells - 1) + *(x_delta + num_r_cells - 1)) / length_scale[*(materials + num_r_cells - 1)]);
            *(mult_mat + arr_size * arr_size - 1) = (1.0 - 1.0 / (length_scale[*(materials + num_r_cells - 1)] * SIGMA_T[*(materials + num_r_cells - 1)])) * exp(-(*(x_arr + num_r_cells - 1) + *(x_delta + num_r_cells - 1)) / length_scale[*(materials + num_r_cells - 1)]);
        } else {
            *(mult_mat + (arr_size - 1) * arr_size - 1) = *(x_arr + num_r_cells - 1) + *(x_delta + num_r_cells - 1) + 1.0 / SIGMA_T[*(materials + num_r_cells - 1)];
            *(mult_mat + arr_size * arr_size - 1) = 1.0;
        }

        dgesv_(&arr_size, &rhs, mult_mat, &arr_size, pivot, alpha_beta_vec, &arr_size, &info);

        if (info != 0) {
            // Write an error message
            printf("dgesv returned error %d\n", info);
        }

        iterations_outer++;

        // Assign values of flux from calculation of integration constants
        c = 0;
        for (i = 0; i < num_r_cells; i++) {
            if ((i % NUM_DIVS == 0) && (i != 0)) c += 2;
            if (sigma_a[*(materials + i)] != 0.0) {
                *(phi_morph + i) = *(alpha_beta_vec + c) * exp((*(x_arr + i) + *(x_delta + i) / 2.0) / length_scale[*(materials + i)]) + *(alpha_beta_vec + c + 1) * exp(-(*(x_arr + i) + *(x_delta + i) / 2.0) / length_scale[*(materials + i)]);
            } else {
                *(phi_morph + i) = *(alpha_beta_vec + c) * (*(x_arr + i) + *(x_delta + i) / 2.0) + *(alpha_beta_vec + c + 1);
            }
        }

        // Obtain material balances
        // Average flux in structured cells
        material_calc(&phi_morph, &x_delta, &materials, &num_r_cells, delta_x, NUM_CELLS, &phi_mat_0, &phi_mat_1);

        for (c = 0; c < NUM_CELLS; c++) {
            if (*(phi_mat_0 + c) != 0.0) {
                push(&*(phi_mat_stat_0 + c), *(phi_mat_0 + c));
            }
            if (*(phi_mat_1 + c) != 0.0) {
                push(&*(phi_mat_stat_1 + c), *(phi_mat_1 + c));
            }
        }  // Cell loop

        if (iterations_outer % NUM_SAY == 0) {
            printf("\rRealization Number: %ld / %ld", iterations_outer, NUM_REALIZATIONS);
            fflush(stdout);
        }

        // Produce a one-realization profile
        if (iterations_outer == 1) {
            fp = fopen("../csv/one_realization.csv", "w");
            fprintf(fp, "distance,flux\n");
            for (i = 0; i < num_r_cells; i++) {
                distance += *(x_delta + i);
                fprintf(fp, "%f,%f\n", distance - *(x_delta + i) / 2.0, *(phi_morph + i));
            }
            fclose(fp);
            printf("First realization data written\n");
        }

        if (iterations_outer > NUM_REALIZATIONS) cont_calc_outer = 0;

        // Cleanup
        CLEANUP_LOOP();
    }  // Outer loop
    printf("\n");

    // Save flux data
    OPEN_FILE();
    fprintf(fp, "distance,flux0,varflux0,flux1,varflux1\n");
    distance = 0.0;
    for (c = 0; c < NUM_CELLS; c++) {
        distance += delta_x;
        fprintf(fp, "%f,%f,%f,%f,%f\n", distance - delta_x / 2.0, mean(&*(phi_mat_stat_0 + c)), variance(&*(phi_mat_stat_0 + c)), mean(&*(phi_mat_stat_1 + c)), variance(&*(phi_mat_stat_1 + c)));
    }

    printf("Calculation done\n");

    // Cleanup
    fclose(fp);
    CLEANUP();
}

#endif
