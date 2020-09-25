#ifndef _S2_H
#define _S2_H

#include <stdio.h>

#include <gsl/gsl_rng.h>

#ifndef _LEGENDRE_GAUSS_H
#include "legendre_gauss.h"
#endif

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
    free(phi_mat_0);\
    free(phi_mat_1);\
    free(mu);\
    free(weights);\
    free(psi_bound_l);\
    free(psi_bound_r);\
} while (0)

#define CLEANUP_LOOP() do{\
    free(psi);\
    free(psi_i_p);\
    free(psi_i_m);\
    free(scat_source);\
    free(spont_source);\
    free(tot_source);\
    free(phi_morph_new);\
    free(phi_morph_old);\
} while (0)

double max_relative_error(double *arr_a, double *arr_b) {
    if (sizeof *arr_a != sizeof *arr_b) {
        printf("ERROR: Array sizes in relative error are dissimilar!\n");
        return 0.0;
    }

    double value = fabs(*(arr_a) - *(arr_b)) / *(arr_b);
    double test;

    for (unsigned int i = 1; i < (sizeof *arr_a / sizeof *(arr_a)); i++) {
        test = fabs(*(arr_a + i) - *(arr_b + i)) / *(arr_b);
        if (test > value) {
            value = test;
        }
    }

    return value;
}

void s2(const gsl_rng* rng, runningstat** phi_mat_stat_0, runningstat** phi_mat_stat_1, double chord_0, double chord_1) {
    if (chord_0 == 0.0 && chord_1 == 0.0) {
        (void)chord_0;
        (void)chord_1;
    }

    // Iterators
    long int c, m;
    long int num_iter_inner = (long int)1e5;
    // Material scattering cross-sections
    double sigma_s[2] = {
        SCAT_COEFF[0] * SIGMA_T[0],
        SCAT_COEFF[1] * SIGMA_T[1]
    };  // cm^-1
    // Geometry variables
    // Initially allocate to 1-length
    double *x_delta = malloc(sizeof *x_delta * NUM_DIVS);
    double *x_arr = malloc(sizeof *x_arr * NUM_DIVS);
    int *materials = malloc(sizeof *materials * NUM_DIVS);
    long num_r_cells = 0;
    double delta_x = END_DIST / (double)NUM_CELLS;
    // Structured flux results from map onto array
    double *phi_mat_0 = malloc(sizeof *phi_mat_0 * NUM_CELLS);
    double *phi_mat_1 = malloc(sizeof *phi_mat_1 * NUM_CELLS);
    // Generic return
    int success = 0;

    // Ordinates and weights
    double *mu = malloc(sizeof *mu * NUM_ORDS);
    double *weights = malloc(sizeof *weights * NUM_ORDS);

    // Boundary conditions: 1/cm^3-s-MeV-strad
    // Left boundary
    double *psi_bound_l = malloc(sizeof *psi_bound_l * NUM_ORDS);
    // Right boundary
    double *psi_bound_r = malloc(sizeof *psi_bound_r * NUM_ORDS);
    // Assignment
    for (m = 0; m < NUM_ORDS; m++) {
        // LEFT
        // Isotropic source:
        *(psi_bound_l + m) = INCIDENT_ANGULAR_FLUX;
        // Vacuum source:
        //*(psi_bound_l + m) = 0.0;
        // RIGHT
        // Isotropic source:
        *(psi_bound_r + m) = REVERSE_ANGULAR_FLUX;
        // Vacuum source:
        //*(psi_bound_r + m) = 0.0;
    }

    // Initialize
    for (c = 0; c < NUM_CELLS; c++) {
        *(phi_mat_0 + c) = 0.0;
        *(phi_mat_1 + c) = 0.0;
    }

    // Legendre Gauss Quadrature over chosen ordinates
    if (NUM_ORDS > 2) {
        success = legendre_gauss_quad(NUM_ORDS, -1.0, 1.0, mu, weights);
        if (!success) {
            printf("Failed in construction of weighted ordinates\n");
            CLEANUP();
        }
    } else if (NUM_ORDS < 2) {
        printf("Cannot construct less than two ordinates for SN model\n");
        CLEANUP();
    } else {
        *(mu + 0) = -1.0;
        *(mu + 1) = 1.0;
        *(weights + 0) = 1.0;
        *(weights + 1) = 1.0;
    }

    // Calculation: iterations
    int cont_calc_outer = 1;
    // Start counter at 0
    long int iterations_outer = 0;

    // Outer loop over overall mixed system
    while (cont_calc_outer) {
        int cont_calc_inner = 1;
        // Start inner iterations at zero
        long int iterations_inner = 0;

        // Fill Markovian geometry
        GEOMETRY_GEN();
        if (!success) {
            CLEANUP();
            printf("Failed to generate geometry\n");
            exit(EXIT_FAILURE);
        }

        // Holder for cell indices
        int first_cell = 0;
        int last_cell = num_r_cells - 1;

        // Calculational arrays
        double *psi = malloc(sizeof *psi * num_r_cells * NUM_ORDS);
        double *psi_i_p = malloc(sizeof *psi_i_p * num_r_cells * NUM_ORDS);
        double *psi_i_m = malloc(sizeof *psi_i_m * num_r_cells * NUM_ORDS);

        // Initial source terms
        double *scat_source = malloc(sizeof *scat_source * num_r_cells);
        double *spont_source = malloc(sizeof *spont_source * num_r_cells);
        double *tot_source = malloc(sizeof *tot_source * num_r_cells);

        // Phi calculations
        double *phi_morph_new = malloc(sizeof *phi_morph_new * num_r_cells);
        double *phi_morph_old = malloc(sizeof *phi_morph_old * num_r_cells);

        // Initialize
        for (c = 0; c < num_r_cells; c++) {
            *(scat_source + c) = 0.0;
            *(spont_source + c) = 0.0;
            *(tot_source + c) = 0.0;
            *(phi_morph_new + c) = 1.0;
            *(phi_morph_old + c) = 0.0;
            for (m = 0; m < NUM_ORDS; m++) {
                *(psi + c * NUM_ORDS + m) = 0.0;
                *(psi_i_p + c * NUM_ORDS + m) = 0.0;
                *(psi_i_m + c * NUM_ORDS + m) = 0.0;
            }
        }

        // Flag to tell whether to use the inner loop data in averaging or not
        int converged = 0;

        while (cont_calc_inner) {
            // Update flux, set source terms
            for (c = 0; c < num_r_cells; c++) {
                *(phi_morph_old + c) = *(phi_morph_new + c);
                *(scat_source + c) = sigma_s[*(materials + c)] / 2.0 * *(phi_morph_new + c);
                *(spont_source + c) = SPONT_SOURCE_CONST[*(materials + c)] / 2.0;
                *(tot_source + c) = *(scat_source + c) + *(spont_source + c);
            }

            // Forward sweep (left to right)
            // First cell (left boundary)
            // Ordinate loop, only consider the pos. ords for forward motion
            for (m = NUM_ORDS / 2; m < NUM_ORDS; m++) {
                *(psi + first_cell * NUM_ORDS + m) = pow((1.0 + (SIGMA_T[*(materials + first_cell)] * *(x_delta + first_cell)) / (2.0 * fabs(*(mu + m)))), -1) * (*(psi_bound_l + m) + (*(tot_source + first_cell) * *(x_delta + first_cell)) / (2.0 * fabs(*(mu + m))));
                *(psi_i_p + first_cell * NUM_ORDS + m) = 2.0 * *(psi + first_cell * NUM_ORDS + m) - *(psi_bound_l + m);
            }
            // Rest of the cells (sans left bounding cell)
            for (c = 1; c < num_r_cells; c++) {
                for (m = NUM_ORDS / 2; m < NUM_ORDS; m++) {
                    // Continuity of boundaries
                    *(psi_i_m + c * NUM_ORDS + m) = *(psi_i_p + (c - 1) * NUM_ORDS + m);
                    *(psi + c * NUM_ORDS + m) = pow((1.0 + (SIGMA_T[*(materials + c)] * *(x_delta + c)) / (2.0 * fabs(*(mu + m)))), -1) * (*(psi_i_m + c * NUM_ORDS + m) + (*(tot_source + c) * *(x_delta + c)) / (2.0 * fabs(*(mu + m))));
                    *(psi_i_p + c * NUM_ORDS + m) = 2.0 * *(psi + c * NUM_ORDS + m) - *(psi_i_m + c * NUM_ORDS + m);
                }
            }

            // Backward sweep (right to left)
            // First cell (right boundary)
            // Ordinate loop, only consider neg. ords for backwards motion
            for (m = 0; m < NUM_ORDS / 2; m++) {
                *(psi + last_cell * NUM_ORDS + m) = pow((1.0 + (SIGMA_T[*(materials + last_cell)] * *(x_delta + last_cell)) / (2.0 * fabs(*(mu + m)))), -1) * (*(psi_bound_r + m) + (*(tot_source + last_cell) * *(x_delta + last_cell)) / (2.0 * fabs(*(mu + m))));
                *(psi_i_m + last_cell * NUM_ORDS + m) = 2.0 * *(psi + last_cell * NUM_ORDS + m) - *(psi_bound_r + m);
            }
            // Rest of the cells (sans right bounding cell)
            for (c = last_cell - 1; c >= 0; c--) {
                for (m = 0; m < NUM_ORDS / 2; m++) {
                    // Continuity of boundaries
                    *(psi_i_p + c * NUM_ORDS + m) = *(psi_i_m + (c + 1) * NUM_ORDS + m);
                    *(psi + c * NUM_ORDS + m) = pow((1.0 + (SIGMA_T[*(materials + c)] * *(x_delta + c)) / (2.0 * fabs(*(mu + m)))), -1) * (*(psi_i_p + c * NUM_ORDS + m) + (*(tot_source + c) * *(x_delta + c)) / (2.0 * fabs(*(mu + m))));
                    *(psi_i_m + c * NUM_ORDS + m) = 2.0 * *(psi + c * NUM_ORDS + m) - *(psi_i_p + c * NUM_ORDS + m);
                }
            }

            // Calculate phi from psi
            for (c = 0; c < num_r_cells; c++) {
                double weighted_sum = 0.0;
                for (m = 0; m < NUM_ORDS; m++) {
                    weighted_sum += *(weights + m) * *(psi + c * NUM_ORDS + m);
                }
                *(phi_morph_new + c) = weighted_sum;
            }

            // Relative error for inner loop
            iterations_inner += 1;
            if (max_relative_error(phi_morph_old, phi_morph_new) <= TOLERANCE) {
                cont_calc_inner = 0;
                converged = 1;
            } else if (iterations_inner > num_iter_inner) {
                cont_calc_inner = 0;
                printf("\nNo convergence on loop number %ld: quit after maximum %ld iterations; data will not be used\n", iterations_outer, iterations_inner);
            }
        }  // Inner loop

        if (converged) {
            iterations_outer++;

            // Obtain material balances
            // Average flux in structured cells
            material_calc(&phi_morph_new, &x_delta, &materials, &num_r_cells, delta_x, NUM_CELLS, &phi_mat_0, &phi_mat_1);

            for (c = 0; c < NUM_CELLS; c++) {
                if (*(phi_mat_0 + c) != 0.0) {
                    push(&*((*phi_mat_stat_0) + c), *(phi_mat_0 + c));
                }
                if (*(phi_mat_1 + c) != 0.0) {
                    push(&*((*phi_mat_stat_1) + c), *(phi_mat_1 + c));
                }
            } // Cell loop

            if (iterations_outer % NUM_SAY == 0) {
                printf("\rRealization Number: %ld / %ld", iterations_outer, NUM_REALIZATIONS);
                fflush(stdout);
            }

            if (iterations_outer > NUM_REALIZATIONS) cont_calc_outer = 0;
        }  // Logical test for inner convergence: don't average nonconverged samples

        // Cleanup
        CLEANUP_LOOP();
    }  // Outer loop
    CLEANUP();
}

#endif
