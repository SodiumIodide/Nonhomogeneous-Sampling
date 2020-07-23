#ifndef _LINEAR_GEOMETRY_GEN_REJECTION_H
#define _LINEAR_GEOMETRY_GEN_REJECTION_H

#include <math.h>

#include <gsl/gsl_rng.h>

#include "linear_chord.h"

// Returns boolean value for success
int get_geometry_linear_rejection(const gsl_rng* rng, double start_value_0, double end_value_0, double start_value_1, double end_value_1, double end_dist, long num_divs, double** x_delta, double** x_arr, int** materials, long* num_cells) {
    // Computational values
    int material_num;
    double rand_num;
    double chord_start, chord_slope;
    double limiting_value_chord_0, limiting_value_chord_1, limiting_value_chord;
    double prob_accept;
    double buffer_chord;

    // Buffer pointers for realloc
    double *buf_x_delta = NULL;
    double *buf_x_arr = NULL;
    int *buf_materials = NULL;

    // Updateable values for sizing
    double cons_dist = 0.0;  // cm
    double dist = 0.0;  // cm
    long curr_size = num_divs;

    // Problem parameters
    double chord_slope_0 = (end_value_0 - start_value_0) / end_dist;
    double chord_slope_1 = (end_value_1 - start_value_1) / end_dist;

    // Return value and reallocation sizing
    *num_cells = 0;

    // Determine first material to use
    // For linear model, initial distance is at 0.0 and so the probability is equivalent to the constant term ratio
    const double prob_0 = start_value_0 / (start_value_0 + start_value_1);
    material_num = (gsl_rng_uniform_pos(rng) < prob_0) ? 0 : 1;

    // The value that the chord possesses for the maximum value computed (for rejection purposes)
    // As exponential distribution is computed using the inverse of the chord-length, the minimum value is chosen
    // Material "a"
    if (chord_slope_0 > 0.0) {
        limiting_value_chord_0 = start_value_0;
    } else {
        limiting_value_chord_0 = end_value_0;
    }
    // Material "b"
    if (chord_slope_1 > 0.0) {
        limiting_value_chord_1 = start_value_1;
    } else {
        limiting_value_chord_1 = end_value_1;
    }

    while (cons_dist < end_dist) {
        dist = 0.0;  // cm
        // Assign a chord length based on material number
        if (material_num == 0) {
            chord_start = start_value_0;  // cm
            chord_slope = chord_slope_0;  // cm
            limiting_value_chord = limiting_value_chord_0;  // cm
        } else {
            chord_start = start_value_1;  // cm
            chord_slope = chord_slope_1;  // cm
            limiting_value_chord = limiting_value_chord_1;  // cm
        }

        // Loop for rejection sampling
        int accepted = 0;
        while (!accepted) {
            // Generate a random number
            rand_num = gsl_rng_uniform_pos(rng);

            // Sample from a homogeneous distribution of intensity equal to maximum chord
            dist -= limiting_value_chord * log(rand_num);
            if (dist > end_dist) break;

            // Maximum value achieved with minimum length chord
            // Or, conversely, the inverse of the chord (therefore inverse division)
            buffer_chord = linear_chord(chord_start, chord_slope, cons_dist + dist);
            prob_accept = limiting_value_chord / buffer_chord;

            rand_num = gsl_rng_uniform_pos(rng);
            if (rand_num < prob_accept) accepted = 1;
        }

        cons_dist += dist;  // cm

        // Check on thickness to not overshoot the boundary
        if (cons_dist > end_dist) {
            dist += end_dist - cons_dist;  // cm
            cons_dist = end_dist;  // cm
        }

        // Further discretize geometry
        for (long i = 0; i < num_divs; i++) {
            *num_cells += 1L;

            // Expand arrays in increments of the number of divisions in each sublayer
            if (*num_cells > curr_size) {
                curr_size += num_divs;
                buf_x_delta = realloc(*x_delta, sizeof *x_delta * curr_size);
                buf_x_arr = realloc(*x_arr, sizeof *x_arr * curr_size);
                buf_materials = realloc(*materials, sizeof *materials * curr_size);
                if (buf_x_delta == NULL || buf_x_arr == NULL || buf_materials == NULL) {
                    return 0;
                }
                *x_delta = buf_x_delta;
                *x_arr = buf_x_arr;
                *materials = buf_materials;
                buf_x_delta = NULL;
                buf_x_arr = NULL;
                buf_materials = NULL;
            }

            // The width of each cell
            if (*x_delta != NULL) {
                *((*x_delta) + *num_cells - 1) = dist / (double)num_divs;
            }

            // The initial x-location of each cell
            if (*x_arr != NULL) {
                *((*x_arr) + *num_cells - 1) = cons_dist - dist + (dist / (double)num_divs * (double)i);
            }

            // The material present in each cell
            if (*materials != NULL) {
                *((*materials) + *num_cells - 1) = material_num;
            }
        }

        // Update material number
        material_num = (material_num == 0) ? 1 : 0;
    }

    return 1;
}

#endif
