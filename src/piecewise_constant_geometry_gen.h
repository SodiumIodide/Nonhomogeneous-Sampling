#ifndef _PIECEWISE_CONSTANT_GEOMETRY_GEN_H
#define _PIECEWISE_CONSTANT_GEOMETRY_GEN_H

#include <math.h>

#include <gsl/gsl_rng.h>

#include "piecewise_constant_chord.h"

// Returns boolean value for success
int get_geometry_piecewise_constant(const gsl_rng* rng, double start_value_0, double end_value_0, double start_value_1, double end_value_1, int num_segments, double end_dist, long num_divs, double** x_delta, double** x_arr, int** materials, long* num_cells) {
    // Computational values
    int material_num;
    double rand_num;
    double chord_start, chord_end;

    // Iteration values
    int accepted, segment, current_segment;
    double first_portion, middle_portion;

    // Buffer pointers for realloc
    double *buf_x_delta = NULL;
    double *buf_x_arr = NULL;
    int *buf_materials = NULL;

    // Updateable values for sizing
    double cons_dist = 0.0;  // cm
    double dist = 0.0;  // cm
    long curr_size = num_divs;

    // Return value and reallocation sizing
    *num_cells = 0;

    // Determine first material to use
    // For piecewise constant model, initial distance is at 0.0 and so the probability is equivalent to the term ratio of the averages of the first segments
    double first_value_0 = piecewise_constant_chord(start_value_0, end_value_0, num_segments, end_dist, 0.0);
    double first_value_1 = piecewise_constant_chord(start_value_1, end_value_1, num_segments, end_dist, 0.0);
    const double prob_0 = first_value_0 / (first_value_0 + first_value_1);
    material_num = (gsl_rng_uniform_pos(rng) < prob_0) ? 0 : 1;

    // Delta distance in piecewise function (required to obtain chords at each segment)
    double delta_dist = end_dist / (double)num_segments;

    while (cons_dist < end_dist) {
        // Generate a random number
        rand_num = gsl_rng_uniform_pos(rng);

        // Assign a chord length based on material number
        if (material_num == 0) {
            chord_start = start_value_0;  // cm
            chord_end = end_value_0;  // cm
        } else {
            chord_start = start_value_1;  // cm
            chord_end = end_value_1;  // cm
        }

        // Loop for sampling distance (required to test each segment of piecewise function)
        // Note that segments are zero-indexed
        accepted = 0;
        // Start at current segment
        current_segment = (int)(cons_dist / delta_dist);
        segment = current_segment;
        while (!accepted) {
            if (segment > current_segment) {
                first_portion = ((current_segment + 1) * delta_dist - cons_dist) / piecewise_constant_chord(chord_start, chord_end, num_segments, end_dist, cons_dist);
                middle_portion = 0.0;
                for (int k = current_segment + 1; k < segment; k++) {
                    middle_portion += delta_dist / piecewise_constant_chord(chord_start, chord_end, num_segments, end_dist, ((double)k + 0.5) * delta_dist);
                }
                dist = segment * delta_dist - piecewise_constant_chord(chord_start, chord_end, num_segments, end_dist, ((double)(segment) + 0.5) * delta_dist) * (log(rand_num) + middle_portion + first_portion) - cons_dist;
            } else {
                dist = -log(rand_num) * piecewise_constant_chord(chord_start, chord_end, num_segments, end_dist, cons_dist);
            }
            if (dist > (end_dist - cons_dist)) dist = end_dist - cons_dist;
            if ((dist + cons_dist > (segment + 1) * delta_dist) || (dist + cons_dist < segment * delta_dist)) {
                segment++;
            } else {
                accepted = 1;
                if (dist < 0.0) {
                    accepted = 0;
                    dist = 0.0;
                    rand_num = gsl_rng_uniform_pos(rng);
                }
            }
        }

        cons_dist += dist;  // cm

        // Check on thickness to not overshoot the boundary
        if (cons_dist > end_dist) {
            dist += end_dist - cons_dist;  // cm
            cons_dist = end_dist;  // cm
        }

        // Further discretize geometry
        for (long i = 0; i < num_divs; i++) {
            *num_cells += 1;

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

    free(buf_x_arr);
    free(buf_x_delta);
    free(buf_materials);

    return 1;
}

#endif
