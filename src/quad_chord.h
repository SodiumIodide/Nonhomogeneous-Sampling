#ifndef _QUAD_CHORD_H
#define _QUAD_CHORD_H

#include "math.h"

static inline double quad_chord(double param_a, double param_b, double param_c, double distance) {
    return param_a * pow(distance, 2) + param_b * distance + param_c;
}

// Midpoint location is 2-fold largest endpoint
static inline double mid_value_0_func(double start_value_0, double end_value_0) {
    return (start_value_0 > end_value_0) ? 2.0 * start_value_0 : 2.0 * end_value_0;
}

// Midpoint location is 0.5-fold smallest endpoint
static inline double mid_value_1_func(double start_value_1, double end_value_1) {
    return (start_value_1 > end_value_1) ? end_value_1 / 2.0 : start_value_1 / 2.0;
}

static inline double calc_c(double start_value) {
    return start_value;
}

static inline double calc_b(double delta_m, double delta_t, double end_dist) {
    return (2.0 * delta_m) / end_dist * (1.0 + sqrt(1.0 - delta_t / delta_m));
}

static inline double calc_a(double delta_t, double param_b, double end_dist) {
    return (delta_t - param_b * end_dist) / pow(end_dist, 2);
}

#endif
