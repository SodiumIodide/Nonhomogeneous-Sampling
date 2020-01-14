#ifndef _PIECEWISE_CONSTANT_QUAD_CHORD_H
#define _PIECEWISE_CONSTANT_QUAD_CHORD_H

#include "quad_chord.h"

static inline double piecewise_constant_quad_chord(double param_a, double param_b, double param_c, int num_segments, double end_dist, double distance) {
    double delta_dist = end_dist / (double)num_segments;
    int bin_no = (int)(distance / delta_dist);
    if (distance >= end_dist) bin_no = num_segments - 1;
    double seg_init_distance = delta_dist * bin_no;
    double seg_end_distance = delta_dist * (bin_no + 1);
    // Preserve quadratic average
    return (param_a / 3.0) * (pow(seg_end_distance + seg_init_distance, 2) - seg_end_distance * seg_init_distance) + (param_b / 2.0) * (seg_end_distance + seg_init_distance) + param_c;
}

#endif
