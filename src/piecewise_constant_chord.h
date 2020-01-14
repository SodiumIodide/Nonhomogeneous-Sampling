#ifndef _PIECEWISE_CONSTANT_CHORD_H
#define _PIECEWISE_CONSTANT_CHORD_H

#include "linear_chord.h"

static inline double piecewise_constant_chord(double start_value, double end_value, int num_segments, double end_dist, double distance) {
    double slope = (end_value - start_value) / end_dist;
    double delta_dist = end_dist / (double)num_segments;
    int bin_no = (int)(distance / delta_dist);
    if (distance >= end_dist) bin_no = num_segments - 1;
    double init_val = linear_chord(start_value, slope, delta_dist * bin_no);
    double end_val = linear_chord(start_value, slope, delta_dist * (bin_no + 1));
    return (init_val + end_val) / 2.0;
}

#endif
