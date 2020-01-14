#ifndef _LINEAR_CHORD_H
#define _LINEAR_CHORD_H

static inline double linear_chord(double chord_initial, double chord_slope, double distance) {
    return chord_initial + chord_slope * distance;
}

#endif
