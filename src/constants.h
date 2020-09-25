#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#ifndef _CONFIG_H
#include "config.h"
#endif

// Generator seed (initial)
const int SEED = 1234;
// Number of ordinates
const long int NUM_ORDS = 2;
// Chords
const double START_VALUE[2] = {
    101.0 / 20.0,
    101.0 / 20.0,
};  // cm
const double END_VALUE[2] = {
    99.0 / 10.0,
    11.0 / 10.0,
};  // cm
const double COX_START_VALUE[2] = {
    1.0 / 10.0,
    1.0 / 10.0,
};
const double COX_END_VALUE[2] = {
    197.0 / 10.0,
    21.0 / 10.0,
};
const double GAUSSIAN_VARIANCE[2] = {
    1.0,
    0.1,
};
// Distances
const double END_DIST = 1e1;  // cm
const int NUM_DIVS = 100;
// Total XS
const double SIGMA_T[2] = {
    10.0 / 99.0,
    100.0 / 11.0,
};  // cm^-1
// Scattering coefficient
const double SCAT_COEFF[2] = {
    0.0,
    0.0,
};
// Spontaneous volumetric source constant
const double SPONT_SOURCE_CONST[2] = {
    0.0,
    0.0,
};  // cm^-3
// Number of piecewise segments in chord function
#ifdef PIECEWISE
const int NUM_SEGMENTS = 1000;
#endif
// Number of space-steps
const long int NUM_CELLS = (long int)1e3L;
// Initial beamline flux
const double INCIDENT_ANGULAR_FLUX = 1.0;
const double REVERSE_ANGULAR_FLUX = 0.0;
// Initial scalar flux value - Pure Absorber model
#ifdef PUREABS
const double FLUX_INIT = 1.0;
#endif
// Number of cycles to time
#ifdef GEOMETRYTIME
const long int NUM_TIME = 1e5L;
#endif
// Number of realizations
const long int NUM_REALIZATIONS = (long int)1e3L;
// Number of outer realizations for overall convergence during the Cox process
const long int NUM_COX_REALIZATIONS_OUTER = (long int)1e3L;
// Number of realizations (geometry timing)
const long int NUM_TIME = (long int)1e7L;
// Number to display progress
const long int NUM_SAY = (long int)1e3L;
// Tolerance
const double TOLERANCE = 1e-6;

#endif
