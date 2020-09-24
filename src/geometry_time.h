#ifndef _GEOMETRY_TIME_H
#define _GEOMETRY_TIME_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

#include <gsl/gsl_rng.h>

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
} while (0)

static inline unsigned long rdtsc() {
#if defined (__ia64)  // Intel Itanium 64 bit
    unsigned long value;
    __asm__ __volatile__("mov &#3\;0=ar.itc" : "=r"(value) :: "memory");
    while (__builtin_expect ((int)value == -1, 0)) {
        __asm__ __volatile__("mov %0=ar.itc" : "=r"(value) :: "memory");
    }
    return value
#else  // Everything else
    unsigned long long value;
    unsigned int low, high;
    __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));
    value = high;
    value = (value << 32) + low;
    return value;
#endif
}

void geometry_time(const gsl_rng* rng) {
    // Initially allocate to 1-length
    double *x_delta = malloc(sizeof *x_delta * NUM_DIVS);
    double *x_arr = malloc(sizeof *x_arr * NUM_DIVS);
    int *materials = malloc(sizeof *materials * NUM_DIVS);
    long num_r_cells = 0;
    // Generic return
    int success;
    // Time values
    uint64_t cycles;
    // Iterator
    unsigned volatile int i;
    /*
    Note here that "volatile" increases overhead due to forcing the counter to
    be memory-accessed, but as both methods utilize the same amount of loops, the
    overhead introduced should be equivalent, allowing for a comparison.
    */

    // Time generation
    cycles = rdtsc();
    for (i = 0; i < NUM_TIME; i++) {
        // Fill Markovian geometry
        GEOMETRY_GEN();
        if (!success) {
            CLEANUP();
            printf("Failed geometry generation\n");
            exit(EXIT_FAILURE);
        }
        // Cycle-less no-op to prevent optimizing loops
        __asm__ __volatile__("");
    }
    cycles = rdtsc() - cycles;

    printf("Clock cycles for %ld iterations of generation: %ld\n", NUM_TIME, cycles);
    printf("Approximate cycles per function call: %ld\n\n", (long)(ceil((double)cycles / NUM_TIME)));

    // Cleanup
    CLEANUP();
}

#endif
