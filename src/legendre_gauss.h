#ifndef _LEGENDRE_GAUSS_H
#define _LEGENDRE_GAUSS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef _LINSPACE_H
#include "linspace.h"
#endif

#ifndef PI
#define PI (4.0 * atan(1.0))
#endif

double maximum_subt(double *arr_a, double *arr_b) {
    // Error check
    if (sizeof *arr_a != sizeof *arr_b) {
        printf("ERROR: Array sizes in maximum subtraction are dissimilar!\n");
        return 0.0;
    }

    double value = fabs(*(arr_a) - *(arr_b));
    double test;

    for (unsigned int i = 1; i < (sizeof *arr_a / sizeof *(arr_a)); i++) {
        test = fabs(*(arr_a + i) - *(arr_b + i));
        if (test > value) {
            value = test;
        }
    }
    
    return value;
}

void reverse_arr(double arr[], int start, int end) {
    double temp;
    while (start < end) {
        temp = arr[start];
        arr[start] = arr[end];
        arr[end] = temp;
        start++;
        end--;
    }
}

// Returns a boolean success value
int legendre_gauss_quad(int order, double intv_a, double intv_b, double values[], double weights[]) {
    // Compute Legendre-Gauss values and weights on an interval
    // [intv_a, intv_b] with an input truncation order

    // Test for truncation order
    if (order < 1) {
        printf("Truncation order must be an integer greater than 0\n");
        return 0;
    }

    // Vector variables
    double *x_space = malloc(sizeof *x_space * order);
    double *y_space = malloc(sizeof *y_space * order);
    double *y_hold = malloc(sizeof *y_hold * order);
    double *prime = malloc(sizeof *prime * order);
    // Scalar variables
    int order_place = order - 1;
    int order_one = order_place + 1;
    int order_two = order_place + 2;
    // Iterator
    int i, j;
    // Initial parameter for constant
    const double y_const = 2.0;
    // Legendre-Gauss Vandermonde matrix and its derivative
    double *legendre = malloc(sizeof *legendre * (order) * (order + 1));

    // Assign the arrays
    linspace(x_space, -1.0, 1.0, order_one);
    for (i = 0; i < order; i++) {
        *(y_space + i) = cos((2.0 * i + 1.0) / (2.0 * order_place + 2.0) * PI) + (0.27 / order_one) * sin((PI * *(x_space + i) * order_place) / order_two);
    }
    for (i = 0; i < (order) * (order + 1); i++) {
        *(legendre + i) = 0.0;
    }

    // Compute the zeros of the N+1 Legendre Polynomial using the recursion
    // relation and the Newton-Raphson method
    for (i = 0; i < order; i++) {
        *(y_hold + i) = y_const;
    }
    while (maximum_subt(y_space, y_hold) > DBL_EPSILON) {
        // First Legendre Polynomial
        for (i = 0; i < order; i++) {
            *(legendre + i * (order + 1) + 0) = 1.0;
        }

        // Second Legendre Polynomial
        for (i = 0; i < order; i++) {
            *(legendre + i * (order + 1) + 1) = *(y_space + i);
        }

        for (i = 0; i < order; i++) {
            for (j = 2; j <= order_one; j++) {
                *(legendre + i * (order + 1) + j) = ((2.0 * j - 1.0) * *(y_space + i) * *(legendre + i * (order + 1) + j - 1) - (j - 1.0) * *(legendre + i * (order + 1) + j - 2)) / (double)j;
            }
        }

        // Derivative of the Legendre Polynomial
        for (i = 0; i < order; i++) {
            *(prime + i) = order_two * (*(legendre + i * (order + 1) + order_one - 1) - *(y_space + i) * *(legendre + i * (order + 1) + order_two - 1)) / (1.0 - pow(*(y_space + i), 2));

            *(y_hold + i) = *(y_space + i);
            *(y_space + i) = *(y_hold + i) - *(legendre + i * (order + 1) + order_two - 1) / *(prime + i);
        }
    }

    for (i = 0; i < order; i++) {
        // Linear map from [-1, 1] to [intv_a, intv_b]
        values[i] = intv_a * (1.0 - *(y_space + i)) / 2.0 + intv_b * (1.0 + *(y_space + i)) / 2.0;

        // Compute the weights
        weights[i] = (double)(intv_b - intv_a) / ((1.0 - pow(*(y_space + i), 2)) * pow(*(prime + i), 2)) * pow((double)order_two / (double)order_one, 2);
    }

    reverse_arr(values, 0, order_place);
    reverse_arr(weights, 0, order_place);

    // Free allocations
    free(x_space);
    free(y_space);
    free(y_hold);
    free(prime);
    free(legendre);

    return 1;
}

#endif
