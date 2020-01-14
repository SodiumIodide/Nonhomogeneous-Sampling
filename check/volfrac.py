#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def calc_c(start_value):
    return start_value

def calc_b(delta_m, delta_t, end_dist):
    return (2.0 * delta_m) / end_dist * (1.0 + np.sqrt(1.0 - delta_t / delta_m))

def calc_a(delta_t, param_b, end_dist):
    return (delta_t - param_b * end_dist) / end_dist**2

def piecewise_constant_quad_chord(param_a, param_b, param_c, num_segments, end_dist, distance):
    delta_dist = end_dist / num_segments
    bin_no = distance // delta_dist
    if (distance >= end_dist):
        bin_no = num_segments - 1
    seg_init_distance = delta_dist * bin_no
    seg_end_distance = delta_dist * (bin_no + 1)
    # Preserve quadratic average
    return (param_a / 3.0) * ((seg_end_distance + seg_init_distance)**2 - seg_end_distance * seg_init_distance) + (param_b / 2.0) * (seg_end_distance + seg_init_distance) + param_c

def quad_chord(param_a, param_b, param_c, distance):
    return param_a * distance**2 + param_b * distance + param_c

def linear_chord(chord_initial, chord_slope, distance):
    return chord_initial + chord_slope * distance

def calc_slope(start_value, end_value, end_dist):
    return (end_value - start_value) / end_dist

def main():
    num_cells = int(1e5)
    end_dist = 10.0
    num_segments = int(1e3)
    delta_dist = end_dist / num_cells
    start_value = np.array([101.0 / 20.0, 101.0 / 20.0])
    end_value = np.array([99.0 / 10.0, 11.0 / 10.0])
    mid_value = np.zeros(2)
    # Material 0 local extrema
    if (start_value[0] > end_value[0]):
        mid_value[0] = 2.0 * start_value[0]
    else:
        mid_value[0] = 2.0 * end_value[0]
    # Material 1 local extrema
    if (start_value[1] > end_value[1]):
        mid_value[1] = end_value[1] / 2.0
    else:
        mid_value[1] = start_value[1] / 2.0
    delta_m = mid_value - start_value
    delta_t = end_value - start_value
    # Quadratic parameters
    param_c = np.zeros(2)
    param_b = np.zeros(2)
    param_a = np.zeros(2)
    # Linear parameters
    slope = np.zeros(2)
    # Material 0 parameters
    param_c[0] = calc_c(start_value[0])
    param_b[0] = calc_b(delta_m[0], delta_t[0], end_dist)
    param_a[0] = calc_a(delta_t[0], param_b[0], end_dist)
    slope[0] = calc_slope(start_value[0], end_value[0], end_dist)
    # Material 1 parameters
    param_c[1] = calc_c(start_value[1])
    param_b[1] = calc_b(delta_m[1], delta_t[1], end_dist)
    param_a[1] = calc_a(delta_t[1], param_b[1], end_dist)
    slope[1] = calc_slope(start_value[1], end_value[1], end_dist)
    # Compute analytic volume fractions
    vol_frac_0_quad = np.zeros(int(num_cells))
    vol_frac_1_quad = np.zeros(int(num_cells))
    vol_frac_0_lin = np.zeros(int(num_cells))
    vol_frac_1_lin = np.zeros(int(num_cells))
    prev_value_quad = start_value[0] / (start_value[0] + start_value[1])
    prev_value_lin = prev_value_quad
    distance = 0.0
    for i in range(int(num_cells)):
        distance += delta_dist
        cor_len_quad = (1.0 / quad_chord(param_a[0], param_b[0], param_c[0], distance) + 1.0 / quad_chord(param_a[1], param_b[1], param_c[1], distance))**(-1)
        cor_len_lin = (1.0 / linear_chord(start_value[0], slope[0], distance) + 1.0 / linear_chord(start_value[1], slope[1], distance))**(-1)
        vol_frac_1_quad[i] = prev_value_quad * np.exp(-delta_dist / cor_len_quad) + cor_len_quad / quad_chord(param_a[0], param_b[0], param_c[0], distance) * (1.0 - np.exp(-delta_dist / cor_len_quad))
        vol_frac_0_quad[i] = 1.0 - vol_frac_1_quad[i]
        vol_frac_1_lin[i] = prev_value_lin * np.exp(-delta_dist / cor_len_lin) + cor_len_lin / linear_chord(start_value[0], slope[0], distance) * (1.0 - np.exp(-delta_dist / cor_len_lin))
        vol_frac_0_lin[i] = 1.0 - vol_frac_1_lin[i]
        prev_value_quad = vol_frac_1_quad[i]
        prev_value_lin = vol_frac_1_lin[i]
    # Plot results
    plot_distance = [i * delta_dist for i in range(int(num_cells))]
    # Quadratic
    plt.plot(plot_distance, vol_frac_0_quad, color='r', label="Material 0")
    plt.plot(plot_distance, vol_frac_1_quad, color='b', label="Material 1")
    plt.grid(b=True, which='both', axis='both')
    plt.title("Analytic Quadratic Volume Fraction")
    plt.xlabel("Distance")
    plt.ylabel("Volume Fraction")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("vf_quad.png")
    plt.cla()
    plt.clf()
    # Linear
    plt.plot(plot_distance, vol_frac_0_lin, color='r', label="Material 0")
    plt.plot(plot_distance, vol_frac_1_lin, color='b', label="Material 1")
    plt.grid(b=True, which='both', axis='both')
    plt.title("Analytic Linear Volume Fraction")
    plt.xlabel("Distance")
    plt.ylabel("Volume Fraction")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("vf_lin.png")
    plt.cla()
    plt.clf()

if __name__ == '__main__':
    main()
