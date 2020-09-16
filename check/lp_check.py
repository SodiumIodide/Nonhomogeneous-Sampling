#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import pandas as pd

def flux_plot():
    lp_data = pd.read_csv("uniform_uniform.csv")
    numerical_data = pd.read_csv("../csv/flux_cox_uniform.csv")
    monosample_data = pd.read_csv("../csv/flux_cox_uniform_monosample.csv")
    plt.plot(lp_data['x'], lp_data['flux_0'], label="LP Material 0", color='r')
    plt.plot(lp_data['x'], lp_data['flux_1'], label="LP Material 1", color='b')
    plt.plot(numerical_data['distance'], numerical_data['flux0'], label="Numerical Material 0", color='r', linestyle=':')
    plt.plot(numerical_data['distance'], numerical_data['flux1'], label="Numerical Material 1", color='b', linestyle=':')
    plt.plot(monosample_data['distance'], monosample_data['flux0'], label="Mono Material 0", color='r', linestyle='--')
    plt.plot(monosample_data['distance'], monosample_data['flux1'], label="Mono Material 1", color='b', linestyle='--')
    plt.xlabel("x")
    plt.ylabel("Flux")
    plt.xscale("linear")
    plt.yscale("log")
    plt.xlim([0, 10])
    plt.legend(loc='best')
    plt.grid(which='both', axis='both')
    plt.tight_layout()
    plt.savefig("lp_comparison.png")
    plt.cla()
    plt.clf()

def vf_plot():
    lp_data = pd.read_csv("uniform_uniform.csv")
    numerical_data = pd.read_csv("../csv/flux_cox_uniform.csv")
    monosample_data = pd.read_csv("../csv/flux_cox_uniform_monosample.csv")
    plt.plot(lp_data['x'], lp_data['vf_0'], label="LP Material 0", color='r')
    plt.plot(lp_data['x'], lp_data['vf_1'], label="LP Material 1", color='b')
    plt.plot(numerical_data['distance'], numerical_data['flux0'], label="Numerical Material 0", color='r', linestyle=':')
    plt.plot(numerical_data['distance'], numerical_data['flux1'], label="Numerical Material 1", color='b', linestyle=':')
    plt.plot(monosample_data['distance'], monosample_data['flux1'], label="Mono Material 0", color='r', linestyle='--')
    plt.plot(monosample_data['distance'], monosample_data['flux1'], label="Mono Material 1", color='b', linestyle='--')
    plt.xlabel("x")
    plt.ylabel("Volume Fraction")
    plt.xscale("linear")
    plt.yscale("linear")
    plt.legend(loc='best')
    plt.xlim([0, 10])
    plt.ylim([0, 1])
    plt.tight_layout()
    plt.savefig("lp_comparison_vf.png")
    plt.cla()
    plt.clf()

def main():
    if len(sys.argv) > 1 and "v" in sys.argv[1]:
        vf_plot()
    else:
        flux_plot()

if __name__ == '__main__':
    main()
