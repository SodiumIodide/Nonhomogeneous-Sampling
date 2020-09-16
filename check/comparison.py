#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

def flux_plot():
    lp_data = pd.read_csv("LP.csv")
    numerical_data = pd.read_csv("../csv/flux_linear.csv")
    plt.plot(lp_data['x'], lp_data['flux_0'], label="LP Material 0", color='r')
    plt.plot(lp_data['x'], lp_data['flux_1'], label="LP Material 1", color='b')
    plt.plot(numerical_data['distance'], numerical_data['flux0'], label="Numerical Material 0", color='r', linestyle=':')
    plt.plot(numerical_data['distance'], numerical_data['flux1'], label="Numerical Material 1", color='b', linestyle=':')
    plt.xlabel("x")
    plt.ylabel("Flux")
    plt.xscale("linear")
    plt.yscale("log")
    plt.xlim(0, 10)
    plt.legend(loc="best")
    plt.grid(which="both", axis="both")
    plt.tight_layout()
    plt.savefig("lp_linear_comparison.png")
    plt.cla()
    plt.clf()

def main():
    flux_plot()

if __name__ == '__main__':
    main()