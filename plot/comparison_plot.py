#!/usr/bin/env python3

'''
For the thinning versus regular problems

Problem Sets -> Scattering Ratio 0, Scattering Ratio 1:
0 -> 0.0, 0.0
1 -> 0.9, 0.9
2 -> 1.0, 0.0
3 -> 0.0, 1.0
'''

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot import XSCALE, YSCALE, DATAPATH, PLOTPATH

THINNINGPATH = f"{DATAPATH}/AllData"

class thinning_data:
    def __init__(self, t_data, r_data):
        self.distance = r_data['distance']
        self.r_flux_0 = r_data['flux0']
        self.r_flux_1 = r_data['flux1']
        self.t_flux_0 = t_data['flux0']
        self.t_flux_1 = t_data['flux1']
    def plot(self, plotname, filename):
        # Plot flux
        plt.plot(self.distance, self.t_flux_0, ls=':', lw='4', color='r', label="Thinning Material 0")
        plt.plot(self.distance, self.t_flux_1, ls=':', lw='4', color='b', label="Thinning Material 1")
        plt.plot(self.distance, self.r_flux_0, color='r', label="Numerical Material 0")
        plt.plot(self.distance, self.r_flux_1, color='b', label="Numerical Material 1")
        #plt.title(plotname)
        plt.xlabel("Distance")
        plt.ylabel("Flux")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.xlim([self.distance.iat[0], self.distance.iat[-1]])
        plt.grid(b=True, which='both', axis='both')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/{filename}.png")
        plt.cla()
        plt.clf()

def main():
    '''Main function'''
    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_linear_thinning_0.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_linear_0.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "linear_0")
    except:
        print("Linear Problem 0 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_linear_thinning_1.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_linear_1.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "linear_1")
    except:
        print("Linear Problem 1 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_linear_thinning_2.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_linear_2.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "linear_2")
    except:
        print("Linear Problem 2 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_linear_thinning_2.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_linear_2.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "linear_3")
    except:
        print("Linear Problem 3 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_linear_0.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_linear_thinning_0.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "piecewise_linear_0")
    except:
        print("Piecewise Linear Problem 0 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_linear_1.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_linear_thinning_1.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "piecewise_linear_1")
    except:
        print("Piecewise Linear Problem 1 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_linear_2.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_linear_thinning_2.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "piecewise_linear_2")
    except:
        print("Piecewise Linear Problem 2 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_linear_3.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_linear_thinning_3.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "piecewise_linear_3")
    except:
        print("Piecewise Linear Problem 3 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_quad_0.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_quad_thinning_0.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "piecewise_quad_0")
    except:
        print("Piecewise Quadratic Problem 0 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_quad_1.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_quad_thinning_1.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "piecewise_quad_1")
    except:
        print("Piecewise Quadratic Problem 1 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_quad_2.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_quad_thinning_2.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "piecewise_quad_2")
    except:
        print("Piecewise Quadratic Problem 2 not plotted")

    try:
        data_t = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_quad_3.csv")
        data_r = pd.read_csv(f"{THINNINGPATH}/flux_piecewise_quad_thinning_3.csv")
        obj = thinning_data(data_t, data_r)
        obj.plot("Flux Comparison", "piecewise_quad_3")
    except:
        print("Piecewise Quadratic Problem 3 not plotted")

if __name__ == '__main__':
    main()
