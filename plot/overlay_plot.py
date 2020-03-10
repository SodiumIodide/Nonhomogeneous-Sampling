#!/usr/bin/env python3

'''
For LP problems

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

LPPATH = f"{DATAPATH}/LPData"
ALLPATH = f"{DATAPATH}/AllData"

class lp_data:
    def __init__(self, in_data_0, in_data_1, r_data):
        self.distance = in_data_0['x']
        self.flux_0 = in_data_0['phi1(x)']
        self.flux_1 = in_data_1['phi2(x)']
        self.r_distance = r_data['distance']
        self.r_flux_0 = r_data['flux0']
        self.r_flux_1 = r_data['flux1']

    def plot(self, plotname, filename):
        # Plot flux
        plt.plot(self.distance, self.flux_0, color='r', label="LP Material 0")
        plt.plot(self.distance, self.flux_1, color='b', label="LP Material 1")
        plt.plot(self.r_distance, self.r_flux_0, ls=':', lw='4', color='r', label="Material 0")
        plt.plot(self.r_distance, self.r_flux_1, ls=':', lw='4', color='b', label="Material 1")
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
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob0_MCLOrder1_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob0_MCLOrder1_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_linear_thinning_0.csv")
        obj = lp_data(data_0, data_1, data_r)
        obj.plot("Flux Comparison", "lp_linear_0")
    except:
        print("Linear Problem 0 not plotted")

    try:
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob1_MCLOrder1_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob1_MCLOrder1_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_linear_thinning_1.csv")
        obj = lp_data(data_0, data_1, data_r)
        obj.plot("Flux Comparison", "lp_linear_1")
    except:
        print("Linear Problem 1 not plotted")

    try:
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob2_MCLOrder1_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob2_MCLOrder1_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_linear_thinning_2.csv")
        obj = lp_data(data_0, data_1, data_r)
        obj.plot("Flux Comparison", "lp_linear_2")
    except:
        print("Linear Problem 2 not plotted")

    try:
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob3_MCLOrder1_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob3_MCLOrder1_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_linear_thinning_3.csv")
        obj = lp_data(data_0, data_1, data_r)
        obj.plot("Flux Comparison", "lp_linear_3")
    except:
        print("Linear Problem 3 not plotted")

    try:
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob0_MCLOrder2_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob0_MCLOrder2_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_quad_thinning_0.csv")
        obj = lp_data(data_0, data_1, data_r)
        obj.plot("Flux Comparison", "lp_quad_0")
    except:
        print("Quadratic Problem 0 not plotted")

    try:
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob1_MCLOrder2_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob1_MCLOrder2_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_quad_thinning_1.csv")
        obj = lp_data(data_0, data_1, data_r)
        obj.plot("Flux Comparison", "lp_quad_1")
    except:
        print("Quadratic Problem 1 not plotted")

    try:
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob2_MCLOrder2_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob2_MCLOrder2_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_quad_thinning_2.csv")
        obj = lp_data(data_0, data_1, data_r)
        obj.plot("Flux Comparison", "lp_quad_2")
    except:
        print("Quadratic Problem 2 not plotted")

    try:
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob3_MCLOrder2_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob3_MCLOrder2_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_quad_thinning_3.csv")
        obj = lp_data(data_0, data_1, data_r)
        obj.plot("Flux Comparison", "lp_quad_3")
    except:
        print("Quadratic Problem 3 not plotted")

if __name__ == '__main__':
    main()
