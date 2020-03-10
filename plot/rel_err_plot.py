#!/usr/bin/env python3

'''
For LP problems, relative error with respect to actual data

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

from plot import DATAPATH, PLOTPATH

LPPATH = f"{DATAPATH}/LPData"
ALLPATH = f"{DATAPATH}/AllData"

XSCALE = "linear"
YSCALE = "linear"

class rel_err_data:
    def __init__(self, in_data_0, in_data_1, r_data):
        self.distance = in_data_0['x']
        self.flux_0 = in_data_0['phi1(x)']
        self.flux_1 = in_data_1['phi2(x)']
        self.r_distance = r_data['distance']
        self.r_flux_0 = r_data['flux0']
        self.r_flux_1 = r_data['flux1']
        self._get_err()

    def _get_err(self):
        interp_flux_0 = np.interp(self.distance, self.r_distance, self.r_flux_0)
        interp_flux_1 = np.interp(self.distance, self.r_distance, self.r_flux_1)
        self.rel_err_0 = np.abs(interp_flux_0 - self.flux_0) / interp_flux_0 * 100
        self.rel_err_1 = np.abs(interp_flux_1 - self.flux_1) / interp_flux_1 * 100

    def plot(self, filename):
        # Plot relative error
        plt.plot(self.distance, self.rel_err_0, color='r', label="Material 0")
        plt.plot(self.distance, self.rel_err_1, color='b', label="Material 1")
        plt.xlabel("Distance")
        plt.ylabel("Relative Error (%)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.xlim([self.distance.iat[0], self.distance.iat[-1]])
        if YSCALE == "linear":
            #plt.ylim(bottom=0.0)
            plt.ylim([0.0, 100.0])
        plt.grid(b=True, which='both', axis='both')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/{filename}.png")
        plt.cla()
        plt.clf()

def trial(pname, pnum, porder):
    try:
        data_0 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob{pnum}_MCLOrder{porder}_1.csv")
        data_1 = pd.read_csv(f"{LPPATH}/LPModelRod_Prob{pnum}_MCLOrder{porder}_2.csv")
        data_r = pd.read_csv(f"{ALLPATH}/flux_{pname}_{pnum}.csv")
        obj = rel_err_data(data_0, data_1, data_r)
        obj.plot(f"{pname}_{pnum}_rel_err")
    except:
        raise

def main():
    '''Main function'''
    for pnum in range(4):
        # Linear problems
        trial("linear_thinning", pnum, 1)
        # Quadratic problems
        trial("quad_thinning", pnum, 2)

if __name__ == '__main__':
    main()
