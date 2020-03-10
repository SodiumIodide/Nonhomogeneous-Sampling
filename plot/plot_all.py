#!/usr/bin/env python3

'''
Plot all data found in ../csv/AllData
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

XSCALE = "linear"
YSCALE = "log"
DATAPATH = "../csv/AllData"
PLOTPATH = "../img"
STD = False

class data:
    def __init__(self, in_data):
        self.std_flux_0 = in_data['varflux0'].apply(np.abs).apply(np.sqrt)
        self.std_flux_1 = in_data['varflux1'].apply(np.abs).apply(np.sqrt)
        self.lb_flux_0 = in_data['flux0'] - self.std_flux_0
        self.ub_flux_0 = in_data['flux0'] + self.std_flux_0
        self.lb_flux_1 = in_data['flux1'] - self.std_flux_1
        self.ub_flux_1 = in_data['flux1'] + self.std_flux_1
        self.distance = in_data['distance']
        self.flux_0 = in_data['flux0']
        self.flux_1 = in_data['flux1']

    def plot(self, filename):
        # Plot flux
        plt.plot(self.distance, self.flux_0, color='r', label="Material 0")
        plt.plot(self.distance, self.flux_1, color='b', label="Material 1")
        if STD:
            plt.fill_between(self.distance, self.lb_flux_0, self.ub_flux_0, color='r', alpha=0.5)
            plt.fill_between(self.distance, self.lb_flux_1, self.ub_flux_1, color='b', alpha=0.5)
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

def trial(csv_name):
    try:
        p_data = pd.read_csv(f"{DATAPATH}/{csv_name}.csv")
        obj = data(p_data)
        obj.plot(f"{csv_name}")
    except FileNotFoundError:
        print(f"Could not process {csv_name}")

def main():
    '''Main function'''
    # Each problem has 4 cases
    for p_num in range(4):
       trial(f"flux_linear_{p_num}")
       trial(f"flux_linear_thinning_{p_num}")
       trial(f"flux_quad_thinning_{p_num}")
       trial(f"flux_piecewise_linear_{p_num}")
       trial(f"flux_piecewise_linear_thinning_{p_num}")
       trial(f"flux_piecewise_quad_{p_num}")
       trial(f"flux_piecewise_quad_thinning_{p_num}")

if __name__ == '__main__':
    main()
