#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

XSCALE = "linear"
YSCALE = "log"
DATAPATH = "../csv"
PLOTPATH = "../img"
XMAX = 1e1
XMIN = 0
YMAX = 1e1
YMIN = 1e-2

class data:
    def __init__(self, in_data):
        self.distance = in_data['distance']
        self.flux = in_data['flux']

    def plot(self, plotname, filename):
        # Plot flux
        plt.plot(self.distance, self.flux, color='k')
        plt.title(plotname)
        plt.xlabel("Distance (cm)")
        plt.ylabel("Scalar Flux")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which='both', axis='both')
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/{filename}.png")
        plt.cla()
        plt.clf()

def main():
    '''Main function'''
    exists = True
    try:
        r_data = pd.read_csv(f"{DATAPATH}/one_realization.csv")
    except FileNotFoundError:
        print("One realization data not found")
        exists = False
    if exists:
        obj = data(r_data)
        obj.plot("Single Realization", "one_realization")

if __name__ == '__main__':
    main()
