#!/usr/bin/env python3

'''
Use -v as an argument for volume fraction profiles.
'''

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

XSCALE = "linear"
#YSCALE = "log"
YSCALE = "linear"
DATAPATH = "../csv"
PLOTPATH = "../img"
#YMAX = 1e1
#YMIN = 1e-3
YMAX = 1.0
YMIN = 0.0

class data:
    def __init__(self, in_data, ylabel):
        # Standard deviation computing
        self.std_flux_0 = in_data['varflux0'].apply(np.abs).apply(np.sqrt)
        self.std_flux_1 = in_data['varflux1'].apply(np.abs).apply(np.sqrt)
        self.lb_flux_0 = in_data['flux0'] - self.std_flux_0
        self.ub_flux_0 = in_data['flux0'] + self.std_flux_0
        self.lb_flux_1 = in_data['flux1'] - self.std_flux_1
        self.ub_flux_1 = in_data['flux1'] + self.std_flux_1
        self.distance = in_data['distance']
        self.flux_0 = in_data['flux0']
        self.flux_1 = in_data['flux1']
        self.ylabel = ylabel

    def plot(self, plotname, filename):
        # Plot flux
        plt.plot(self.distance, self.flux_0, color='r', label="Material 0")
        plt.plot(self.distance, self.flux_1, color='b', label="Material 1")
        #plt.fill_between(self.distance, self.lb_flux_0, self.ub_flux_0, color='r', alpha=0.5)
        #plt.fill_between(self.distance, self.lb_flux_1, self.ub_flux_1, color='b', alpha=0.5)
        #plt.title(plotname)
        plt.xlabel("Distance (cm)")
        plt.ylabel(self.ylabel)
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.ylim([YMIN, YMAX])
        plt.xlim([self.distance.iat[0], self.distance.iat[-1]])
        plt.grid(b=True, which='both', axis='both')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/{filename}.png")
        plt.cla()
        plt.clf()

def main():
    '''Main function'''
    flux = True
    if len(sys.argv) > 1 and "v" in sys.argv[1]:
        flux = False
        ylabel = "Volume Fraction"
    else:
        ylabel = "Flux"
    r_exists = l_exists = q_r_exists = pw_exists = pw_r_exists = pw_q_exists = pw_q_r_exists = c_u_exists = c_e_exists = c_g_exists = c_u_m_exists = c_e_m_exists = c_g_m_exists = con_exists = True
    try:
        r_data = pd.read_csv(f"{DATAPATH}/flux_linear_thinning.csv")
    except FileNotFoundError:
        print("Linear thinning not found")
        r_exists = False
    try:
        l_data = pd.read_csv(f"{DATAPATH}/flux_linear.csv")
    except FileNotFoundError:
        print("Linear not found")
        l_exists = False
    try:
        q_r_data = pd.read_csv(f"{DATAPATH}/flux_quad_thinning.csv")
    except FileNotFoundError:
        print("Quadratic thinning not found")
        q_r_exists = False
    try:
        pw_data = pd.read_csv(f"{DATAPATH}/flux_piecewise_linear.csv")
    except FileNotFoundError:
        print("Piecewise linear not found")
        pw_exists = False
    try:
        pw_r_data = pd.read_csv(f"{DATAPATH}/flux_piecewise_linear_thinning.csv")
    except FileNotFoundError:
        print("Piecewise linear thinning not found")
        pw_r_exists = False
    try:
        pw_q_data = pd.read_csv(f"{DATAPATH}/flux_piecewise_quad.csv")
    except FileNotFoundError:
        print("Piecewise quadratic not found")
        pw_q_exists = False
    try:
        pw_q_r_data = pd.read_csv(f"{DATAPATH}/flux_piecewise_quad_thinning.csv")
    except FileNotFoundError:
        print("Piecewise quadratic thinning not found")
        pw_q_r_exists = False
    try:
        c_u_data = pd.read_csv(f"{DATAPATH}/flux_cox_uniform.csv")
    except FileNotFoundError:
        print("Uniform Cox not found")
        c_u_exists = False
    try:
        c_e_data = pd.read_csv(f"{DATAPATH}/flux_cox_exponential.csv")
    except FileNotFoundError:
        print("Exponential Cox not found")
        c_e_exists = False
    try:
        c_g_data = pd.read_csv(f"{DATAPATH}/flux_cox_gaussian.csv")
    except FileNotFoundError:
        print("Gaussian Cox not found")
        c_g_exists = False
    try:
        c_u_m_data = pd.read_csv(f"{DATAPATH}/flux_cox_uniform_monosample.csv")
    except FileNotFoundError:
        print("Uniform Cox monosample not found")
        c_u_m_exists = False
    try:
        c_e_m_data = pd.read_csv(f"{DATAPATH}/flux_cox_exponential_monosample.csv")
    except FileNotFoundError:
        print("Exponential Cox monosample not found")
        c_e_m_exists = False
    try:
        c_g_m_data = pd.read_csv(f"{DATAPATH}/flux_cox_gaussian_monosample.csv")
    except FileNotFoundError:
        print("Gaussian Cox monosample not found")
        c_g_m_exists = False
    try:
        con_data = pd.read_csv(f"{DATAPATH}/flux_constant.csv")
    except FileNotFoundError:
        print("Constant not found")
        con_exists = False
    if r_exists:
        r_obj = data(r_data, ylabel)
        if flux:
            r_obj.plot("Linear Thinning Flux", "flux_linear_thinning")
        else:
            r_obj.plot("Linear Thinning Volume Fraction", "vf_linear_thinning")
    if l_exists:
        l_obj = data(l_data, ylabel)
        if flux:
            l_obj.plot("Linear Direct Flux", "flux_linear")
        else:
            l_obj.plot("Linear Direct Volume Fraction", "vf_linear")
    if q_r_exists:
        q_r_obj = data(q_r_data, ylabel)
        if flux:
            q_r_obj.plot("Quadratic Thinning Flux", "flux_quad_thinning")
        else:
            q_r_obj.plot("Quadratic Thinning Volume Fraction", "vf_quad_thinning")
    if pw_exists:
        pw_obj = data(pw_data, ylabel)
        if flux:
            pw_obj.plot("Linear Piecewise Constant Flux", "flux_piecewise_linear")
        else:
            pw_obj.plot("Linear Piecewise Constant Volume Fraction", "vf_piecewise_linear")
    if pw_r_exists:
        pw_r_obj = data(pw_r_data, ylabel)
        if flux:
            pw_r_obj.plot("Linear Piecewise Constant Thinning Flux", "flux_piecewise_linear_thinning")
        else:
            pw_r_obj.plot("Linear Piecewise Constant Thinning Volume Fraction", "vf_piecewise_linear_thinning")
    if pw_q_exists:
        pw_q_obj = data(pw_q_data, ylabel)
        if flux:
            pw_q_obj.plot("Quadratic Piecewise Constant Flux", "flux_piecewise_quad")
        else:
            pw_q_obj.plot("Quadratic Piecewise Constant Volume Fraction", "vf_piecewise_quad")
    if pw_q_r_exists:
        pw_q_r_obj = data(pw_q_r_data, ylabel)
        if flux:
            pw_q_r_obj.plot("Quadratic Piecewise Constant Thinning Flux", "flux_piecewise_quad_thinning")
        else:
            pw_q_r_obj.plot("Quadratic Piecewise Constant Thinning Volume Fraction", "vf_piecewise_quad_thinning")
    if c_u_exists:
        c_u_obj = data(c_u_data, ylabel)
        if flux:
            c_u_obj.plot("Uniform Cox Flux", "cox_uniform")
        else:
            c_u_obj.plot("Uniform Cox Volume Fraction", "vf_cox_uniform")
    if c_e_exists:
        c_e_obj = data(c_e_data, ylabel)
        if flux:
            c_e_obj.plot("Exponential Cox Flux", "cox_exponential")
        else:
            c_e_obj.plot("Exponential Cox Volume Fraction", "vf_cox_exponential")
    if c_g_exists:
        c_g_obj = data(c_g_data, ylabel)
        if flux:
            c_g_obj.plot("Gaussian Cox Flux", "cox_gaussian")
        else:
            c_g_obj.plot("Gaussian Cox Volume Fraction", "vf_cox_gaussian")
    if c_u_m_exists:
        c_u_m_obj = data(c_u_m_data, ylabel)
        if flux:
            c_u_m_obj.plot("Uniform Cox Monosample Flux", "cox_uniform_monosample")
        else:
            c_u_m_obj.plot("Uniform Cox Monosample Volume Fraction", "vf_cox_uniform_monosample")
    if c_e_m_exists:
        c_e_m_obj = data(c_e_m_data, ylabel)
        if flux:
            c_e_m_obj.plot("Exponential Cox Monosample Flux", "cox_exponential_monosample")
        else:
            c_e_m_obj.plot("Exponential Cox Volume Fraction", "vf_cox_exponential_monosample")
    if c_g_m_exists:
        c_g_m_obj = data(c_g_m_data, ylabel)
        if flux:
            c_g_m_obj.plot("Gaussian Cox Monosample Flux", "cox_gaussian_monosample")
        else:
            c_g_m_obj.plot("Gaussian Cox Volume Fraction", "vf_cox_gaussian_monosample")
    if con_exists:
        con_obj = data(con_data, ylabel)
        if flux:
            con_obj.plot("Constant Flux", "constant")
        else:
            con_obj.plot("Constat Volume Fraction", "vf_constant")

if __name__ == '__main__':
    main()
