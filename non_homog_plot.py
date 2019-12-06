#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

XSCALE = "log"
YSCALE = "log"
PLOTPATH = "./out/nonlinear"
DATAPATH = f"{PLOTPATH}/data"

class data:
    def __init__(self, in_data):
        # Standard deviation computing
        std_intensity_1 = in_data['varintensity1'].apply(np.abs).apply(np.sqrt)
        std_intensity_2 = in_data['varintensity2'].apply(np.abs).apply(np.sqrt)
        std_temp_1 = in_data['vartemperature1'].apply(np.abs).apply(np.sqrt)
        std_temp_2 = in_data['vartemperature2'].apply(np.abs).apply(np.sqrt)
        std_opacity_1 = in_data['varopacity1'].apply(np.abs).apply(np.sqrt)
        std_opacity_2 = in_data['varopacity2'].apply(np.abs).apply(np.sqrt)
        self.intensity_1 = in_data['intensity1']
        self.intensity_2 = in_data['intensity2']
        self.temp_1 = in_data['temperature1']
        self.temp_2 = in_data['temperature2']
        self.opacity_1 = in_data['opacity1']
        self.opacity_2 = in_data['opacity2']
        self.lb_intensity_1 = self.intensity_1 - std_intensity_1
        self.ub_intensity_1 = self.intensity_1 + std_intensity_1
        self.lb_intensity_2 = self.intensity_2 - std_intensity_2
        self.ub_intensity_2 = self.intensity_2 + std_intensity_2
        self.lb_temp_1 = self.temp_1 - std_temp_1
        self.ub_temp_1 = self.temp_1 + std_temp_1
        self.lb_temp_2 = self.temp_2 - std_temp_2
        self.ub_temp_2 = self.temp_2 + std_temp_2
        self.lb_opacity_1 = self.opacity_1 - std_opacity_1
        self.ub_opacity_1 = self.opacity_1 + std_opacity_1
        self.lb_opacity_2 = self.opacity_2 - std_opacity_2
        self.ub_opacity_2 = self.opacity_2 + std_opacity_2
        self.time = in_data['time']

    def plot(self, plotname, filename):
        # Plot intensity
        plt.plot(self.time, self.intensity_1, color='b', label="Material 1")
        plt.plot(self.time, self.intensity_2, color='r', label="Material 2")
        plt.fill_between(self.time, self.lb_intensity_1, self.ub_intensity_1, color='b', alpha=0.5)
        plt.fill_between(self.time, self.lb_intensity_2, self.ub_intensity_2, color='r', alpha=0.5)
        plt.title(f"{plotname} Intensity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel(r"Intensity (erg/cm$^2$-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which='both', axis='both')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/{filename}_intensity.png")
        plt.cla()
        plt.clf()
        # Plot temperature
        plt.plot(self.time, self.temp_1, color='b', label="Material 1")
        plt.plot(self.time, self.temp_2, color='r', label="Material 2")
        plt.fill_between(self.time, self.lb_temp_1, self.ub_temp_1, color='b', alpha=0.5)
        plt.fill_between(self.time, self.lb_temp_2, self.ub_temp_2, color='r', alpha=0.5)
        plt.title(f"{plotname} Temperature")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which='both', axis='both')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/{filename}_temperature.png")
        plt.cla()
        plt.clf()
        # Plot opacity
        plt.plot(self.time, self.opacity_1, color='b', label="Material 1")
        plt.plot(self.time, self.opacity_2, color='r', label="Material 2")
        plt.fill_between(self.time, self.lb_opacity_1, self.ub_opacity_1, color='b', alpha=0.5)
        plt.fill_between(self.time, self.lb_opacity_2, self.ub_opacity_2, color='r', alpha=0.5)
        plt.title(f"{plotname} Opacity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel(r"Opacity (cm$^{-1}$)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which='both', axis='both')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/{filename}_opacity.png")
        plt.cla()
        plt.clf()


def main():
    '''Main function'''
    mc_quad_exists = mc_linear_exists = quad_exists = linear_exists = True
    try:
        mc_quad_data = pd.read_csv(f"{DATAPATH}/nonlinearmc_nonhomog_quad.csv")
    except FileNotFoundError:
        print("Monte Carlo Nonhomogeneous Quadratic not found")
        mc_quad_exists = False
    try:
        mc_linear_data = pd.read_csv(f"{DATAPATH}/nonlinearmc_nonhomog_linear.csv")
    except FileNotFoundError:
        print("Monte Carlo Nonhomogeneous Linear not found")
        mc_linear_exists = False
    try:
        quad_data = pd.read_csv(f"{DATAPATH}/nonlinear_nonhomog_quad.csv")
    except FileNotFoundError:
        print("Nonhomogeneous Quadratic not found")
        quad_exists = False
    try:
        linear_data = pd.read_csv(f"{DATAPATH}/nonlinear_nonhomog_linear.csv")
    except FileNotFoundError:
        print("Nonhomogeneous Linear not found")
        linear_exists = False
    if mc_quad_exists:
        mc_quad_obj = data(mc_quad_data)
        mc_quad_obj.plot("Monte Carlo Nonhomogeneous Quadratic", "mc_quad")
    if mc_linear_exists:
        mc_linear_obj = data(mc_linear_data)
        mc_linear_obj.plot("Monte Carlo Nonhomogeneous Linear", "mc_linear")
    if quad_exists:
        quad_obj = data(quad_data)
        quad_obj.plot("Nonhomogeneous Quadratic", "num_quad")
    if linear_exists:
        linear_obj = data(linear_data)
        linear_obj.plot("Nonhomogeneous Linear", "num_linear")

if __name__ == '__main__':
    main()
