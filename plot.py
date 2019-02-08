#!/usr/bin/env python3
'''
Read data and create figures
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas

PLOTPATH = "./out/nonlinear"
CSVPATH = f"{PLOTPATH}/data"

def main():
    '''Main function'''
    h_exists, he_exists, exists, l_exists, mc_exists = True, True, True, True, True
    #time,intensity,temperature
    try:
        homog_nonlinear = pandas.read_csv(f"{CSVPATH}/homog_nonlinear.csv")
    except FileNotFoundError:
        h_exists = False
    #time,intensity,temperature
    try:
        homog_nonlinear_exp = pandas.read_csv(f"{CSVPATH}/homog_nonlinear_exp.csv")
    except FileNotFoundError:
        he_exists = False
    #time,intensity1,varintensity1,temperature1,vartemperature1,intensity2,varintensity2,temperature2,vartemperature2
    try:
        nonlinear = pandas.read_csv(f"{CSVPATH}/nonlinear.csv")
    except FileNotFoundError:
        exists = False
    #time,intensity1,temperature1,intensity2,temperature2
    try:
        linear = pandas.read_csv(f"{CSVPATH}/linear.csv")
    except FileNotFoundError:
        l_exists = False
    #time,intensity1,varintensity1,temperature1,vartemperature1,intensity2,varintensity2,temperature2,vartemperature2
    try:
        nonlinear_mc = pandas.read_csv(f"{CSVPATH}/nonlinearmc.csv")
    except FileNotFoundError:
        mc_exists = False

    # IMPLICIT PLOTS
    if h_exists:
        # Homogeneous Intensity
        plt.plot(homog_nonlinear['time'], homog_nonlinear['intensity'], color='m')
        plt.title("Homogeneous Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/homog_intensity.png")
        plt.cla()
        plt.clf()

        # Homogeneous Temperature
        plt.plot(homog_nonlinear['time'], homog_nonlinear['temperature'], color='g')
        plt.title("Homogeneous Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/homog_temperature.png")
        plt.cla()
        plt.clf()

    # EXPLICIT PLOTS
    if he_exists:
        # Explicit Homogeneous Intensity
        plt.plot(homog_nonlinear_exp['time'], homog_nonlinear_exp['intensity'], color='m')
        plt.title("Explicit Homogeneous Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/exp_homog_intensity.png")
        plt.cla()
        plt.clf()

        # Explicit Homogeneous Temperature
        plt.plot(homog_nonlinear_exp['time'], homog_nonlinear_exp['temperature'], color='g')
        plt.title("Explicit Homogeneous Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/exp_homog_temperature.png")
        plt.cla()
        plt.clf()

    # IMPLICIT AND EXPLICIT COMPARISON
    if h_exists and he_exists:
        # Intensity
        plt.plot(homog_nonlinear['time'], homog_nonlinear['intensity'], color='m', label="Implicit")
        plt.plot(homog_nonlinear_exp['time'], homog_nonlinear_exp['intensity'], color='m', linestyle=':', label="Explicit")
        plt.title("Implicit and Explicit Homogeneous Intensity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/imp_exp_homog_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(homog_nonlinear['time'], homog_nonlinear['temperature'], color='g', label="Implicit")
        plt.plot(homog_nonlinear_exp['time'], homog_nonlinear_exp['temperature'], color='g', linestyle=':', label="Implicit")
        plt.title("Implicit and Explicit Homogeneous Temperature")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/imp_exp_homog_temperature.png")
        plt.cla()
        plt.clf()

    # REALIZATIONS PLOTS
    if exists:
        # Standard deviation computing
        nl_std_intensity_1 = nonlinear['varintensity1'].apply(np.sqrt)
        nl_std_intensity_2 = nonlinear['varintensity2'].apply(np.sqrt)
        nl_std_temp_1 = nonlinear['vartemperature1'].apply(np.sqrt)
        nl_std_temp_2 = nonlinear['vartemperature2'].apply(np.sqrt)
        nl_lb_intensity_1 = nonlinear['intensity1'] - nl_std_intensity_1
        nl_ub_intensity_1 = nonlinear['intensity1'] + nl_std_intensity_1
        nl_lb_intensity_2 = nonlinear['intensity2'] - nl_std_intensity_2
        nl_ub_intensity_2 = nonlinear['intensity2'] + nl_std_intensity_2
        nl_lb_temp_1 = nonlinear['temperature1'] - nl_std_temp_1
        nl_ub_temp_1 = nonlinear['temperature1'] + nl_std_temp_1
        nl_lb_temp_2 = nonlinear['temperature2'] - nl_std_temp_2
        nl_ub_temp_2 = nonlinear['temperature2'] + nl_std_temp_2

        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2")
        plt.plot(nonlinear['time'], nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_intensity_1, nl_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(nonlinear['time'], nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_intensity_2, nl_ub_intensity_2, color='r', alpha=0.5)
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2")
        plt.plot(nonlinear['time'], nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_temp_1, nl_ub_temp_1, color='b', alpha=0.5)
        plt.plot(nonlinear['time'], nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_temp_2, nl_ub_temp_2, color='r', alpha=0.5)
        plt.title("Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature.png")
        plt.cla()
        plt.clf()

    # MONTE CARLO PLOTS
    if mc_exists:
        # Standard deviation computing
        mc_nl_std_intensity_1 = nonlinear_mc['varintensity1'].apply(np.sqrt)
        mc_nl_std_intensity_2 = nonlinear_mc['varintensity2'].apply(np.sqrt)
        mc_nl_std_temp_1 = nonlinear_mc['vartemperature1'].apply(np.sqrt)
        mc_nl_std_temp_2 = nonlinear_mc['vartemperature2'].apply(np.sqrt)
        mc_nl_lb_intensity_1 = nonlinear_mc['intensity1'] - mc_nl_std_intensity_1
        mc_nl_ub_intensity_1 = nonlinear_mc['intensity1'] + mc_nl_std_intensity_1
        mc_nl_lb_intensity_2 = nonlinear_mc['intensity2'] - mc_nl_std_intensity_2
        mc_nl_ub_intensity_2 = nonlinear_mc['intensity2'] + mc_nl_std_intensity_2
        mc_nl_lb_temp_1 = nonlinear_mc['temperature1'] - mc_nl_std_temp_1
        mc_nl_ub_temp_1 = nonlinear_mc['temperature1'] + mc_nl_std_temp_1
        mc_nl_lb_temp_2 = nonlinear_mc['temperature2'] - mc_nl_std_temp_2
        mc_nl_ub_temp_2 = nonlinear_mc['temperature2'] + mc_nl_std_temp_2

        # Intensity
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity1'], color='b', label="Material 1")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity2'], color='r', label="Material 2")
        plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear_mc['time'], mc_nl_lb_intensity_1, mc_nl_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear_mc['time'], mc_nl_lb_intensity_2, mc_nl_ub_intensity_2, color='r', alpha=0.5)
        plt.title("Monte Carlo Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity_mc.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature1'], color='b', label="Material 1")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature2'], color='r', label="Material 2")
        plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear_mc['time'], mc_nl_lb_temp_1, mc_nl_ub_temp_1, color='b', alpha=0.5)
        plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear_mc['time'], mc_nl_lb_temp_2, mc_nl_ub_temp_2, color='r', alpha=0.5)
        plt.title("Monte Carlo Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_mc.png")
        plt.cla()
        plt.clf()

    if exists and mc_exists:
        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1 - Realization")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2 - Realization")
        plt.plot(nonlinear['time'], nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity1'], color='c', label="Material 1 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity2'], color='m', label="Material 2 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_2, color='m', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_2, color='m', linestyle=':', label=None)
        plt.title("Intensity Plot - Realizations & Monte Carlo")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity_std_comp.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1 - Realizations")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2 - Realizations")
        plt.plot(nonlinear['time'], nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature1'], color='c', label="Material 1 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature2'], color='m', label="Material 2 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_2, color='m', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_2, color='m', linestyle=':', label=None)
        plt.title("Temperature Plot - Realizations & Monte Carlo")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_std_comp.png")
        plt.cla()
        plt.clf()

    # LINEAR PLOTS
    if l_exists:
        # Intensity
        plt.plot(linear['time'], linear['intensity1'], color='b', label="Material 1")
        plt.plot(linear['time'], linear['intensity2'], color='r', label="Material 2")
        plt.title("Linear Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/linear_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(linear['time'], linear['temperature1'], color='b', label="Material 1")
        plt.plot(linear['time'], linear['temperature2'], color='r', label="Material 2")
        plt.title("Linear Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/linear_temperature.png")
        plt.cla()
        plt.clf()

    # REALIZATIONS VS MC COMPARISON
    if exists and h_exists:
        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1 - Realizations")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2 - Realizations")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity1'], color='b', linestyle='--', label="Material 1 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity2'], color='r', linestyle='--', label="Material 2 - Monte Carlo")
        plt.title("Realizations and Monte Carlo Intensity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/realization_mc_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1 - Realizations")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2 - Realizations")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature1'], color='b', linestyle='--', label="Material 1 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature2'], color='r', linestyle='--', label="Material 2 - Monte Carlo")
        plt.title("Realizations and Monte Carlo Temperature")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/realization_mc_temperature.png")
        plt.cla()
        plt.clf()

    # REALIZATIONS COMPARISON
    if exists and h_exists:
        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2")
        plt.plot(homog_nonlinear['time'], homog_nonlinear['intensity'], color='m', linestyle='--', label="Homogeneous")
        plt.plot(nonlinear['time'], nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_intensity_1, nl_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(nonlinear['time'], nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_intensity_2, nl_ub_intensity_2, color='r', alpha=0.5)
        plt.title("Intensity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity_comp.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2")
        plt.plot(homog_nonlinear['time'], homog_nonlinear['temperature'], color='g', linestyle='--', label="Homogeneous")
        plt.plot(nonlinear['time'], nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_temp_1, nl_ub_temp_1, color='b', alpha=0.5)
        plt.plot(nonlinear['time'], nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_temp_2, nl_ub_temp_2, color='r', alpha=0.5)
        plt.title("Temperature")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_comp.png")
        plt.cla()
        plt.clf()

    # FULL COMPARISON
    if exists and h_exists and l_exists:
        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1 - Nonlinear")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2 - Nonlinear")
        plt.plot(homog_nonlinear['time'], homog_nonlinear['intensity'], color='m', linestyle='--', label="Homogeneous - Nonlinear")
        plt.plot(linear['time'], linear['intensity1'], color='b', linestyle='-.', label="Material 1 - Linear")
        plt.plot(linear['time'], linear['intensity2'], color='r', linestyle='-.', label="Material 2 - Linear")
        plt.plot(nonlinear['time'], nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_intensity_1, nl_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(nonlinear['time'], nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_intensity_2, nl_ub_intensity_2, color='r', alpha=0.5)
        plt.title("Intensity Comparison")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity_full_comp.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1 - Nonlinear")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2 - Nonlinear")
        plt.plot(homog_nonlinear['time'], homog_nonlinear['temperature'], color='g', linestyle='--', label="Homogeneous - Nonlinear")
        plt.plot(linear['time'], linear['temperature1'], color='b', linestyle='-.', label="Material 1 - Linear")
        plt.plot(linear['time'], linear['temperature2'], color='r', linestyle='-.', label="Material 2 - Linear")
        plt.plot(nonlinear['time'], nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_temp_1, nl_ub_temp_1, color='b', alpha=0.5)
        plt.plot(nonlinear['time'], nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_temp_2, nl_ub_temp_2, color='r', alpha=0.5)
        plt.title("Temperature Comparison")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_full_comp.png")
        plt.cla()
        plt.clf()

if __name__ == '__main__':
    main()
