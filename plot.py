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
    h_exists, he_exists, exists, l_exists = True, True, True, True
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
        plt.savefig(f"{PLOTPATH}/temperature.png")
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
        plt.savefig(f"{PLOTPATH}/linear_temperature.png")
        plt.cla()
        plt.clf()

    # REALIZATIONS COMPARISON
    if exists and h_exists:
        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2")
        plt.plot(homog_nonlinear['time'], homog_nonlinear['intensity'], color='m', linestyle=':', label="Homogeneous")
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
        plt.savefig(f"{PLOTPATH}/intensity_comp.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2")
        plt.plot(homog_nonlinear['time'], homog_nonlinear['temperature'], color='g', linestyle=':', label="Homogeneous")
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
        plt.savefig(f"{PLOTPATH}/temperature_comp.png")
        plt.cla()
        plt.clf()

    # FULL COMPARISON
    if exists and h_exists and l_exists:
        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1 - Nonlinear")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2 - Nonlinear")
        plt.plot(homog_nonlinear['time'], homog_nonlinear['intensity'], color='m', linestyle=':', label="Homogeneous - Nonlinear")
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
        plt.savefig(f"{PLOTPATH}/intensity_full_comp.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1 - Nonlinear")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2 - Nonlinear")
        plt.plot(homog_nonlinear['time'], homog_nonlinear['temperature'], color='g', linestyle=':', label="Homogeneous - Nonlinear")
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
        plt.savefig(f"{PLOTPATH}/temperature_full_comp.png")
        plt.cla()
        plt.clf()

if __name__ == '__main__':
    main()
