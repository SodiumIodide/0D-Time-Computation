#!/usr/bin/env python3
'''
Read data and create figures
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas

PLOTPATH = "./out/nonlinear"
CSVPATH = f"{PLOTPATH}/data"
HISTPATH = f"{PLOTPATH}/pdf_data"

def main():
    '''Main function'''
    h_exists, he_exists, exists, ex_exists, l_exists, mc_exists, mc_hist, re_hist = True, True, True, True, True, True, True, True
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
    #time,intensity1,varintensity1,maxintensity1,minintensity1,temperature1,vartemperature1,maxtemperature1,mintemperature1,intensity2,varintensity2,maxintensity2,minintensity2,temperature2,vartemperature2,maxtemperature2,mintemperature2
    try:
        nonlinear = pandas.read_csv(f"{CSVPATH}/nonlinear.csv")
    except FileNotFoundError:
        exists = False
    #time,intensity1,varintensity1,maxintensity1,minintensity1,temperature1,vartemperature1,maxtemperature1,mintemperature1,intensity2,varintensity2,maxintensity2,minintensity2,temperature2,vartemperature2,maxtemperature2,mintemperature2
    try:
        nonlinear_exp = pandas.read_csv(f"{CSVPATH}/nonlinear_exp.csv")
    except FileNotFoundError:
        ex_exists = False
    #time,intensity1,temperature1,intensity2,temperature2
    try:
        linear = pandas.read_csv(f"{CSVPATH}/linear.csv")
    except FileNotFoundError:
        l_exists = False
    #time,intensity1,varintensity1,maxintensity1,minintensity1,temperature1,vartemperature1,maxtemperature1,mintemperature1,intensity2,varintensity2,maxintensity2,minintensity2,temperature2,vartemperature2,maxtemperature2,mintemperature2
    try:
        nonlinear_mc = pandas.read_csv(f"{CSVPATH}/nonlinearmc.csv")
    except FileNotFoundError:
        mc_exists = False
    #time,intensity1arr,freqintensity1,intensity2arr,freqintensity2,temperature1arr,freqtemperature1,temperature2arr,freqtemperature2,opacity1arr,freqopacity1,opacity2arr,freqopacity2
    try:
        mc_hist_data = pandas.read_csv(f"{HISTPATH}/mc_pdf.csv")
    except FileNotFoundError:
        mc_hist = False
    #time,intensity1arr,freqintensity1,intensity2arr,freqintensity2,temperature1arr,freqtemperature1,temperature2arr,freqtemperature2,opacity1arr,opacity2arr
    try:
        re_hist_data = pandas.read_csv(f"{HISTPATH}/realizations_pdf.csv")
    except FileNotFoundError:
        re_hist = False

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
        plt.ylabel("Intensity (erg/cm^2-s)")
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

    # REALIZATIONS EXPLICIT PLOTS
    if ex_exists:
        # Standard deviation computing
        ex_nl_std_intensity_1 = nonlinear_exp['varintensity1'].apply(np.sqrt)
        ex_nl_std_intensity_2 = nonlinear_exp['varintensity2'].apply(np.sqrt)
        ex_nl_std_temp_1 = nonlinear_exp['vartemperature1'].apply(np.sqrt)
        ex_nl_std_temp_2 = nonlinear_exp['vartemperature2'].apply(np.sqrt)
        ex_nl_lb_intensity_1 = nonlinear_exp['intensity1'] - ex_nl_std_intensity_1
        ex_nl_ub_intensity_1 = nonlinear_exp['intensity1'] + ex_nl_std_intensity_1
        ex_nl_lb_intensity_2 = nonlinear_exp['intensity2'] - ex_nl_std_intensity_2
        ex_nl_ub_intensity_2 = nonlinear_exp['intensity2'] + ex_nl_std_intensity_2
        ex_nl_lb_temp_1 = nonlinear_exp['temperature1'] - ex_nl_std_temp_1
        ex_nl_ub_temp_1 = nonlinear_exp['temperature1'] + ex_nl_std_temp_1
        ex_nl_lb_temp_2 = nonlinear_exp['temperature2'] - ex_nl_std_temp_2
        ex_nl_ub_temp_2 = nonlinear_exp['temperature2'] + ex_nl_std_temp_2

        # Intensity
        plt.plot(nonlinear_exp['time'], nonlinear_exp['intensity1'], color='b', label="Material 1")
        plt.plot(nonlinear_exp['time'], nonlinear_exp['intensity2'], color='r', label="Material 2")
        plt.plot(nonlinear_exp['time'], ex_nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear_exp['time'], ex_nl_lb_intensity_1, ex_nl_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(nonlinear_exp['time'], ex_nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear_exp['time'], ex_nl_lb_intensity_2, ex_nl_ub_intensity_2, color='r', alpha=0.5)
        plt.title("Explicit Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity_exp.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear_exp['time'], nonlinear_exp['temperature1'], color='b', label="Material 1")
        plt.plot(nonlinear_exp['time'], nonlinear_exp['temperature2'], color='r', label="Material 2")
        plt.plot(nonlinear_exp['time'], ex_nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear_exp['time'], ex_nl_lb_temp_1, ex_nl_ub_temp_1, color='b', alpha=0.5)
        plt.plot(nonlinear_exp['time'], ex_nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear_exp['time'], ex_nl_lb_temp_2, ex_nl_ub_temp_2, color='r', alpha=0.5)
        plt.title("Explicit Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_exp.png")
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
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1 - Realizations")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2 - Realizations")
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

    if ex_exists and mc_exists:
        # Intensity
        plt.plot(nonlinear_exp['time'], nonlinear_exp['intensity1'], color='b', label="Material 1 - Realizations Explicit")
        plt.plot(nonlinear_exp['time'], nonlinear_exp['intensity2'], color='r', label="Material 2 - Realizations Explicit")
        plt.plot(nonlinear_exp['time'], ex_nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity1'], color='c', label="Material 1 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity2'], color='m', label="Material 2 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_2, color='m', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_2, color='m', linestyle=':', label=None)
        plt.title("Intensity Plot - Realizations Explicit & Monte Carlo")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity_std_comp_exp.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear_exp['time'], nonlinear_exp['temperature1'], color='b', label="Material 1 - Realizations Explicit")
        plt.plot(nonlinear_exp['time'], nonlinear_exp['temperature2'], color='r', label="Material 2 - Realizations Explicit")
        plt.plot(nonlinear_exp['time'], ex_nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature1'], color='c', label="Material 1 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature2'], color='m', label="Material 2 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_2, color='m', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_2, color='m', linestyle=':', label=None)
        plt.title("Temperature Plot - Realizations Explicit & Monte Carlo")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_std_comp_exp.png")
        plt.cla()
        plt.clf()

    if exists and ex_exists:
        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1 - Realizations")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2 - Realizations")
        plt.plot(nonlinear['time'], nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], nonlinear_exp['intensity1'], color='c', label="Material 1 - Realizations Explicit")
        plt.plot(nonlinear_exp['time'], nonlinear_exp['intensity2'], color='m', label="Material 2 - Realizations Explicit")
        plt.plot(nonlinear_exp['time'], ex_nl_lb_intensity_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_intensity_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_lb_intensity_2, color='m', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_intensity_2, color='m', linestyle=':', label=None)
        plt.title("Intensity Plot - Realizations Implicit & Explicit")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity_imp_exp_std_comp.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1 - Realizations")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2 - Realizations")
        plt.plot(nonlinear['time'], nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], nonlinear_exp['temperature1'], color='c', label="Material 1 - Realizations Explicit")
        plt.plot(nonlinear_exp['time'], nonlinear_exp['temperature2'], color='m', label="Material 2 - Realizations Explicit")
        plt.plot(nonlinear_exp['time'], ex_nl_lb_temp_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_temp_1, color='c', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_lb_temp_2, color='m', linestyle=':', label=None)
        plt.plot(nonlinear_exp['time'], ex_nl_ub_temp_2, color='m', linestyle=':', label=None)
        plt.title("Temperature Plot - Realizations Implicit & Explicit")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_imp_exp_std_comp.png")
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

    # HISTOGRAM PLOTS
    if mc_hist:
        time_ss = mc_hist_data['time'][0]

        # Intensity Material 1
        plt.plot(mc_hist_data['intensity1arr'][1:len(mc_hist_data['intensity1arr'])-1], mc_hist_data['freqintensity1'][1:len(mc_hist_data['freqintensity1'])-1])
        plt.title(f"Intensity 1 Histogram - Monte Carlo (Time={time_ss:.4E} ct)")
        plt.xlabel("Intensity (erg/cm^2-s)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/mc_intensity_1_hist.png")
        plt.cla()
        plt.clf()

        # Intensity Material 2
        plt.plot(mc_hist_data['intensity2arr'][1:len(mc_hist_data['intensity2arr'])-1], mc_hist_data['freqintensity2'][1:len(mc_hist_data['freqintensity2'])-1])
        plt.title(f"Intensity 2 Histogram - Monte Carlo (Time={time_ss:.4E} ct)")
        plt.xlabel("Intensity (erg/cm^2-s)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/mc_intensity_2_hist.png")
        plt.cla()
        plt.clf()

        # Temperature Material 1
        plt.plot(mc_hist_data['temperature1arr'][1:len(mc_hist_data['temperature1arr'])-1], mc_hist_data['freqtemperature1'][1:len(mc_hist_data['freqtemperature1'])-1])
        plt.title(f"Temperature 1 Histogram - Monte Carlo (Time={time_ss:.4E} ct)")
        plt.xlabel("Temperature (eV)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/mc_temperature_1_hist.png")
        plt.cla()
        plt.clf()

        # Temperature Material 2
        plt.plot(mc_hist_data['temperature2arr'][1:len(mc_hist_data['temperature2arr'])-1], mc_hist_data['freqtemperature2'][1:len(mc_hist_data['freqtemperature2'])-1])
        plt.title(f"Temperature 2 Histogram - Monte Carlo (Time={time_ss:.4E} ct)")
        plt.xlabel("Temperature (eV)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/mc_temperature_2_hist.png")
        plt.cla()
        plt.clf()

        # Opacity Material 1
        plt.plot(mc_hist_data['opacity1arr'][1:len(mc_hist_data['opacity1arr'])-1], mc_hist_data['freqopacity1'][1:len(mc_hist_data['freqopacity1'])-1])
        plt.title(f"Opacity 1 Histogram - Monte Carlo (Time={time_ss:.4E} ct)")
        plt.xlabel("Opacity (cm^-1)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/mc_opacity_1_hist.png")
        plt.cla()
        plt.clf()

        # Opacity Material 2
        plt.plot(mc_hist_data['opacity2arr'][1:len(mc_hist_data['opacity2arr'])-1], mc_hist_data['freqopacity2'][1:len(mc_hist_data['freqopacity2'])-1])
        plt.title(f"Opacity 2 Histogram - Monte Carlo (Time={time_ss:.4E} ct)")
        plt.xlabel("Opacity (cm^-1)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/mc_opacity_2_hist.png")
        plt.cla()
        plt.clf()

    if re_hist:
        # Intensity Material 1
        plt.plot(re_hist_data['intensity1arr'], re_hist_data['freqintensity1'])
        plt.title("Intensity 1 Histogram - Realizations")
        plt.xlabel("Intensity (erg/cm^2-s)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/realizations_intensity_1_hist.png")
        plt.cla()
        plt.clf()

        # Intensity Material 2
        plt.plot(re_hist_data['intensity2arr'], re_hist_data['freqintensity2'])
        plt.title("Intensity 2 Histogram - Realizations")
        plt.xlabel("Intensity (erg/cm^2-s)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/realizations_intensity_2_hist.png")
        plt.cla()
        plt.clf()

        # Temperature Material 1
        plt.plot(re_hist_data['temperature1arr'], re_hist_data['freqtemperature1'])
        plt.title("Temperature 1 Histogram - Realizations")
        plt.xlabel("Temperature (eV)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/realizations_temperature_1_hist.png")
        plt.cla()
        plt.clf()

        # Temperature Material 2
        plt.plot(re_hist_data['temperature2arr'], re_hist_data['freqtemperature2'])
        plt.title("Temperature 2 Histogram - Realizations")
        plt.xlabel("Temperature (eV)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/realizations_temperature_2_hist.png")
        plt.cla()
        plt.clf()

        # Opacity Material 1
        plt.plot(re_hist_data['opacity1arr'], re_hist_data['freqopacity1'])
        plt.title("Opacity 1 Histogram - Realizations")
        plt.xlabel("Opacity (cm^-1)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/realizations_opacity_1_hist.png")
        plt.cla()
        plt.clf()

        # Opacity Material 2
        plt.plot(re_hist_data['opacity2arr'], re_hist_data['freqopacity2'])
        plt.title("Opacity 2 Histogram - Realizations")
        plt.xlabel("Opacity (cm^-1)")
        plt.ylabel("Frequency")
        #plt.xscale("log")
        #plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/realizations_opacity_2_hist.png")
        plt.cla()
        plt.clf()

if __name__ == '__main__':
    main()
