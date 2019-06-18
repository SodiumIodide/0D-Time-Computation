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
    exists = mc_exists = heur_exists = True
    #time,intensity1,varintensity1,maxintensity1,minintensity1,temperature1,vartemperature1,maxtemperature1,mintemperature1,intensity2,varintensity2,maxintensity2,minintensity2,temperature2,vartemperature2,maxtemperature2,mintemperature2
    try:
        nonlinear = pandas.read_csv(f"{CSVPATH}/nonlinear.csv")
    except FileNotFoundError:
        exists = False
    #time,intensity1,varintensity1,maxintensity1,minintensity1,temperature1,vartemperature1,maxtemperature1,mintemperature1,intensity2,varintensity2,maxintensity2,minintensity2,temperature2,vartemperature2,maxtemperature2,mintemperature2
    try:
        nonlinear_mc = pandas.read_csv(f"{CSVPATH}/nonlinearmc.csv")
    except FileNotFoundError:
        mc_exists = False
    #time,intensity1,temperature1,intensity2,temperature2
    try:
        heur_data = pandas.read_csv(f"{CSVPATH}/heuristic.csv")
    except FileNotFoundError:
        heur_exists = False

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

    # REALIZATIONS VS HEURISTIC MODEL
    if exists and heur_exists:
        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1 - Realizations")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2 - Realizations")
        plt.plot(heur_data['time'], heur_data['intensity1'], color='c', linestyle=':', label="Material 1 - Heuristic")
        plt.plot(heur_data['time'], heur_data['intensity2'], color='m', linestyle=':', label="Material 2 - Heuristic")
        plt.title("Intensity - Realizations vs. Heuristic")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intenisty (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/heur_real_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1 - Realizations")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2 - Realizations")
        plt.plot(heur_data['time'], heur_data['temperature1'], color='c', linestyle=':', label="Material 1 - Heuristic")
        plt.plot(heur_data['time'], heur_data['temperature2'], color='m', linestyle=':', label="Material 2 - Heuristic")
        plt.title("Temperature - Realizations vs. Heuristic")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/heur_real_temperature.png")
        plt.cla()
        plt.clf()

        # Compute relative error
        if len(nonlinear['time']) == len(heur_data['time']):
            intensity_1_err = np.zeros(len(nonlinear['time']))
            intensity_2_err = np.zeros(len(nonlinear['time']))
            temp_1_err = np.zeros(len(nonlinear['time']))
            temp_2_err = np.zeros(len(nonlinear['time']))
            for index, _ in enumerate(nonlinear['time']):
                intensity_1_err[index] = np.abs(nonlinear['intensity1'][index] - heur_data['intensity1'][index]) / nonlinear['intensity1'][index]
                intensity_2_err[index] = np.abs(nonlinear['intensity2'][index] - heur_data['intensity2'][index]) / nonlinear['intensity1'][index]
                temp_1_err[index] = np.abs(nonlinear['temperature1'][index] - heur_data['temperature1'][index]) / nonlinear['temperature1'][index]
                temp_2_err[index] = np.abs(nonlinear['temperature2'][index] - heur_data['temperature2'][index]) / nonlinear['temperature2'][index]
            plt.plot(nonlinear['time'], intensity_1_err, label="Intensity 1")
            plt.plot(nonlinear['time'], intensity_2_err, label="Intensity 2")
            plt.plot(nonlinear['time'], temp_1_err, label="Temperature 1")
            plt.plot(nonlinear['time'], temp_2_err, label="Temperature 2")
            plt.title("Relative Error - Realizations vs. Heuristic")
            plt.xlabel("Time - ct (cm)")
            plt.ylabel("Relative Error")
            plt.xscale("log")
            plt.yscale("log")
            plt.grid(b=True, which="both", axis="both")
            plt.legend(loc="best")
            plt.tight_layout()
            plt.savefig(f"{PLOTPATH}/heur_real_rel_error.png")
            plt.cla()
            plt.clf()

            print("Realizations vs. Heuristic:")
            print(f"Intensity 1 - Mean Relative Error = {np.mean(intensity_1_err)}")
            print(f"Intensity 1 - Max Relative Error = {np.max(intensity_1_err)}")
            print(f"Intensity 2 - Mean Relative Error = {np.mean(intensity_2_err)}")
            print(f"Intensity 2 - Max Relative Error = {np.max(intensity_2_err)}")
            print(f"Temperature 1 - Mean Relative Error = {np.mean(temp_1_err)}")
            print(f"Temperature 1 - Max Relative Error = {np.max(temp_1_err)}")
            print(f"Temperature 2 - Mean Relative Error = {np.mean(temp_2_err)}")
            print(f"Temperature 2 - Max Relative Error = {np.max(temp_2_err)}")
            print("\n")

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

    # MONTE CARLO VS HEURISTIC MODEL
    if exists and heur_exists:
        # Intensity
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity1'], color='b', label="Material 1 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity2'], color='r', label="Material 2 - Monte Carlo")
        plt.plot(heur_data['time'], heur_data['intensity1'], color='c', linestyle=':', label="Material 1 - Heuristic")
        plt.plot(heur_data['time'], heur_data['intensity2'], color='m', linestyle=':', label="Material 2 - Heuristic")
        plt.title("Intensity - Monte Carlo vs. Heuristic")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intenisty (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/heur_mc_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature1'], color='b', label="Material 1 - Monte Carlo")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature2'], color='r', label="Material 2 - Monte Carlo")
        plt.plot(heur_data['time'], heur_data['temperature1'], color='c', linestyle=':', label="Material 1 - Heuristic")
        plt.plot(heur_data['time'], heur_data['temperature2'], color='m', linestyle=':', label="Material 2 - Heuristic")
        plt.title("Temperature - Monte Carlo vs. Heuristic")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/heur_mc_temperature.png")
        plt.cla()
        plt.clf()

        # Compute relative error
        if len(nonlinear_mc['time']) == len(heur_data['time']):
            intensity_1_err = np.zeros(len(nonlinear_mc['time']))
            intensity_2_err = np.zeros(len(nonlinear_mc['time']))
            temp_1_err = np.zeros(len(nonlinear_mc['time']))
            temp_2_err = np.zeros(len(nonlinear_mc['time']))
            for index, _ in enumerate(nonlinear_mc['time']):
                intensity_1_err[index] = np.abs(nonlinear_mc['intensity1'][index] - heur_data['intensity1'][index]) / nonlinear_mc['intensity1'][index]
                intensity_2_err[index] = np.abs(nonlinear_mc['intensity2'][index] - heur_data['intensity2'][index]) / nonlinear_mc['intensity1'][index]
                temp_1_err[index] = np.abs(nonlinear_mc['temperature1'][index] - heur_data['temperature1'][index]) / nonlinear_mc['temperature1'][index]
                temp_2_err[index] = np.abs(nonlinear_mc['temperature2'][index] - heur_data['temperature2'][index]) / nonlinear_mc['temperature2'][index]
            plt.plot(nonlinear['time'], intensity_1_err, label="Intensity 1")
            plt.plot(nonlinear['time'], intensity_2_err, label="Intensity 2")
            plt.plot(nonlinear['time'], temp_1_err, label="Temperature 1")
            plt.plot(nonlinear['time'], temp_2_err, label="Temperature 2")
            plt.title("Relative Error - Monte Carlo vs. Heuristic")
            plt.xlabel("Time - ct (cm)")
            plt.ylabel("Relative Error")
            plt.xscale("log")
            plt.yscale("log")
            plt.grid(b=True, which="both", axis="both")
            plt.legend(loc="best")
            plt.tight_layout()
            plt.savefig(f"{PLOTPATH}/heur_mc_rel_error.png")
            plt.cla()
            plt.clf()

            print("Monte Carlo vs. Heuristic:")
            print(f"Intensity 1 - Mean Relative Error = {np.mean(intensity_1_err)}")
            print(f"Intensity 1 - Max Relative Error = {np.max(intensity_1_err)}")
            print(f"Intensity 2 - Mean Relative Error = {np.mean(intensity_2_err)}")
            print(f"Intensity 2 - Max Relative Error = {np.max(intensity_2_err)}")
            print(f"Temperature 1 - Mean Relative Error = {np.mean(temp_1_err)}")
            print(f"Temperature 1 - Max Relative Error = {np.max(temp_1_err)}")
            print(f"Temperature 2 - Mean Relative Error = {np.mean(temp_2_err)}")
            print(f"Temperature 2 - Max Relative Error = {np.max(temp_2_err)}")
            print("\n")

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

    # REALIZATIONS VS MC COMPARISON
    if exists and mc_exists:
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

if __name__ == '__main__':
    main()
