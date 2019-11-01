#!/usr/bin/env python3
'''
Produce histogram figures for the ANS report from MC data
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLOTPATH = "./out/nonlinear"
HISTPATH = f"{PLOTPATH}/pdf_data"
DATAPATH = f"{PLOTPATH}/data"

PLOT_STD_DEV = False

XSCALE = "log"
YSCALE = "linear"

HISTXSCALE = "linear"
HISTYSCALE = "linear"

def main():
    '''Main function'''
    data_exists = num_exists = mc_exists = am_exists = lp_exists = True

    #time,intensity1arr,freqintensity1,intensity2arr,freqintensity2,temperature1arr,freqtemperature1,temperature2arr,freqtemperature2,opacity1arr,freqopacity1,opacity2arr,freqopacity2
    try:
        data = pd.read_csv(f"{HISTPATH}/mc_pdf.csv")
    except FileNotFoundError:
        data_exists = False
    try:
        nonlinear = pd.read_csv(f"{DATAPATH}/nonlinear.csv")
    except FileNotFoundError:
        num_exists = False
    try:
        nonlinear_mc = pd.read_csv(f"{DATAPATH}/nonlinearmc.csv")
    except FileNotFoundError:
        mc_exists = False
    try:
        am = pd.read_csv(f"{DATAPATH}/atomicmix.csv")
    except FileNotFoundError:
        am_exists = False
    try:
        lp = pd.read_csv(f"{DATAPATH}/levermorepomraning.csv")
    except FileNotFoundError:
        lp_exists = False

    if data_exists:
        #time_ss = data['time'][0]

        # Intensity Material 1 and 2
        plt.plot(data['intensity1arr'][1:len(data['intensity1arr'])-1], data['freqintensity1'][1:len(data['freqintensity1'])-1], color='b', label="Material 1")
        plt.plot(data['intensity2arr'][1:len(data['intensity2arr'])-1], data['freqintensity2'][1:len(data['freqintensity2'])-1], color='r', label="Material 2")
        plt.xlabel("Intensity (ergs/cm$^2$-s)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/intensity_mc_hist.png")
        plt.cla()
        plt.clf()

        # Temperature Material 1 and 2
        plt.plot(data['temperature1arr'][1:len(data['temperature1arr'])-1], data['freqtemperature1'][1:len(data['freqtemperature1'])-1], color='b', label="Material 1")
        plt.plot(data['temperature2arr'][1:len(data['temperature2arr'])-1], data['freqtemperature2'][1:len(data['freqtemperature2'])-1], color='r', label="Material 2")
        plt.xlabel("Temperature (eV)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_mc_hist.png")
        plt.cla()
        plt.clf()

        # Opacity Material 1 and 2
        plt.plot(data['opacity1arr'][1:len(data['opacity1arr'])-1], data['freqopacity1'][1:len(data['freqopacity1'])-1], color='b', label="Material 1")
        plt.plot(data['opacity2arr'][1:len(data['opacity2arr'])-1], data['freqopacity2'][1:len(data['freqopacity2'])-1], color='r', label="Material 2")
        plt.xlabel("Opacity (cm$^{-1}$)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/opacity_mc_hist.png")
        plt.cla()
        plt.clf()

    if data_exists and am_exists and lp_exists:
        #time_ss = data['time'][0]
        #ts_df = am.query(f"time == {time_ss}")
        ts_df = am.query("time == 0.0037024368563")
        ts_lp = lp.query("time == 0.0037024368563")
        am_intensity, am_temperature, am_opacity = (ts_df['intensityarr'].iat[0], ts_df['temperaturearr'].iat[0], ts_df['opacityarr'].iat[0])
        lp_intensity_1, lp_intensity_2 = (ts_lp['intensity1'].iat[0], ts_lp['intensity2'].iat[0])
        lp_temperature_1, lp_temperature_2 = (ts_lp['temperature1'].iat[0], ts_lp['temperature2'].iat[0])
        lp_opacity_1, lp_opacity_2 = (ts_lp['opacity1'].iat[0], ts_lp['opacity2'].iat[0])
        lp_intensity, lp_temperature, lp_opacity = (ts_lp['intensityarr'].iat[0], ts_lp['temperaturearr'].iat[0], ts_lp['opacityarr'].iat[0])
        print(f"AM Intensity: {am_intensity}")
        print(f"AM Temperature: {am_temperature}")
        print(f"AM Opacity: {am_opacity}")
        print(f"Heuristic Intensity: {lp_intensity}")
        print(f"Heuristic Temperature: {lp_temperature}")
        print(f"Heuristic Opacity: {lp_opacity}")

        # Intensity Material 1 and 2
        plt.plot(data['intensity1arr'][1:len(data['intensity1arr'])-1], data['freqintensity1'][1:len(data['freqintensity1'])-1], color='b', label="Material 1")
        plt.plot(data['intensity2arr'][1:len(data['intensity2arr'])-1], data['freqintensity2'][1:len(data['freqintensity2'])-1], color='r', label="Material 2")
        plt.axvline(x=am_intensity, color='k', label="Atomic Mix")
        plt.xlabel("Intensity (ergs/cm$^2$-s)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/am_intensity_mc_hist.png")
        plt.cla()
        plt.clf()

        # Temperature Material 1 and 2
        plt.plot(data['temperature1arr'][1:len(data['temperature1arr'])-1], data['freqtemperature1'][1:len(data['freqtemperature1'])-1], color='b', label="Material 1")
        plt.plot(data['temperature2arr'][1:len(data['temperature2arr'])-1], data['freqtemperature2'][1:len(data['freqtemperature2'])-1], color='r', label="Material 2")
        plt.axvline(x=am_temperature, color='k', label="Atomic Mix")
        plt.xlabel("Temperature (eV)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/am_temperature_mc_hist.png")
        plt.cla()
        plt.clf()

        # Opacity Material 1 and 2
        plt.plot(data['opacity1arr'][1:len(data['opacity1arr'])-1], data['freqopacity1'][1:len(data['freqopacity1'])-1], color='b', label="Material 1")
        plt.plot(data['opacity2arr'][1:len(data['opacity2arr'])-1], data['freqopacity2'][1:len(data['freqopacity2'])-1], color='r', label="Material 2")
        plt.axvline(x=am_opacity, color='k', label="Atomic Mix")
        plt.xlabel("Opacity (cm$^{-1}$)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/am_opacity_mc_hist.png")
        plt.cla()
        plt.clf()

    if num_exists and mc_exists and am_exists and lp_exists:
        # Standard deviation computing
        nl_std_intensity_1 = nonlinear['varintensity1'].apply(np.abs).apply(np.sqrt)
        nl_std_intensity_2 = nonlinear['varintensity2'].apply(np.abs).apply(np.sqrt)
        nl_std_temp_1 = nonlinear['vartemperature1'].apply(np.abs).apply(np.sqrt)
        nl_std_temp_2 = nonlinear['vartemperature2'].apply(np.abs).apply(np.sqrt)
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
        if PLOT_STD_DEV:
            plt.plot(nonlinear['time'], nl_lb_intensity_1, color='b', linestyle=':', label=None)
            plt.plot(nonlinear['time'], nl_ub_intensity_1, color='b', linestyle=':', label=None)
            plt.fill_between(nonlinear['time'], nl_lb_intensity_1, nl_ub_intensity_1, color='b', alpha=0.5)
            plt.plot(nonlinear['time'], nl_lb_intensity_2, color='r', linestyle=':', label=None)
            plt.plot(nonlinear['time'], nl_ub_intensity_2, color='r', linestyle=':', label=None)
            plt.fill_between(nonlinear['time'], nl_lb_intensity_2, nl_ub_intensity_2, color='r', alpha=0.5)
        plt.plot(am['time'], am['intensityarr'], color='k', label="Atomic Mix")
        plt.plot(lp['time'], lp['intensityarr'], color='g', label="Heuristic")
        #plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm$^2$-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/am_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2")
        if PLOT_STD_DEV:
            plt.plot(nonlinear['time'], nl_lb_temp_1, color='b', linestyle=':', label=None)
            plt.plot(nonlinear['time'], nl_ub_temp_1, color='b', linestyle=':', label=None)
            plt.fill_between(nonlinear['time'], nl_lb_temp_1, nl_ub_temp_1, color='b', alpha=0.5)
            plt.plot(nonlinear['time'], nl_lb_temp_2, color='r', linestyle=':', label=None)
            plt.plot(nonlinear['time'], nl_ub_temp_2, color='r', linestyle=':', label=None)
            plt.fill_between(nonlinear['time'], nl_lb_temp_2, nl_ub_temp_2, color='r', alpha=0.5)
        plt.plot(am['time'], am['temperaturearr'], color='k', label="Atomic Mix")
        plt.plot(lp['time'], lp['temperaturearr'], color='g', label="Heuristic")
        #plt.title("Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/am_temperature.png")
        plt.cla()
        plt.clf()

        # Standard deviation computing
        mc_nl_std_intensity_1 = nonlinear_mc['varintensity1'].apply(np.abs).apply(np.sqrt)
        mc_nl_std_intensity_2 = nonlinear_mc['varintensity2'].apply(np.abs).apply(np.sqrt)
        mc_nl_std_temp_1 = nonlinear_mc['vartemperature1'].apply(np.abs).apply(np.sqrt)
        mc_nl_std_temp_2 = nonlinear_mc['vartemperature2'].apply(np.abs).apply(np.sqrt)
        mc_nl_lb_intensity_1 = nonlinear_mc['intensity1'] - mc_nl_std_intensity_1
        mc_nl_ub_intensity_1 = nonlinear_mc['intensity1'] + mc_nl_std_intensity_1
        mc_nl_lb_intensity_2 = nonlinear_mc['intensity2'] - mc_nl_std_intensity_2
        mc_nl_ub_intensity_2 = nonlinear_mc['intensity2'] + mc_nl_std_intensity_2
        mc_nl_lb_temp_1 = nonlinear_mc['temperature1'] - mc_nl_std_temp_1
        mc_nl_ub_temp_1 = nonlinear_mc['temperature1'] + mc_nl_std_temp_1
        mc_nl_lb_temp_2 = nonlinear_mc['temperature2'] - mc_nl_std_temp_2
        mc_nl_ub_temp_2 = nonlinear_mc['temperature2'] + mc_nl_std_temp_2

        # Intensity
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensityarr'], color='c', label="Exact")
        #plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity1'], color='b', label="Material 1")
        #plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity2'], color='r', label="Material 2")
        if PLOT_STD_DEV:
            plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_1, color='b', linestyle=':', label=None)
            plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_1, color='b', linestyle=':', label=None)
            plt.fill_between(nonlinear_mc['time'], mc_nl_lb_intensity_1, mc_nl_ub_intensity_1, color='b', alpha=0.5)
            plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_2, color='r', linestyle=':', label=None)
            plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_2, color='r', linestyle=':', label=None)
            plt.fill_between(nonlinear_mc['time'], mc_nl_lb_intensity_2, mc_nl_ub_intensity_2, color='r', alpha=0.5)
        plt.plot(am['time'], am['intensityarr'], color='k', label="Atomic Mix")
        plt.plot(lp['time'], lp['intensityarr'], color='g', label="Heuristic")
        #plt.title("Monte Carlo Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm$^2$-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/am_intensity_mc.png")
        plt.cla()
        plt.clf()

        # Temperature
        #plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature1'], color='b', label="Material 1")
        #plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature2'], color='r', label="Material 2")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperaturearr'], color='c', label="Exact")
        if PLOT_STD_DEV:
            plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_1, color='b', linestyle=':', label=None)
            plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_1, color='b', linestyle=':', label=None)
            plt.fill_between(nonlinear_mc['time'], mc_nl_lb_temp_1, mc_nl_ub_temp_1, color='b', alpha=0.5)
            plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_2, color='r', linestyle=':', label=None)
            plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_2, color='r', linestyle=':', label=None)
            plt.fill_between(nonlinear_mc['time'], mc_nl_lb_temp_2, mc_nl_ub_temp_2, color='r', alpha=0.5)
        plt.plot(am['time'], am['temperaturearr'], color='k', label="Atomic Mix")
        plt.plot(lp['time'], lp['temperaturearr'], color='g', label="Heuristic")
        #plt.title("Monte Carlo Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/am_temperature_mc.png")
        plt.cla()
        plt.clf()

if __name__ == '__main__':
    main()
