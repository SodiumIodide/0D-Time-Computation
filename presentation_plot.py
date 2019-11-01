#!/usr/bin/env python3
'''
Read data and create figures
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

HISTPATH = "out/nonlinear/pdf_data"
DATAPATH = "out/nonlinear/data"
PLOTPATH = "out/nonlinear"

XSCALE = "log"
YSCALE = "linear"
HISTXSCALE = "linear"
HISTYSCALE = "linear"

def main():
    '''Main function'''
    serial_exists = gpu_exists = hist_exists = True
    try:
        serial_data = pd.read_csv(f"./{DATAPATH}/output.csv")
    except FileNotFoundError:
        serial_exists = False
    try:
        gpu_data = pd.read_csv(f"./{DATAPATH}/gpu_output.csv")
    except FileNotFoundError:
        gpu_exists = False
    try:
        hist_data = pd.read_csv(f"./{HISTPATH}/mc_pdf.csv")
    except FileNotFoundError:
        hist_exists = False

    if not os.path.exists(f"./{PLOTPATH}"):
        os.makedirs(f"./{PLOTPATH}")

    if serial_exists:
        # Standard deviation computing
        sr_std_intensity_1 = serial_data['varintensity1'].apply(np.abs).apply(np.sqrt)
        sr_std_intensity_2 = serial_data['varintensity2'].apply(np.abs).apply(np.sqrt)
        sr_std_temp_1 = serial_data['vartemperature1'].apply(np.abs).apply(np.sqrt)
        sr_std_temp_2 = serial_data['vartemperature2'].apply(np.abs).apply(np.sqrt)
        sr_std_opac_1 = serial_data['varopacity1'].apply(np.abs).apply(np.sqrt)
        sr_std_opac_2 = serial_data['varopacity2'].apply(np.abs).apply(np.sqrt)
        sr_lb_intensity_1 = serial_data['intensity1'] - sr_std_intensity_1
        sr_ub_intensity_1 = serial_data['intensity1'] + sr_std_intensity_1
        sr_lb_intensity_2 = serial_data['intensity2'] - sr_std_intensity_2
        sr_ub_intensity_2 = serial_data['intensity2'] + sr_std_intensity_2
        sr_lb_temp_1 = serial_data['temperature1'] - sr_std_temp_1
        sr_ub_temp_1 = serial_data['temperature1'] + sr_std_temp_1
        sr_lb_temp_2 = serial_data['temperature2'] - sr_std_temp_2
        sr_ub_temp_2 = serial_data['temperature2'] + sr_std_temp_2
        sr_lb_opac_1 = serial_data['opacity1'] - sr_std_opac_1
        sr_ub_opac_1 = serial_data['opacity1'] + sr_std_opac_1
        sr_lb_opac_2 = serial_data['opacity2'] - sr_std_opac_2
        sr_ub_opac_2 = serial_data['opacity2'] + sr_std_opac_2

        # Intensity
        plt.plot(serial_data['time'], serial_data['intensity1'], color='b', label="Material 1")
        plt.plot(serial_data['time'], serial_data['intensity2'], color='r', label="Material 2")
        plt.plot(serial_data['time'], sr_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(serial_data['time'], sr_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(serial_data['time'], sr_lb_intensity_1, sr_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(serial_data['time'], sr_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(serial_data['time'], sr_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(serial_data['time'], sr_lb_intensity_2, sr_ub_intensity_2, color='r', alpha=0.5)
        plt.title("Serial Intensity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm$^2$-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/serial_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(serial_data['time'], serial_data['temperature1'], color='b', label="Material 1")
        plt.plot(serial_data['time'], serial_data['temperature2'], color='r', label="Material 2")
        plt.plot(serial_data['time'], sr_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(serial_data['time'], sr_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(serial_data['time'], sr_lb_temp_1, sr_ub_temp_1, color='b', alpha=0.5)
        plt.plot(serial_data['time'], sr_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(serial_data['time'], sr_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(serial_data['time'], sr_lb_temp_2, sr_ub_temp_2, color='r', alpha=0.5)
        plt.title("Serial Temperature")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/serial_temperature.png")
        plt.cla()
        plt.clf()

        # Opacity
        plt.plot(serial_data['time'], serial_data['opacity1'], color='b', label="Material 1")
        plt.plot(serial_data['time'], serial_data['opacity2'], color='r', label="Material 2")
        plt.plot(serial_data['time'], sr_lb_opac_1, color='b', linestyle=':', label=None)
        plt.plot(serial_data['time'], sr_ub_opac_1, color='b', linestyle=':', label=None)
        plt.fill_between(serial_data['time'], sr_lb_opac_1, sr_ub_opac_1, color='b', alpha=0.5)
        plt.plot(serial_data['time'], sr_lb_opac_2, color='r', linestyle=':', label=None)
        plt.plot(serial_data['time'], sr_ub_opac_2, color='r', linestyle=':', label=None)
        plt.fill_between(serial_data['time'], sr_lb_opac_2, sr_ub_opac_2, color='r', alpha=0.5)
        plt.title("Serial Opacity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Opacity (cm$^{-1}$)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/serial_opacity.png")
        plt.cla()
        plt.clf()

    if gpu_exists:
        # Standard deviation computing
        gp_std_intensity_1 = gpu_data['varintensity1'].apply(np.abs).apply(np.sqrt)
        gp_std_intensity_2 = gpu_data['varintensity2'].apply(np.abs).apply(np.sqrt)
        gp_std_temp_1 = gpu_data['vartemperature1'].apply(np.abs).apply(np.sqrt)
        gp_std_temp_2 = gpu_data['vartemperature2'].apply(np.abs).apply(np.sqrt)
        gp_std_opac_1 = gpu_data['varopacity1'].apply(np.abs).apply(np.sqrt)
        gp_std_opac_2 = gpu_data['varopacity2'].apply(np.abs).apply(np.sqrt)
        gp_lb_intensity_1 = gpu_data['intensity1'] - gp_std_intensity_1
        gp_ub_intensity_1 = gpu_data['intensity1'] + gp_std_intensity_1
        gp_lb_intensity_2 = gpu_data['intensity2'] - gp_std_intensity_2
        gp_ub_intensity_2 = gpu_data['intensity2'] + gp_std_intensity_2
        gp_lb_temp_1 = gpu_data['temperature1'] - gp_std_temp_1
        gp_ub_temp_1 = gpu_data['temperature1'] + gp_std_temp_1
        gp_lb_temp_2 = gpu_data['temperature2'] - gp_std_temp_2
        gp_ub_temp_2 = gpu_data['temperature2'] + gp_std_temp_2
        gp_lb_opac_1 = gpu_data['opacity1'] - gp_std_opac_1
        gp_ub_opac_1 = gpu_data['opacity1'] + gp_std_opac_1
        gp_lb_opac_2 = gpu_data['opacity2'] - gp_std_opac_2
        gp_ub_opac_2 = gpu_data['opacity2'] + gp_std_opac_2

        # Intensity
        plt.plot(gpu_data['time'], gpu_data['intensity1'], color='b', label="Material 1")
        plt.plot(gpu_data['time'], gpu_data['intensity2'], color='r', label="Material 2")
        plt.plot(gpu_data['time'], gp_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(gpu_data['time'], gp_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(gpu_data['time'], gp_lb_intensity_1, gp_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(gpu_data['time'], gp_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(gpu_data['time'], gp_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(gpu_data['time'], gp_lb_intensity_2, gp_ub_intensity_2, color='r', alpha=0.5)
        #plt.title("GPU Intensity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm$^2$-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/gpu_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(gpu_data['time'], gpu_data['temperature1'], color='b', label="Material 1")
        plt.plot(gpu_data['time'], gpu_data['temperature2'], color='r', label="Material 2")
        plt.plot(gpu_data['time'], gp_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(gpu_data['time'], gp_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(gpu_data['time'], gp_lb_temp_1, gp_ub_temp_1, color='b', alpha=0.5)
        plt.plot(gpu_data['time'], gp_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(gpu_data['time'], gp_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(gpu_data['time'], gp_lb_temp_2, gp_ub_temp_2, color='r', alpha=0.5)
        #plt.title("GPU Temperature")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/gpu_temperature.png")
        plt.cla()
        plt.clf()

        # Opacity
        plt.plot(gpu_data['time'], gpu_data['opacity1'], color='b', label="Material 1")
        plt.plot(gpu_data['time'], gpu_data['opacity2'], color='r', label="Material 2")
        plt.plot(gpu_data['time'], gp_lb_opac_1, color='b', linestyle=':', label=None)
        plt.plot(gpu_data['time'], gp_ub_opac_1, color='b', linestyle=':', label=None)
        plt.fill_between(gpu_data['time'], gp_lb_opac_1, gp_ub_opac_1, color='b', alpha=0.5)
        plt.plot(gpu_data['time'], gp_lb_opac_2, color='r', linestyle=':', label=None)
        plt.plot(gpu_data['time'], gp_ub_opac_2, color='r', linestyle=':', label=None)
        plt.fill_between(gpu_data['time'], gp_lb_opac_2, gp_ub_opac_2, color='r', alpha=0.5)
        #plt.title("GPU Opacity")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Opacity (cm$^{-1}$)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/gpu_opacity.png")
        plt.cla()
        plt.clf()

    if hist_exists:
        time_ss = hist_data['time'][0]

        # Intensity Material 1 and 2
        plt.plot(hist_data['intensity1arr'][1:len(hist_data['intensity1arr'])-1], hist_data['freqintensity1'][1:len(hist_data['freqintensity1'])-1], color='b', label="Material 1")
        plt.plot(hist_data['intensity2arr'][1:len(hist_data['intensity2arr'])-1], hist_data['freqintensity2'][1:len(hist_data['freqintensity2'])-1], color='r', label="Material 2")
        plt.axvline(x=21848239083.052482, color='k', label="Atomic Mix")
        plt.axvline(x=2.115966666e10, color='c', label="Ensemble Average")
        plt.axvline(x=20484661056.158676, color='g', label="Heuristic")
        plt.title(f"Intensity PDF - GPU Monte Carlo (Time = {time_ss:.4E} ct)")
        plt.xlabel("Intensity (erg/cm$^2$-s)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/intensity_hist.png")
        plt.cla()
        plt.clf()

        # Temperature Material 1 and 2
        plt.plot(hist_data['temperature1arr'][1:len(hist_data['temperature1arr'])-1], hist_data['freqtemperature1'][1:len(hist_data['freqtemperature1'])-1], color='b', label="Material 1")
        plt.plot(hist_data['temperature2arr'][1:len(hist_data['temperature2arr'])-1], hist_data['freqtemperature2'][1:len(hist_data['freqtemperature2'])-1], color='r', label="Material 2")
        plt.axvline(x=0.2712211898921925, color='k', label="Atomic Mix")
        plt.axvline(x=0.2906166666, color='c', label="Ensemble Average")
        plt.axvline(x=0.31670525696951857, color='g', label="Heuristic")
        plt.title(f"Temperature PDF - GPU Monte Carlo (Time = {time_ss:.4E} ct)")
        plt.xlabel("Temperature (eV)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/temperature_hist.png")
        plt.cla()
        plt.clf()

        # Opacity Material 1 and 2
        plt.plot(hist_data['opacity1arr'][1:len(hist_data['opacity1arr'])-1], hist_data['freqopacity1'][1:len(hist_data['freqopacity1'])-1], color='b', label="Material 1")
        plt.plot(hist_data['opacity2arr'][1:len(hist_data['opacity2arr'])-1], hist_data['freqopacity2'][1:len(hist_data['freqopacity2'])-1], color='r', label="Material 2")
        plt.axvline(x=217.17419755097302, color='k', label="Atomic Mix")
        plt.axvline(x=188.9405, color='c', label="Ensemble Average")
        plt.axvline(x=181.94386078272075, color='g', label="Heuristic")
        plt.title(f"Opacity PDF - GPU Monte Carlo (Time = {time_ss:.4E} ct)")
        plt.xlabel("Opacity (cm$^{-1}$)")
        plt.ylabel("Probability Density")
        plt.xscale(HISTXSCALE)
        plt.yscale(HISTYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"./{PLOTPATH}/opacity_hist.png")
        plt.cla()
        plt.clf()

        print("Output:")
        average_temp_1 = np.sum([hist_data['temperature1arr'].iat[i] * hist_data['freqtemperature1'].iat[i] for i in range(hist_data.shape[0])])
        average_temp_2 = np.sum([hist_data['temperature2arr'].iat[i] * hist_data['freqtemperature2'].iat[i] for i in range(hist_data.shape[0])])
        average_opac_1 = np.sum([hist_data['opacity1arr'].iat[i] * hist_data['freqopacity1'].iat[i] for i in range(hist_data.shape[0])])
        average_opac_2 = np.sum([hist_data['opacity2arr'].iat[i] * hist_data['freqopacity2'].iat[i] for i in range(hist_data.shape[0])])
        average_int_1 = np.sum([hist_data['intensity1arr'].iat[i] * hist_data['freqintensity1'].iat[i] for i in range(hist_data.shape[0])])
        average_int_2 = np.sum([hist_data['intensity2arr'].iat[i] * hist_data['freqintensity2'].iat[i] for i in range(hist_data.shape[0])])
        print(f"Average Material 1 Temperature: {average_temp_1:.4f}")
        print(f"Average Material 2 Temperature: {average_temp_2:.4f}")
        print(f"Average Material 1 Opacity: {average_opac_1:.4f}")
        print(f"Average Material 2 Opacity: {average_opac_2:.4f}")
        print(f"Average Material 1 Intensity: {average_int_1:.4e}")
        print(f"Average Material 2 Intensity: {average_int_2:.4e}")


if __name__ == '__main__':
    print("Plotting...")
    main()
    print("Done")
