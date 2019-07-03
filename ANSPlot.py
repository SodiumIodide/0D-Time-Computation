#!/usr/bin/env python3
'''
Produce histogram figures for the ANS report from MC data
'''

import matplotlib.pyplot as plt
import pandas as pd

PLOTPATH = "./out/nonlinear"
HISTPATH = f"{PLOTPATH}/pdf_data"

HISTXSCALE = "linear"
HISTYSCALE = "linear"

def main():
    '''Main function'''
    data_exists = True

    #time,intensity1arr,freqintensity1,intensity2arr,freqintensity2,temperature1arr,freqtemperature1,temperature2arr,freqtemperature2,opacity1arr,freqopacity1,opacity2arr,freqopacity2
    try:
        data = pd.read_csv(f"{HISTPATH}/mc_pdf.csv")
    except FileNotFoundError:
        data_exists = False

    if data_exists:
        time_ss = data['time'][0]

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


if __name__ == '__main__':
    main()
