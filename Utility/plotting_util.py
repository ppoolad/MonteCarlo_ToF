"""
A collection of utilities to help with plotting
Feel free to add utilities as needed.

Author: Dylan Vogel
Last Modified: 2021-04-23
"""

# ============================================================================
# IMPORTS
# ============================================================================

import matplotlib.pyplot as plt
import numpy as np


# ============================================================================
# FUNCTIONS
# ============================================================================

def plot_histogram(hist, period, nbins, num_data, outname, xlabel=False):
    # merged_sipm = np.sum(SiPM_hists_a, axis=0)
    #print(np.shape(merged_sipm))

    # get a random figure number
    plt.figure(np.random.randint(0, 1000), figsize=(7, 0.5 + num_data * 1.5), dpi=300)
    period_ns = period * 1e9  # convert to ns for convenience
    hst_time = np.linspace(0, period_ns * (1-1/nbins), num=nbins)

    for i in range(num_data):
        plt.subplot(num_data, 1, i + 1)

        plt.rcParams['axes.axisbelow'] = True
        plt.grid(b=True, which='major', alpha=0.5)
        plt.grid(b=True, which='minor', linestyle='-', alpha=0.1)
        plt.minorticks_on()
        # just to make it easier to process
        hst_values = np.sum(np.reshape(hist[i], (-1, nbins), order='F'), axis=0)
        plt.bar(hst_time, hst_values, width=period_ns/nbins)  # 1?
        # plt.bar(time_steps, SiPM_hists_a[i], width=const.sim_res)
        plt.xlim((0, period_ns))

    # plt.xticks(np.arange(0, 100e-9, 10e-9))
    if xlabel:
        # add marign and text for the x-axis
        plt.subplots_adjust(bottom=0.25)
        plt.xlabel('Time [ns]')
    else:
        plt.subplots_adjust(bottom=0.15)  # just enough space so labels fit
    plt.ylabel('Counts')

    plt.savefig(outname)

    plt.show()


