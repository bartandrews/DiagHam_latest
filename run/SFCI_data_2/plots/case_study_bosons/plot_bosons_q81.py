import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import csv
from matplotlib.patches import Polygon
from math import gcd
from fractions import Fraction
import os
from matplotlib.patches import ConnectionPatch

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def can_convert_to_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


if __name__ == "__main__":

    # construct many-body spectrum 1
    mb_file = "bosons_hofstadter_X_9_Y_9_q_1_n_10_x_4_y_5_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    lin_K = []
    mb_E = []
    with open(os.path.join("q81/n10/nonauto", mb_file), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K.append(float(row[0])*5 + float(row[1]))
                mb_E.append(float(row[2]))

    # construct many-body spectrum 2
    mb_file2 = "bosons_hofstadter_X_9_Y_9_q_1_n_10_x_4_y_5_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    lin_K2 = []
    mb_E2 = []
    with open(os.path.join("q81/n10/nonauto", mb_file2), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K2.append(float(row[0])*5 + float(row[1]))
                mb_E2.append(float(row[2]))

    # construct entanglement spectrum 1
    ent_file = "bosons_hofstadter_X_9_Y_9_q_1_n_10_x_4_y_5_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_5.parentspec"
    lin_K_ent = []
    ent = []
    with open(os.path.join("q81/n10/ent2", ent_file), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K_ent.append(float(row[3]))
                ent.append(2*float(row[5]))
    comb = np.array([[lin_K_ent[i], ent[i]] for i in range(len(lin_K_ent))])
    comb = sorted(comb, key=lambda x: x[1])
    sorted_lin_K_ent = [comb[i][0] for i in range(len(lin_K_ent))]
    sorted_ent = [comb[i][1] for i in range(len(lin_K_ent))]
    lin_K_ent_below = sorted_lin_K_ent[:4004]
    lin_K_ent_above = sorted_lin_K_ent[4004:]
    ent_below = sorted_ent[:4004]
    ent_above = sorted_ent[4004:]

    # construct entanglement spectrum 2
    ent_file2 = "bosons_hofstadter_X_9_Y_9_q_1_n_10_x_4_y_5_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_5.parentspec"
    lin_K_ent2 = []
    ent2 = []
    with open(os.path.join(f"q81/n10/ent2", ent_file2), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K_ent2.append(float(row[3]))
                ent2.append(2*float(row[5]))
    comb2 = np.array([[lin_K_ent2[i], ent2[i]] for i in range(len(lin_K_ent2))])
    comb2 = sorted(comb2, key=lambda x: x[1])
    sorted_lin_K_ent2 = [comb2[i][0] for i in range(len(lin_K_ent2))]
    sorted_ent2 = [comb2[i][1] for i in range(len(lin_K_ent2))]
    lin_K_ent2_below = sorted_lin_K_ent2[:4004]
    lin_K_ent2_above = sorted_lin_K_ent2[4004:]
    ent2_below = sorted_ent2[:4004]
    ent2_above = sorted_ent2[4004:]

    ##########
    # figure #
    ##########

    fig = plt.figure(figsize=(6, 6))
    outer_grid = gridspec.GridSpec(2, 1, hspace=0.3, height_ratios=[3, 1])

    upper_cell = outer_grid[0]
    lower_cell = outer_grid[1]

    upper_inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2, upper_cell, hspace=0, wspace=0.3)
    lower_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, lower_cell, hspace=0, wspace=0.55)

    #####################
    # many-body spectra #
    #####################

    ax0 = plt.subplot(upper_inner_grid[0])
    ax0.set_title("$\mathrm{(a)}$ $\\langle \\mathcal{T} \\rangle = 0.0140$")
    scaled_mb_E = np.subtract(mb_E, min(mb_E))*1e3
    ax0.plot(lin_K, scaled_mb_E, '+', c='k', markersize=5)
    ax0.set_ylabel('$(E_\\mathrm{m.b.} - E_{\\mathrm{m.b.},0})/10^{-3}$')
    ax0.set_xlabel('$k_x L_y + k_y$')
    ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.annotate(s='', xy=(2, sorted(scaled_mb_E)[2]), xytext=(2, -0.1), arrowprops=dict(arrowstyle='<->'))
    ax0.text(4, (sorted(scaled_mb_E)[2]/2)*0.9, "$q\\Delta_\\mathrm{m.b.}=0.604$")
    ax0.xaxis.set_visible(False)

    ax1 = plt.subplot(upper_inner_grid[1])
    ax1.set_title("$\mathrm{(b)}$ $\\langle \\mathcal{T} \\rangle = 60.5$")
    scaled_mb_E2 = np.subtract(mb_E2, min(mb_E2))*1e3
    ax1.plot(lin_K2, scaled_mb_E2, '+', c='k', markersize=5)
    ax1.set_ylabel('$(E_\\mathrm{m.b.} - E_{\\mathrm{m.b.},0})/10^{-3}$')
    ax1.set_xlabel('$k_x L_y + k_y$')
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.annotate(s='', xy=(2, sorted(scaled_mb_E2)[2]), xytext=(2, -0.1), arrowprops=dict(arrowstyle='<->'))
    ax1.text(4, (sorted(scaled_mb_E2)[2]/2)*0.9, "$q \\Delta_\\mathrm{m.b.} = 0.356$")
    ax1.xaxis.set_visible(False)

    ########################
    # entanglement spectra #
    ########################

    ax2 = plt.subplot(upper_inner_grid[2])
    ax2.plot(lin_K_ent_below, ent_below, '+', c='r', markersize=2)
    ax2.plot(lin_K_ent_above, ent_above, '+', c='k', markersize=2)
    ax2.set_ylabel('$\\xi$')
    ax2.set_xlabel('$k_x L_y + k_y$')
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.annotate(s='', xy=(2, ent_above[0]), xytext=(2, ent_below[-1]), arrowprops=dict(arrowstyle='<->'))
    ax2.text(4, (ent_below[-1]+(ent_above[0] - ent_below[-1])/2)*0.9, "$\\Delta_\\xi=22.1$")

    ax3 = plt.subplot(upper_inner_grid[3])
    ax3.plot(lin_K_ent2_below, ent2_below, '+', c='r', markersize=2)
    ax3.plot(lin_K_ent2_above, ent2_above, '+', c='k', markersize=2)
    ax3.set_ylabel('$\\xi$')
    ax3.set_xlabel('$k_x L_y + k_y$')
    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # ax3.annotate(s='', xy=(2, 2 * 23.05), xytext=(2, 2 * 12.5), arrowprops=dict(arrowstyle='<->'))
    # ax3.text(4, 2 * 17, "$\\Delta_\\xi=22.0$")

    #######################
    # finite-size scaling #
    #######################

    # extract data
    fs_file = f"data_bosons_q81.txt"
    invertN = []
    scaled_gap1, scaled_gap2 = [], []
    scaled_ent1, scaled_ent2 = [], []
    with open(fs_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                invertN.append(1/float(row[0]))
                scaled_gap1.append(float(row[1]))
                scaled_gap2.append(float(row[2]))
                scaled_ent1.append(float(row[3]))
                scaled_ent2.append(float(row[4]))

    color1 = 'C0'
    ax4 = plt.subplot(lower_inner_grid[0])
    ax4.plot(invertN, scaled_gap1, '.-', color=color1)
    ax4.set_xlabel('$1/N$')
    ax4.set_ylabel('$q \\Delta_\\mathrm{m.b.}$', color=color1)
    ax4.tick_params(axis='y', labelcolor=color1)
    ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.axhline(0.616, ls='--')
    ax4.set_ylim([0, 0.8])
    ax4.axvline(0.1, ls='-', zorder=-5, c='gray')
    #
    color2 = 'C3'
    ax4_2 = ax4.twinx()
    ax4_2.plot(invertN, scaled_ent1, '.-', zorder=5, color=color2, marker='x')
    ax4_2.set_ylabel(f'$\\Delta_\\xi$', color=color2)
    ax4_2.tick_params(axis='y', labelcolor=color2)
    ax4_2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4_2.set_ylim([0, 30])

    ax5 = plt.subplot(lower_inner_grid[1])
    ax5.plot(invertN, scaled_gap2, '.-', color=color1)
    ax5.set_xlabel('$1/N$')
    ax5.set_ylabel('$q \\Delta_\\mathrm{m.b.}$', color=color1)
    ax5.tick_params(axis='y', labelcolor=color1)
    ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5.axhline(0.616, ls='--')
    ax5.set_ylim([0, 0.8])
    ax5.axvline(0.1, ls='-', zorder=-5, c='gray')
    #
    ax5_2 = ax5.twinx()
    ax5_2.plot(invertN, scaled_ent2, '.-', zorder=5, color=color2, marker='x')
    ax5_2.set_ylabel(f'$\\Delta_\\xi$', color=color2)
    ax5_2.tick_params(axis='y', labelcolor=color2)
    ax5_2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5_2.set_ylim([0, 30])

    ###

    # (-1, 12)
    left_con = ConnectionPatch(xyA=(0.1, 0.8), xyB=(0, 0), coordsA="data", coordsB="axes points",
                               axesA=ax4, axesB=ax2, connectionstyle="angle3,angleA=150,angleB=235", arrowstyle='->',
                               facecolor='k', edgecolor='k', zorder=2)
    fig.add_artist(left_con)
    # (-1, 10.5)
    right_con = ConnectionPatch(xyA=(0.1, 0.8), xyB=(0, 0), coordsA="data", coordsB="axes points",
                                axesA=ax5, axesB=ax3, connectionstyle="angle3,angleA=150,angleB=235", arrowstyle='->',
                                facecolor='k', edgecolor='k', zorder=3)
    fig.add_artist(right_con)

    ###

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/case_study_bosons/bosons_q81_3.png",
                bbox_inches='tight', dpi=300)
    plt.show()
