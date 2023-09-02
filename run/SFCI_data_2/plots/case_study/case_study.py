import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import csv
from matplotlib.patches import Polygon
from math import gcd
from fractions import Fraction

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def can_convert_to_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


if __name__ == "__main__":

    ent_factor = 2

    # construct many-body spectrum 1
    mb_file = "fermions_hofstadter_X_12_Y_10_q_1_n_10_x_5_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    lin_K = []
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K.append(float(row[0])*6+float(row[1]))
                mb_E.append(float(row[2])/2)

    # construct many-body spectrum 2
    mb_file2 = "fermions_hofstadter_X_12_Y_10_q_1_n_10_x_5_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    lin_K2 = []
    mb_E2 = []
    with open(mb_file2, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K2.append(float(row[0]) * 6 + float(row[1]))
                mb_E2.append(float(row[2]) / 2)

    # construct entanglement spectrum 1
    ent_file = "ent/fermions_hofstadter_X_12_Y_10_q_1_n_10_x_5_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_5.parentspec"
    lin_K_ent = []
    ent = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K_ent.append(float(row[3]))
                ent.append(ent_factor*float(row[5]))

    # construct entanglement spectrum 2
    ent_file2 = "ent/fermions_hofstadter_X_12_Y_10_q_1_n_10_x_5_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_5.parentspec"
    lin_K_ent2 = []
    ent2 = []
    with open(ent_file2, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K_ent2.append(float(row[3]))
                ent2.append(ent_factor*float(row[5]))

    ##########
    # figure #
    ##########

    fig = plt.figure(figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2, hspace=0, wspace=0.3)

    #####################
    # many-body spectra #
    #####################

    ax0 = plt.subplot(gs[0])
    ax0.set_title("(a) $\\langle \\mathcal{T} \\rangle = 0.0141$")
    ax0.plot(lin_K, np.subtract(mb_E, min(mb_E))*1e3, '+', c='k', markersize=5)
    ax0.set_ylabel('$(E_\\mathrm{m.b.} - E_{\\mathrm{m.b.},0})/10^{-3}$')
    ax0.set_xlabel('$k_x L_y + k_y$')
    ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.annotate(s='', xy=(2, 0.085), xytext=(2, 0), arrowprops=dict(arrowstyle='<->'))
    ax0.text(4, 0.035, "$q^2\\Delta_\\mathrm{m.b.}=1.26$")
    ax0.xaxis.set_visible(False)

    ax1 = plt.subplot(gs[1])
    ax1.set_title("(b) $\\langle \\mathcal{T} \\rangle = 72.1$")
    ax1.plot(lin_K2, np.subtract(mb_E2, min(mb_E2)) * 1e3, '+', c='k', markersize=5)
    ax1.set_ylabel('$(E_\\mathrm{m.b.} - E_{\\mathrm{m.b.},0})/10^{-3}$')
    ax1.set_xlabel('$k_x L_y + k_y$')
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.annotate(s='', xy=(2, 0.25), xytext=(2, 0), arrowprops=dict(arrowstyle='<->'))
    ax1.text(4, 0.12, "$q^2 \\Delta_\\mathrm{m.b.} = 3.62$")
    ax1.xaxis.set_visible(False)

    ########################
    # entanglement spectra #
    ########################

    ax2 = plt.subplot(gs[2])
    ax2.plot(lin_K_ent, ent, '+', c='k', markersize=2)
    ax2.set_ylabel('$\\xi$')
    ax2.set_xlabel('$k_x L_y + k_y$')
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.annotate(s='', xy=(2, ent_factor*19.5), xytext=(2, ent_factor*12.5), arrowprops=dict(arrowstyle='<->'))
    ax2.text(4, ent_factor*15, "$\\Delta_\\xi=14.6$")

    ax3 = plt.subplot(gs[3])
    ax3.plot(lin_K_ent2, ent2, '+', c='k', markersize=2)
    ax3.set_ylabel('$\\xi$')
    ax3.set_xlabel('$k_x L_y + k_y$')
    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.annotate(s='', xy=(2, ent_factor*23.05), xytext=(2, ent_factor*12.5), arrowprops=dict(arrowstyle='<->'))
    ax3.text(4, ent_factor*17, "$\\Delta_\\xi=22.0$")

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/case_study/case_study.png",
                bbox_inches='tight', dpi=600)
    plt.show()
