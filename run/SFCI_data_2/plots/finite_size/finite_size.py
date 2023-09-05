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

    q = 100
    numbs = [6, 7, 8, 9]
    gaps_left, gaps_right = [], []

    # extract many-body gap
    mb_file = f"q{q}/n6/fermions_hofstadter_X_20_Y_5_q_1_n_6_x_2_y_9_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2])/2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3]-mb_E[2])

    # extract many-body gap
    mb_file = f"q{q}/n7/fermions_hofstadter_X_50_Y_2_q_1_n_7_x_1_y_21_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract many-body gap
    mb_file = f"q{q}/n8/fermions_hofstadter_X_25_Y_4_q_1_n_8_x_2_y_12_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract many-body gap
    mb_file = f"q{q}/n9/fermions_hofstadter_X_50_Y_2_q_1_n_9_x_1_y_27_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract many-body gap
    # mb_file = f"q{q}/n10/fermions_hofstadter_X_30_Y_5_q_1_n_10_x_2_y_15_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    # mb_E = []
    # with open(mb_file, 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter=' ')
    #     for i, row in enumerate(plots):
    #         if can_convert_to_float(row[0]):
    #             mb_E.append(float(row[2]) / 2)
    # mb_E = sorted(mb_E)
    # gaps_left.append(mb_E[3] - mb_E[2])

    # extract many-body gap
    mb_file = f"q{q}/n6/fermions_hofstadter_X_20_Y_5_q_1_n_6_x_2_y_9_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract many-body gap
    mb_file = f"q{q}/n7/fermions_hofstadter_X_50_Y_2_q_1_n_7_x_1_y_21_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract many-body gap
    mb_file = f"q{q}/n8/fermions_hofstadter_X_25_Y_4_q_1_n_8_x_2_y_12_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract many-body gap
    mb_file = f"q{q}/n9/fermions_hofstadter_X_50_Y_2_q_1_n_9_x_1_y_27_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract many-body gap
    # mb_file = f"q{q}/n10/fermions_hofstadter_X_30_Y_5_q_1_n_10_x_2_y_15_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    # mb_E = []
    # with open(mb_file, 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter=' ')
    #     for i, row in enumerate(plots):
    #         if can_convert_to_float(row[0]):
    #             mb_E.append(float(row[2]) / 2)
    # mb_E = sorted(mb_E)
    # gaps_right.append(mb_E[3] - mb_E[2])

    ##########
    # figure #
    ##########

    fig = plt.figure(figsize=(4, 2))
    gs = gridspec.GridSpec(1, 2, wspace=0.5)

    #####################
    # many-body spectra #
    #####################

    ax0 = plt.subplot(gs[0])
    ax0.plot([1/i for i in numbs], np.multiply(q**2, gaps_left), '.-')
    ax0.set_ylabel('$q^2 \\Delta_\\mathrm{m.b.}$')
    ax0.set_xlabel('$1/N$')
    ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.set_xlim(0)
    ax0.set_ylim([0, 1.5])

    ax1 = plt.subplot(gs[1])
    ax1.plot([1 / i for i in numbs], np.multiply(q**2, gaps_right), '.-')
    ax1.set_ylabel('$q^2 \\Delta_\\mathrm{m.b.}$')
    ax1.set_xlabel('$1/N$')
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.set_xlim(0)
    ax1.set_ylim(0)

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/finite_size_q{q}.png",
                bbox_inches='tight', dpi=300)
    plt.show()
