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

    q = 96
    numbs = [6, 7, 8, 9, 10, 11]
    gaps_left, gaps_right = [], []
    ent_left, ent_right = [], []

    # extract many-body gap
    mb_file = f"q{q}/n6/fermions_hofstadter_X_12_Y_8_q_1_n_6_x_3_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2])/2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3]-mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n6/ent/fermions_hofstadter_X_12_Y_8_q_1_n_6_x_3_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[762] - ent_E[761])  # 762

    # extract many-body gap
    mb_file = f"q{q}/n7/fermions_hofstadter_X_16_Y_6_q_1_n_7_x_3_y_7_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n7/ent/fermions_hofstadter_X_16_Y_6_q_1_n_7_x_3_y_7_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[637] - ent_E[636])  # 637

    # extract many-body gap
    mb_file = f"q{q}/n8/fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n8/ent/fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[2730] - ent_E[2729])  # 2730

    # extract many-body gap
    mb_file = f"q{q}/n9/fermions_hofstadter_X_16_Y_6_q_1_n_9_x_3_y_9_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n9/ent/fermions_hofstadter_X_16_Y_6_q_1_n_9_x_3_y_9_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[5508] - ent_E[5507])  # 5508

    # extract many-body gap
    mb_file = f"q{q}/n10/fermions_hofstadter_X_16_Y_6_q_1_n_10_x_3_y_10_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n10/ent/fermions_hofstadter_X_16_Y_6_q_1_n_10_x_3_y_10_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_5.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[23256] - ent_E[23255])  # 23256

    # extract many-body gap
    mb_file = f"q{q}/n11/fermions_hofstadter_X_16_Y_6_q_1_n_11_x_3_y_11_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    ####################################################################################################################

    # extract many-body gap
    mb_file = f"q{q}/n6/fermions_hofstadter_X_12_Y_8_q_1_n_6_x_3_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n6/ent/fermions_hofstadter_X_12_Y_8_q_1_n_6_x_3_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[762] - ent_E[761])  # 762

    # extract many-body gap
    mb_file = f"q{q}/n7/fermions_hofstadter_X_16_Y_6_q_1_n_7_x_3_y_7_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n7/ent/fermions_hofstadter_X_16_Y_6_q_1_n_7_x_3_y_7_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[637] - ent_E[636])  # 637

    # extract many-body gap
    mb_file = f"q{q}/n8/fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n8/ent/fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[2730] - ent_E[2729])  # 2730

    # extract many-body gap
    mb_file = f"q{q}/n9/fermions_hofstadter_X_16_Y_6_q_1_n_9_x_3_y_9_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n9/ent/fermions_hofstadter_X_16_Y_6_q_1_n_9_x_3_y_9_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[5508] - ent_E[5507])  # 5508

    # extract many-body gap
    mb_file = f"q{q}/n10/fermions_hofstadter_X_16_Y_6_q_1_n_10_x_3_y_10_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n10/ent/fermions_hofstadter_X_16_Y_6_q_1_n_10_x_3_y_10_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_5.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[23256] - ent_E[23255])  # 23256

    # extract many-body gap
    mb_file = f"q{q}/n11/fermions_hofstadter_X_16_Y_6_q_1_n_11_x_3_y_11_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    ##########
    # figure #
    ##########

    fig = plt.figure(figsize=(4, 2))
    gs = gridspec.GridSpec(2, 2, wspace=0.5)

    #####################
    # many-body spectra #
    #####################

    ent_left.append(0)
    ent_right.append(0)
    ratio = [0.75, 1.1428571428571428, 1.0, 0.8888888888888888, 0.8, 0.7272727272727273]
    f = open("finite_size.txt", "w")
    for idx, i in enumerate(numbs):
        f.write(f"{i}\t{ratio[idx]}\t{q**2 * gaps_left[idx]}\t{q**2 * gaps_right[idx]}\t"
                f"{ent_left[idx]}\t{ent_right[idx]}\n")
    f.close()

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

    ax2 = plt.subplot(gs[2])
    ax2.plot([1 / i for i in numbs], ent_left, '.-')

    ax3 = plt.subplot(gs[3])
    ax3.plot([1 / i for i in numbs], ent_right, '.-')

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/finite_size_q{q}.png",
                bbox_inches='tight', dpi=300)
    plt.show()
