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
    mb_file = f"q{q}/n6/ent3/fermions_hofstadter_X_12_Y_8_q_1_n_6_x_3_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2])/2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3]-mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n6/ent3/fermions_hofstadter_X_12_Y_8_q_1_n_6_x_3_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[330] - ent_E[329])  # 330

    # extract many-body gap
    mb_file = f"q{q}/n7/nonauto/fermions_hofstadter_X_16_Y_6_q_1_n_7_x_3_y_7_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n7/ent3/fermions_hofstadter_X_16_Y_6_q_1_n_7_x_3_y_7_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[630] - ent_E[629])  # 630

    # extract many-body gap
    mb_file = f"q{q}/n8/nonauto/fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n8/ent2/fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[2730] - ent_E[2729])  # 2730

    # extract many-body gap
    mb_file = f"q{q}/n9/ent3/fermions_hofstadter_X_16_Y_6_q_1_n_9_x_3_y_9_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n9/ent3/fermions_hofstadter_X_16_Y_6_q_1_n_9_x_3_y_9_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[5508] - ent_E[5507])  # 5508

    # extract many-body gap
    mb_file = f"q{q}/n10/nonauto/fermions_hofstadter_X_16_Y_6_q_1_n_10_x_3_y_10_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
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

    # extract ent gap
    ent_file = f"q{q}/n11/ent/fermions_hofstadter_X_16_Y_6_q_1_n_11_x_3_y_11_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_5.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[48279] - ent_E[48278])  # 48279

    ####################################################################################################################

    # extract many-body gap
    mb_file = f"q{q}/n6/ent3/fermions_hofstadter_X_12_Y_8_q_1_n_6_x_3_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n6/ent3/fermions_hofstadter_X_12_Y_8_q_1_n_6_x_3_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[330] - ent_E[329])  # 330

    # extract many-body gap
    mb_file = f"q{q}/n7/nonauto/fermions_hofstadter_X_16_Y_6_q_1_n_7_x_3_y_7_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n7/ent3/fermions_hofstadter_X_16_Y_6_q_1_n_7_x_3_y_7_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[630] - ent_E[629])  # 630

    # extract many-body gap
    mb_file = f"q{q}/n8/nonauto/fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n8/ent2/fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[2730] - ent_E[2729])  # 2730

    # extract many-body gap
    mb_file = f"q{q}/n9/ent3/fermions_hofstadter_X_16_Y_6_q_1_n_9_x_3_y_9_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]) / 2)
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[3] - mb_E[2])

    # extract ent gap
    ent_file = f"q{q}/n9/ent3/fermions_hofstadter_X_16_Y_6_q_1_n_9_x_3_y_9_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[5508] - ent_E[5507])  # 5508

    # extract many-body gap
    mb_file = f"q{q}/n10/nonauto/fermions_hofstadter_X_16_Y_6_q_1_n_10_x_3_y_10_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
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

    # extract ent gap
    ent_file = f"q{q}/n11/ent/fermions_hofstadter_X_16_Y_6_q_1_n_11_x_3_y_11_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_5.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[48279] - ent_E[48278])  # 48279

    ####################################################################################################################

    f = open(f"data_fermions_q{q}.txt", "w")
    for idx, i in enumerate(numbs):
        f.write(f"{i}\t{q**2 * gaps_left[idx]}\t{q**2 * gaps_right[idx]}\t{ent_left[idx]}\t{ent_right[idx]}\n")
    f.close()
