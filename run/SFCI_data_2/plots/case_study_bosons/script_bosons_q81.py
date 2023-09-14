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

    q = 81
    numbs = [6, 7, 8, 9, 10]
    gaps_left, gaps_right = [], []
    ent_left, ent_right = [], []

    # extract many-body gap
    mb_file = f"q{q}/n6/nonauto/bosons_hofstadter_X_9_Y_9_q_1_n_6_x_3_y_4_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[2]-mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n6/ent2/bosons_hofstadter_X_9_Y_9_q_1_n_6_x_3_y_4_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[112] - ent_E[111])  # 112

    # extract many-body gap
    mb_file = f"q{q}/n7/nonauto/bosons_hofstadter_X_27_Y_3_q_1_n_7_x_1_y_14_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n7/ent3/bosons_hofstadter_X_27_Y_3_q_1_n_7_x_1_y_14_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[210] - ent_E[209])  # 210

    # extract many-body gap
    mb_file = f"q{q}/n8/nonauto/bosons_hofstadter_X_9_Y_9_q_1_n_8_x_4_y_4_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n8/ent2/bosons_hofstadter_X_9_Y_9_q_1_n_8_x_4_y_4_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[660] - ent_E[659])  # 660

    # extract many-body gap
    mb_file = f"q{q}/n9/nonauto/bosons_hofstadter_X_9_Y_9_q_1_n_9_x_3_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n9/ent2/bosons_hofstadter_X_9_Y_9_q_1_n_9_x_3_y_6_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[1287] - ent_E[1286])  # 1287

    # extract many-body gap
    mb_file = f"q{q}/n10/nonauto/bosons_hofstadter_X_9_Y_9_q_1_n_10_x_4_y_5_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_left.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n10/ent2/bosons_hofstadter_X_9_Y_9_q_1_n_10_x_4_y_5_t3_1.31_t6_-0.25_t9_-0.25_u_1_gx_0_gy_0.na_5.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_left.append(ent_E[4004] - ent_E[4003])  # 4004

    ####################################################################################################################

    # extract many-body gap
    mb_file = f"q{q}/n6/nonauto/bosons_hofstadter_X_9_Y_9_q_1_n_6_x_3_y_4_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n6/ent2/bosons_hofstadter_X_9_Y_9_q_1_n_6_x_3_y_4_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[112] - ent_E[111])  # 112

    # extract many-body gap
    mb_file = f"q{q}/n7/nonauto/bosons_hofstadter_X_27_Y_3_q_1_n_7_x_1_y_14_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n7/ent3/bosons_hofstadter_X_27_Y_3_q_1_n_7_x_1_y_14_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_3.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[210] - ent_E[209])  # 210

    # extract many-body gap
    mb_file = f"q{q}/n8/nonauto/bosons_hofstadter_X_9_Y_9_q_1_n_8_x_4_y_4_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n8/ent2/bosons_hofstadter_X_9_Y_9_q_1_n_8_x_4_y_4_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[660] - ent_E[659])  # 660

    # extract many-body gap
    mb_file = f"q{q}/n9/nonauto/bosons_hofstadter_X_9_Y_9_q_1_n_9_x_3_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n9/ent2/bosons_hofstadter_X_9_Y_9_q_1_n_9_x_3_y_6_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_4.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[1287] - ent_E[1286])  # 1287

    # extract many-body gap
    mb_file = f"q{q}/n10/nonauto/bosons_hofstadter_X_9_Y_9_q_1_n_10_x_4_y_5_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0_ext.dat"
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                mb_E.append(float(row[2]))
    mb_E = sorted(mb_E)
    gaps_right.append(mb_E[2] - mb_E[1])

    # extract ent gap
    ent_file = f"q{q}/n10/ent2/bosons_hofstadter_X_9_Y_9_q_1_n_10_x_4_y_5_t3_-1.81_t6_0.25_t9_0.25_u_1_gx_0_gy_0.na_5.parentspec"
    ent_E = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                ent_E.append(float(row[5]) * 2)
    ent_E = sorted(ent_E)
    ent_right.append(ent_E[4004] - ent_E[4003])  # 4004

    ####################################################################################################################

    f = open(f"data_bosons_q{q}.txt", "w")
    for idx, i in enumerate(numbs):
        f.write(f"{i}\t{q * gaps_left[idx]}\t{q * gaps_right[idx]}\t{ent_left[idx]}\t{ent_right[idx]}\n")
    f.close()
