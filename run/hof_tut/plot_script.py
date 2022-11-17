import os
import numpy as np
import csv
import matplotlib.pyplot as plt

if __name__ == "__main__":

    gap_executable = "/home/bart/DiagHam_Nimit/scripts_bart/PlotHofstadterSpectrum.pl"

    os.chdir("fermions")  # enter the parent directory

    q_list = np.arange(3, 7)  # integer list from 3 to 6
    data = np.zeros((len(q_list), 2))  # initialize a data array with q rows and 2 columns (q and Delta)

    for i, q in enumerate(q_list):  # loop over integer q from 3 to 6
        x, y = 3, 6  # number of MUCs
        os.chdir(f"x_{x}_y_{y}_X_{q}_Y_3")  # enter the directory
        gap_flags = f" -s fermions_hofstadter_X_{q}_Y_3_q_1_n_6_x_{x}_y_{y}_u_1_gx_0_gy_0.dat"
        gap_command = gap_executable+gap_flags  # define the command
        os.system(gap_command)  # execute the command

        with open(f"fermions_hofstadter_X_{q}_Y_3_q_1_n_6_x_{x}_y_{y}_u_1_gx_0_gy_0_ext.dat", 'r') as csvfile:
            file = csv.reader(csvfile, delimiter=' ')  # rows of file are defined by new lines, columns by spaces
            energies = []  # initialize an empty energies list
            for j, row in enumerate(file):
                if j > 1:  # skip first two lines of file
                    energies.append(float(row[2]))  # append the energies to the list
            energies = sorted(energies)  # sort the energies in ascending order
            data[i, 0] = q  # save the q values in the first column of data array
            data[i, 1] = np.abs(energies[1]-energies[0])  # save the Delta values in the second column of data array

        os.chdir("..")  # leave the directory

    print(data)

    plt.figure()
    ax = plt.subplot(111)
    ax.plot(data[:, 0], data[:, 1], '.-')
    ax.set_title("fermions at 1/3 filling")
    ax.set_xlabel("q")
    ax.set_ylabel("Delta")
    plt.savefig("fermions.png", bbox_inches='tight', dpi=300)
    plt.show()
