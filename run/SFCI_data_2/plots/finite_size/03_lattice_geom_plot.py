import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

# (lx, ly) = MUC dimensions
# (Lx, Ly) = system dimensions in units of MUC

# conditions
# nu = N / (Lx*Ly) = 1/3
# nphi = 1 / (lx * ly)


def can_convert_to_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def find_optimal_parameters(q_target, q_err, square_err, N):

    data = []

    for Lx in range(1, 100):
        for Ly in range(1, 100):
            if Lx*Ly/N == 3 and Lx <= Ly:  # check filling
                for lx in range(1, 100):
                    for ly in range(1, 100):
                        if 0 <= lx*ly-q_target <= q_err and lx >= ly:  # check nphi
                            squareness = np.abs(1 - (Lx*lx)/(Ly*ly))
                            if squareness <= square_err:  # check square total system
                                print(f"{lx}\t{ly}\t{Lx}\t{Ly}\t{Lx*lx}\t{Ly*ly}\t{lx*ly}\t{squareness}")
                                data.append([lx, ly, Lx, Ly, squareness])
    data = np.array(data)
    data = sorted(data, key=lambda x:(x[4], -x[2]))
    print(data)

    return int(data[0][0]), int(data[0][1]), int(data[0][2]), int(data[0][3])


if __name__ == "__main__":

    q_list = [96]
    N_list = [6, 7, 8, 9, 10, 11]

    for q in q_list:
        gaps_left, gaps_right = [], []
        numbs = []
        for N in N_list:
            print(f"(q, N) = ({q}, {N})")

            # change directory
            path = f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q{q:g}/n{N:g}"
            os.chdir(path)

            X, Y, x, y = find_optimal_parameters(q, 0, 0.5, N)

            # continue to next q value if no optimal geometry found
            if X == Y == x == y == 99:
                break

            print(X, Y, x, y)

            # plot the energy spectra
            PlotHofstadterSpectrum = "/home/bart/DiagHam_latest/scripts_bart/PlotHofstadterSpectrum.pl -s "
            for file in os.listdir("."):
                if f"t6_-0.25_t9_-0.25_" in file and file.endswith("_0.dat"):
                    given_file = file
            os.system(PlotHofstadterSpectrum+given_file)
            for file in os.listdir("."):
                if f"t6_0.25_t9_0.25_" in file and file.endswith("_0.dat"):
                    given_file = file
            os.system(PlotHofstadterSpectrum+given_file)

            # extract the gaps
            for file in os.listdir("."):
                if "t6_-0.25_t9_-0.25_" in file and file.endswith("_ext.dat"):
                    given_ext_file = file
            with open(given_ext_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                E = []
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        E.append(float(row[2]))
                E = sorted(E)
                gaps_left.append(np.abs(E[3] - E[2]))
            for file in os.listdir("."):
                if "t6_0.25_t9_0.25_" in file and file.endswith("_ext.dat"):
                    given_ext_file = file
            with open(given_ext_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                E = []
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        E.append(float(row[2]))
                E = sorted(E)
                gaps_right.append(np.abs(E[3] - E[2]))
            numbs.append(N)

        # plot the figure
        fig = plt.figure(figsize=(4, 2))
        gs = gridspec.GridSpec(1, 2, wspace=0.5)
        #
        ax0 = plt.subplot(gs[0])
        ax0.plot([1 / i for i in numbs], np.multiply(q**2 / 2, gaps_left), '.-')
        ax0.set_ylabel('$q^2 \\Delta_\\mathrm{m.b.}$')
        ax0.set_xlabel('$1/N$')
        ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax0.set_xlim(0)
        ax0.set_ylim(0)
        #
        ax1 = plt.subplot(gs[1])
        ax1.plot([1 / i for i in numbs], np.multiply(q**2 / 2, gaps_right), '.-')
        ax1.set_ylabel('$q^2 \\Delta_\\mathrm{m.b.}$')
        ax1.set_xlabel('$1/N$')
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax1.set_xlim(0)
        ax1.set_ylim(0)

        plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/finite_size_q{q}.png",
                    bbox_inches='tight', dpi=300)
        plt.show()
