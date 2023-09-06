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
    N_list = [6, 7, 8, 9, 10]

    for q in q_list:
        gaps_left, gaps_right = [], []
        ents_left, ents_right = [], []
        numbs = []
        for N in N_list:
            print(f"(q, N) = ({q}, {N})")

            # make directory
            path = f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output_2/q_{q:g}/n{N:g}"
            os.makedirs(path, exist_ok=True)
            os.chdir(path)

            X, Y, x, y = find_optimal_parameters(q, 0, 0.5, N)

            # continue to next q value if no optimal geometry found
            if X == Y == x == y == 99:
                break

            print(X, Y, x, y)

            # generate the energy spectra
            FCIHofstadterModel = "/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel"
            os.system(f"{FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} "
                      f"--t3hop 1.31 --t6hop -0.25 --t9hop -0.25 -m 32000 "
                      f"-S --processors 6 -n 5 --lanczos-precision 1e-10 "
                      f"--u-potential 1 --v-potential 0 "
                      f"--v2-potential 0 --v3-potential 0 --auto-addprojector > /dev/null")
            os.system(f"{FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} "
                      f"--t3hop -1.81 --t6hop 0.25 --t9hop 0.25 -m 32000 "
                      f"-S --processors 6 -n 5 --lanczos-precision 1e-10 "
                      f"--u-potential 1 --v-potential 0 "
                      f"--v2-potential 0 --v3-potential 0 --auto-addprojector > /dev/null")

            # plot the energy spectra
            PlotHofstadterSpectrum = "/home/bart/DiagHam_latest/scripts_bart/PlotHofstadterSpectrum.pl -s "
            for file in os.listdir("."):
                if f"t6_-0.25_t9_-0.25_" in file and file.endswith("_0.dat"):
                    given_file1 = file
            os.system(PlotHofstadterSpectrum+given_file1)
            for file in os.listdir("."):
                if f"t6_0.25_t9_0.25_" in file and file.endswith("_0.dat"):
                    given_file2 = file
            os.system(PlotHofstadterSpectrum+given_file2)

            # extract the gaps
            for file in os.listdir("."):
                if f"t6_-0.25_t9_-0.25_" in file and file.endswith("_ext.dat"):
                    given_ext_file = file
            with open(given_ext_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                E = []
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        E.append(float(row[2]))
                E = sorted(E)
                # if np.abs(E[s]-E[s-1]) > np.abs(E[s-1]-E[s-2]):  # ensure m.b. gap > g.s. degeneracy
                #     gap = np.abs(E[s]-E[s-1])
                # else:
                #     gap = np.nan
                gaps_left.append(np.abs(E[3] - E[2]))
            for file in os.listdir("."):
                if f"t6_0.25_t9_0.25_" in file and file.endswith("_ext.dat"):
                    given_ext_file = file
            with open(given_ext_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                E = []
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        E.append(float(row[2]))
                E = sorted(E)
                # if np.abs(E[s]-E[s-1]) > np.abs(E[s-1]-E[s-2]):  # ensure m.b. gap > g.s. degeneracy
                #     gap = np.abs(E[s]-E[s-1])
                # else:
                #     gap = np.nan
                gaps_right.append(np.abs(E[3] - E[2]))
            numbs.append(N)

            ent_path = os.path.join(path, "ent")
            os.makedirs(ent_path, exist_ok=True)
            os.chdir(ent_path)

            # extract the g.s. momentum sectors and create taskfile
            for file in os.listdir(".."):
                if f"t6_-0.25_t9_-0.25_" in file and file.endswith("_0_ext.dat"):
                    given_ext_file = file
            with open("../" + given_ext_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                kx, ky, E = [], [], []
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        kx.append(int(row[0]))
                        ky.append(int(row[1]))
                        E.append(float(row[2]))
                # sortedE = sorted(E)  # ensure m.b. gap > g.s. degeneracy
                # assert np.abs(sortedE[s] - sortedE[s - 1]) > np.abs(sortedE[s - 1] - sortedE[s - 2])
            kx_min1, ky_min1, E_min1 = [], [], []
            for i in range(3):
                kx_min1.append(kx[np.argsort(E)[i]])
                ky_min1.append(ky[np.argsort(E)[i]])
                E_min1.append(E[np.argsort(E)[i]])
                print(f"(t6, t9) = (-0.25, -0.25) ==> ({kx_min1[i]}, {ky_min1[i]}): {E_min1[i]}")
            #
            # given_file1 = given_ext_file.replace("_ext", "")
            string = f"{given_file1}\t3\t"
            for i in range(3):
                string += f"{kx_min1[i]}\t{ky_min1[i]}\t"
            taskfile_name1 = f"taskfile_t6_-0.25_t9_-0.25"
            with open(taskfile_name1, 'w') as file:
                file.write(f"{string}")

            for file in os.listdir(".."):
                if f"t6_0.25_t9_0.25_" in file and file.endswith("_0_ext.dat"):
                    given_ext_file2 = file
            with open("../" + given_ext_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                kx, ky, E = [], [], []
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        kx.append(int(row[0]))
                        ky.append(int(row[1]))
                        E.append(float(row[2]))
                # sortedE = sorted(E)  # ensure m.b. gap > g.s. degeneracy
                # assert np.abs(sortedE[s] - sortedE[s - 1]) > np.abs(sortedE[s - 1] - sortedE[s - 2])
            kx_min2, ky_min2, E_min2 = [], [], []
            for i in range(3):
                kx_min2.append(kx[np.argsort(E)[i]])
                ky_min2.append(ky[np.argsort(E)[i]])
                E_min2.append(E[np.argsort(E)[i]])
                print(f"(t6, t9) = (0.25, 0.25) ==> ({kx_min2[i]}, {ky_min2[i]}): {E_min2[i]}")
            #
            # given_file2 = given_ext_file.replace("_ext", "")
            string = f"{given_file2}\t3\t"
            for i in range(3):
                string += f"{kx_min2[i]}\t{ky_min2[i]}\t"
            taskfile_name2 = f"taskfile_t6_0.25_t9_0.25"
            with open(taskfile_name2, 'w') as file:
                file.write(f"{string}")

            # get the g.s. eigenvectors
            for i in range(3):
                if i > 0:  # don't compute the same eigenvectors again
                    if kx_min1[i] == kx_min1[i - 1] and ky_min1[i] == ky_min1[i - 1]:
                        continue
                command = f"{FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} " \
                          f"--t3hop 1.31 --t6hop -0.25 --t9hop -0.25 -m 32000 " \
                          f"-S --processors 6 -n 5 --lanczos-precision 1e-10 " \
                          f"--u-potential 1 --v-potential 0 " \
                          f"--v2-potential 0 --v3-potential 0 " \
                          f"--only-kx {kx_min1[i]} --only-ky {ky_min1[i]} --eigenstate > /dev/null;"
                os.system(command)
            for i in range(3):
                if i > 0:  # don't compute the same eigenvectors again
                    if kx_min2[i] == kx_min2[i - 1] and ky_min2[i] == ky_min2[i - 1]:
                        continue
                command = f"{FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} " \
                          f"--t3hop -1.81 --t6hop 0.25 --t9hop 0.25 -m 32000 " \
                          f"-S --processors 6 -n 5 --lanczos-precision 1e-10 " \
                          f"--u-potential 1 --v-potential 0 " \
                          f"--v2-potential 0 --v3-potential 0 " \
                          f"--only-kx {kx_min2[i]} --only-ky {ky_min2[i]} --eigenstate > /dev/null;"
                os.system(command)

            # create the entanglement spectrum
            GetHofstadterEntanglementSpectra = \
                "/home/bart/DiagHam_latest/scripts_bart/entanglement/GetHofstadterEntanglementSpectra.pl "
            os.system(f"eval $({GetHofstadterEntanglementSpectra + taskfile_name1})")
            os.system(f"eval $({GetHofstadterEntanglementSpectra + taskfile_name2})")

            # plot the entanglement spectrum
            PlotHofstadterFTIEntanglementSpectrum = \
                "/home/bart/DiagHam_latest/scripts_bart/entanglement/PlotHofstadterFTIEntanglementSpectrum.pl "
            os.system(PlotHofstadterFTIEntanglementSpectrum + given_file1.replace(".dat", ".full.parent"))
            os.system(PlotHofstadterFTIEntanglementSpectrum + given_file2.replace(".dat", ".full.parent"))

            # extract the ent spec file
            for file in os.listdir("."):
                if f"t6_-0.25_t9_-0.25_" in file and file.endswith(".parentspec"):
                    spec_file1 = file
            for file in os.listdir("."):
                if f"t6_0.25_t9_0.25_" in file and file.endswith(".parentspec"):
                    spec_file2 = file

            # compute the spectral gap
            ents = []
            with open(spec_file1, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        ents.append(float(row[5]))
            ents_left.append(ents[2730] - ents[2729])
            ents = []
            with open(spec_file2, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        ents.append(float(row[5]))
            ents_right.append(ents[2730] - ents[2729])

        # plot the figure
        fig = plt.figure(figsize=(4, 2))
        gs = gridspec.GridSpec(1, 2, wspace=0.5)
        #
        ax0 = plt.subplot(gs[0])
        ax0.plot([1 / i for i in numbs], np.multiply(q**2, gaps_left), '.-')
        ax0.set_ylabel('$q^2 \\Delta_\\mathrm{m.b.}$')
        ax0.set_xlabel('$1/N$')
        ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax0.set_xlim(0)
        ax0.set_ylim([0, 1.5])
        #
        ax1 = plt.subplot(gs[1])
        ax1.plot([1 / i for i in numbs], np.multiply(q**2, gaps_right), '.-')
        ax1.set_ylabel('$q^2 \\Delta_\\mathrm{m.b.}$')
        ax1.set_xlabel('$1/N$')
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax1.set_xlim(0)
        ax1.set_ylim(0)

        plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output_2/finite_size_q{q}.png",
                    bbox_inches='tight', dpi=300)
        plt.show()
