import os
import csv
import numpy as np


def can_convert_to_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


if __name__ == "__main__":

    program = "turn"  # "tune" or "turn"
    stats = "fermions"  # "fermions" or "bosons"

    if stats == "fermions":
        s = 3  # g.s. degeneracy
        N_list = [6]
        p_list = [9]  # 9, 19, 33, 51, 73, 99, 129, 163
    else:
        s = 2  # g.s. degeneracy
        N_list = [6]
        p_list = [13, 28, 49, 76, 109, 148, 193]  # 13, 28, 49, 76, 109, 148, 193

    ind_var = "alpha" if program == "tune" else "theta"
    out_file = open(f"{program}/{stats}/{stats}_{ind_var}_gap.dat", "w")

    for N in N_list:
        for p in p_list:

            os.chdir(f"/home/bart/DiagHam_latest/run/hof_int/{program}/{stats}/p_{p}_N_{N}")

            file_list = []
            for file in os.listdir("."):
                if file.endswith("_0.dat"):
                    file_list.append(file)
            file_list = sorted(file_list)

            PlotHofstadterSpectrum = "/home/bart/DiagHam_latest/scripts_bart/PlotHofstadterSpectrum.pl -s "
            for i in file_list:
                os.system(PlotHofstadterSpectrum+i)

            ext_file_list = []
            for file in os.listdir("."):
                if file.endswith("_ext.dat"):
                    ext_file_list.append(file)
            print(ext_file_list)
            ext_file_list = sorted(ext_file_list)

            gaps = []
            for i in ext_file_list:
                print(i)
                E = []
                with open(i, 'r') as csvfile:
                    plots = csv.reader(csvfile, delimiter=' ')
                    for j, row in enumerate(plots):
                        if can_convert_to_float(row[0]):  # if value is a number
                            E.append(float(row[2]))
                    E = sorted(E)
                    gaps.append(np.abs(E[s]-E[s-1]))

            for i, gap in enumerate(gaps):
                out_file.write(f"{N}\t{p}\t{i/(len(gaps)-1)}\t{gap}\n")
