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

    stats = "bosons"  # "fermions" or "bosons"
    q = 49

    if stats == "fermions":
        s = 3  # g.s. degeneracy
    else:  # stats == "bosons"
        s = 2  # g.s. degeneracy

    for alpha in [0]:
        os.chdir(f"/home/bart/DiagHam_latest/run/SFCI_data/{stats}_alpha_{alpha}/q_{q:g}")
        out_file = open(f"mb_ener_q_{q:g}.txt", "w")

        for t3hop in np.linspace(-0.25, 0, 26):

            for file in os.listdir("."):
                if f"t3_{t3hop:.2f}_" in file and file.endswith("_0.dat"):
                    given_file = file

            PlotHofstadterSpectrum = "/home/bart/DiagHam_latest/scripts_bart/PlotHofstadterSpectrum.pl -s "
            os.system(PlotHofstadterSpectrum+given_file)

            for file in os.listdir("."):
                if f"t3_{t3hop:.2f}_" in file and file.endswith("_ext.dat"):
                    given_ext_file = file

            with open(given_ext_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                E = []
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        E.append(float(row[2]))
                E = sorted(E)
                assert np.abs(E[s]-E[s-1]) > np.abs(E[s-1]-E[s-2])  # ensure m.b. gap > g.s. degeneracy
                gap = np.abs(E[s]-E[s-1])

            out_file.write(f"{t3hop:.2f}\t{gap}\n")

        out_file.close()
