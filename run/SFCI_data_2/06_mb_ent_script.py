import os
import csv
import numpy as np
import shutil


def can_convert_to_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def pot(distance):
    return np.exp(-np.abs(distance)**4)


if __name__ == "__main__":

    stats = "fermions"  # "fermions" or "bosons"
    alpha = 0
    p, X, Y, x, y = 8, 15, 10, 4, 6
    q = X*Y
    ts = np.linspace(-0.25, 0.25, 11)

    if stats == "fermions":
        s = 3  # g.s. degeneracy
    else:  # stats == "bosons"
        s = 2  # g.s. degeneracy

    ener_path = f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats}_alpha_{alpha}/q_{q:g}"
    ent_path = os.path.join(ener_path, "ent")
    os.makedirs(ent_path, exist_ok=True)
    os.chdir(ent_path)

    for t6hop in ts:
        for t9hop in ts:
            t3hop = (-9 * t6hop - 16 * t9hop - 1) / 4  # quartic plane

            # extract the g.s. momentum sectors
            for file in os.listdir(".."):
                if f"t6_{t6hop:.2f}_t9_{t9hop:.2f}_" in file and file.endswith("_0_ext.dat"):
                    given_ext_file = file

            with open("../"+given_ext_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                kx, ky, E = [], [], []
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        kx.append(int(row[0]))
                        ky.append(int(row[1]))
                        E.append(float(row[2]))
                # sortedE = sorted(E)  # ensure m.b. gap > g.s. degeneracy
                # assert np.abs(sortedE[s] - sortedE[s - 1]) > np.abs(sortedE[s - 1] - sortedE[s - 2])

            kx_min, ky_min, E_min = [], [], []
            for i in range(s):
                kx_min.append(kx[np.argsort(E)[i]])
                ky_min.append(ky[np.argsort(E)[i]])
                E_min.append(E[np.argsort(E)[i]])
                print(f"(t6, t9) = ({t6hop:.2f}, {t9hop:.2f}) ==> ({kx_min[i]}, {ky_min[i]}): {E_min[i]}")

            # create the taskfile
            given_file = given_ext_file.replace("_ext", "")
            string = f"{given_file}\t{s}\t"
            for i in range(s):
                string += f"{kx_min[i]}\t{ky_min[i]}\t"

            taskfile_name = f"taskfile_t6_{t6hop:.2f}_t9_{t9hop:.2f}"
            with open(taskfile_name, 'w') as file:
                file.write(f"{string}")

            # get the g.s. eigenvectors
            FCIHofstadterModel = "/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel"
            for i in range(s):

                if i > 0:  # don't compute the same eigenvectors again
                    if kx_min[i] == kx_min[i-1] and ky_min[i] == ky_min[i-1]:
                        continue

                if stats == "fermions":
                    command = f"{FCIHofstadterModel} -p {p} -X {X} -Y {Y} -x {x} -y {y} " \
                              f"--t3hop {t3hop:.2f} --t6hop {t6hop:.2f} --t9hop {t9hop:.2f} -m 32000 " \
                              f"-S --processors 6 -n 5 --lanczos-precision 1e-10 " \
                              f"--u-potential {(1-alpha)+alpha*pot(1):g} --v-potential {alpha*pot(np.sqrt(2)):g} " \
                              f"--v2-potential {alpha*pot(2):g} --v3-potential {alpha*pot(np.sqrt(5)):g} " \
                              f"--only-kx {kx_min[i]} --only-ky {ky_min[i]} --eigenstate > /dev/null;"
                else:  # stats == "bosons"
                    command = f"{FCIHofstadterModel} --boson -p {p} -X {X} -Y {Y} -x {x} -y {y} " \
                              f"--t3hop {t3hop:.2f} --t6hop {t6hop:.2f} --t9hop {t9hop:.2f} -m 32000 " \
                              f"-S --processors 6 -n 5 --lanczos-precision 1e-10 " \
                              f"--u-potential {1:g} --v-potential {alpha*pot(1):g} " \
                              f"--v2-potential {alpha*pot(np.sqrt(2)):g} --v3-potential {alpha*pot(2):g} " \
                              f"--v4-potential {alpha*pot(np.sqrt(5)):g} " \
                              f"--only-kx {kx_min[i]} --only-ky {ky_min[i]} --eigenstate > /dev/null;"
                os.system(command)

            # create the entanglement spectrum
            GetHofstadterEntanglementSpectra = \
                "/home/bart/DiagHam_latest/scripts_bart/entanglement/GetHofstadterEntanglementSpectra.pl "
            os.system(f"eval $({GetHofstadterEntanglementSpectra+taskfile_name})")

            # plot the entanglement spectrum
            PlotHofstadterFTIEntanglementSpectrum = \
                "/home/bart/DiagHam_latest/scripts_bart/entanglement/PlotHofstadterFTIEntanglementSpectrum.pl "
            os.system(PlotHofstadterFTIEntanglementSpectrum+given_file.replace(".dat", ".full.parent"))

    # clean up
    # ent_path = os.path.join(ener_path, "ent")
    # os.makedirs(ent_path, exist_ok=True)
    # for file in os.listdir("."):
    #     if "task" in file or "par" in file or file.endswith(".vec") or file.endswith(".gs"):
    #         shutil.move(os.path.join(ener_path, file), os.path.join(ent_path, file))
