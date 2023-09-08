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
    # print(data)

    return int(data[0][0]), int(data[0][1]), int(data[0][2]), int(data[0][3])


if __name__ == "__main__":

    file = open(f"02_lattice_geom_script.sh", "w")
    file.write("#!/bin/bash\n\n")

    q_list = [128]
    N_list = [6, 7, 8, 9]

    for q in q_list:
        gaps_left, gaps_right = [], []
        numbs = []
        for N in N_list:
            print(f"(q, N) = ({q}, {N})")

            # make directory
            path = f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q{q:g}/n{N:g}"
            file.write(f"mkdir -p {path}\n")
            file.write(f"cd {path}\n")

            X, Y, x, y = find_optimal_parameters(q, 0, 0.5, N)

            # continue to next q value if no optimal geometry found
            if X == Y == x == y == 99:
                break

            print(X, Y, x, y)

            # write the commands
            FCIHofstadterModel = "/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel "
            file.write(f"{FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} "
                       f"--t3hop 1.31 --t6hop -0.25 --t9hop -0.25 -m 96000 "
                       f"-S --processors 6 -n 5 --lanczos-precision 1e-10 "
                       f"--u-potential 1 --v-potential 0 "
                       f"--v2-potential 0 --v3-potential 0 --auto-addprojector;\n")
            file.write(f"{FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} "
                       f"--t3hop -1.81 --t6hop 0.25 --t9hop 0.25 -m 96000 "
                       f"-S --processors 6 -n 5 --lanczos-precision 1e-10 "
                       f"--u-potential 1 --v-potential 0 "
                       f"--v2-potential 0 --v3-potential 0 --auto-addprojector;\n")
