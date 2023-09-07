import numpy as np

# (lx, ly) = MUC dimensions
# (Lx, Ly) = system dimensions in units of MUC

# conditions
# nu = N / (Lx*Ly) = 1/3
# nphi = 1 / (lx * ly)


def find_optimal_config(q_target, N):
    for square_err in np.linspace(0, 0.2, 21):
        for Lx in range(1, 100):
            for Ly in range(1, 100):
                if Lx * Ly / N == 3 and Lx <= Ly:  # check filling
                    for lx in range(1, 100):
                        for ly in range(1, 100):
                            if 0 <= lx * ly - q_target <= q_err and lx >= ly:  # check nphi
                                squareness = np.abs(1 - (Lx * lx) / (Ly * ly))
                                if squareness <= square_err:  # check square total system
                                    # print(f"{lx}\t{ly}\t{Lx}\t{Ly}\t{Lx * lx}\t{Ly * ly}\t{lx * ly}\t{squareness}")
                                    return square_err
    return None


if __name__ == "__main__":

    q_err = 0
    # square_err = 0.2  # within 20%
    # N = 12

    for q_target in range(96, 151):
        value_arr = []
        for N in range(6, 11):
            value = find_optimal_config(q_target, N)
            value_arr.append(value)
        if None not in value_arr:
            print(f"q={q_target} has optimal {np.max(value_arr):g}")
        # else:
        #     print(f"q={q_target} has no optimal parameter with square_err<=0.5")

# 1) continuum limit for each N = 6, 7, 8, 9, 10
# 2) then 1/N scaling for each cont(q^2 * Delta)
# 3) potentially repeat for ent gap
