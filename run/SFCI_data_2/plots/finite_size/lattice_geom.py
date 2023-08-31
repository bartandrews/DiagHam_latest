import numpy as np

# (lx, ly) = MUC dimensions
# (Lx, Ly) = system dimensions in units of MUC

# conditions
# nu = N / (Lx*Ly) = 1/3
# nphi = 1 / (lx * ly)

if __name__ == "__main__":

    q_target = 120
    q_err = 0
    square_err = 0.2
    N = 6

    for Lx in range(1, 100):
        for Ly in range(1, 100):
            if Lx*Ly/N == 3 and Lx <= Ly:  # check filling
                for lx in range(1, 100):
                    for ly in range(1, 100):
                        if 0 <= lx*ly-q_target <= q_err and lx >= ly:  # check nphi
                            squareness = np.abs(1 - (Lx*lx)/(Ly*ly))
                            if squareness <= square_err:  # check square total system
                                print(f"{lx}\t{ly}\t{Lx}\t{Ly}\t{Lx*lx}\t{Ly*ly}\t{lx*ly}\t{squareness}")