import numpy as np
import math

if __name__ == "__main__":

    Nx = 4
    Ny = 6
    N = 4

    number = Nx * Ny * math.factorial(Nx * Ny - 2*N - 1) / (math.factorial(N) * math.factorial(Nx*Ny - 3*N))

    print(number)
