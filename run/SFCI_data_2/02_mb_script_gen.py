import numpy as np


def pot(distance):
    return np.exp(-np.abs(distance)**4)


if __name__ == "__main__":

    stats = "fermions"  # "fermions" or "bosons"
    alpha = 0
    p, X, Y, x, y = 8, 12, 8, 4, 6
    q = X * Y
    ts = np.linspace(-0.25, 0.25, 11)

    FCIHofstadterModel = "~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel"

    file = open(f"03_mb_ener_script.sh", "w")
    file.write("#!/bin/bash\n\n")

    if stats == "fermions":
        file.write(f"mkdir -p fermions_alpha_{alpha:g}\n")
        file.write(f"cd fermions_alpha_{alpha:g}\n")
        file.write(f"mkdir -p q_{q:g}\n")
        file.write(f"cd q_{q:g}\n")
        for t6hop in ts:
            for t9hop in ts:
                t3hop = (-9*t6hop - 16*t9hop - 1)/4  # quartic plane
                file.write(f"{FCIHofstadterModel} -p {p} -X {X} -Y {Y} -x {x} -y {y} "
                           f"--t3hop {t3hop:.2f} --t6hop {t6hop:.2f} --t9hop {t9hop:.2f} -m 32000 "
                           f"-S --processors 6 -n 5 --lanczos-precision 1e-10 "
                           f"--u-potential {(1-alpha)+alpha*pot(1):g} --v-potential {alpha*pot(np.sqrt(2)):g} "
                           f"--v2-potential {alpha*pot(2):g} --v3-potential {alpha*pot(np.sqrt(5)):g} > /dev/null;\n")
        file.write(f"cd ../..\n")
    else:  # stats == "bosons"
        file.write(f"mkdir -p bosons_alpha_{alpha:g}\n")
        file.write(f"cd bosons_alpha_{alpha:g}\n")
        file.write(f"mkdir -p q_{q:g}\n")
        file.write(f"cd q_{q:g}\n")
        for t6hop in ts:
            for t9hop in ts:
                t3hop = (-9*t6hop - 16*t9hop - 1)/4  # quartic plane
                file.write(f"{FCIHofstadterModel} --boson -p {p} -X {X} -Y {Y} -x {x} -y {y} "
                           f"--t3hop {t3hop:.2f} --t6hop {t6hop:g} --t9hop {t9hop:g} -m 32000 "
                           f"-S --processors 6 -n 5 --lanczos-precision 1e-10 "
                           f"--u-potential {1:g} --v-potential {alpha*pot(1):g} "
                           f"--v2-potential {alpha*pot(np.sqrt(2)):g} --v3-potential {alpha*pot(2):g} "
                           f"--v4-potential {alpha*pot(np.sqrt(5)):g} > /dev/null;\n")
        file.write(f"cd ../..\n")

    file.close()
