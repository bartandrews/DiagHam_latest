import numpy as np


def pot(distance):
    return np.exp(-np.abs(distance)**4)


if __name__ == "__main__":

    stats = "bosons"  # "fermions" or "bosons"
    q = 81

    FCIHofstadterModel = "~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel"

    file = open(f"02_mb_script.sh", "w")
    file.write("#!/bin/bash\n\n")

    if stats == "fermions":
        for alpha in [1]:
            file.write(f"mkdir -p fermions_alpha_{alpha:g}\n")
            file.write(f"cd fermions_alpha_{alpha:g}\n")
            file.write(f"mkdir -p q_{q:g}\n")
            file.write(f"cd q_{q:g}\n")
            for t3hop in np.linspace(-0.25, 0, 26):
                file.write(f"{FCIHofstadterModel} -p 8 -X 12 -Y 8 -x 4 -y 6 --t3hop {t3hop:.2f} -m 32000 "
                           f"-S --processors 6 -n 5 --lanczos-precision 1e-10 "
                           f"--u-potential {(1-alpha)+alpha*pot(1):g} --v-potential {alpha*pot(np.sqrt(2)):g} "
                           f"--v2-potential {alpha*pot(2):g} --v3-potential {alpha*pot(np.sqrt(5)):g} > /dev/null;\n")
            file.write(f"cd ../..\n")
    else:  # stats == "bosons"
        for alpha in [1]:
            file.write(f"mkdir -p bosons_alpha_{alpha:g}\n")
            file.write(f"cd bosons_alpha_{alpha:g}\n")
            file.write(f"mkdir -p q_{q:g}\n")
            file.write(f"cd q_{q:g}\n")
            for t3hop in np.linspace(-0.25, 0, 26):
                file.write(f"{FCIHofstadterModel} --boson -p 8 -X 9 -Y 9 -x 4 -y 4 --t3hop {t3hop:.2f} -m 32000 "
                           f"-S --processors 6 -n 5 --lanczos-precision 1e-10 "
                           f"--u-potential {1:g} --v-potential {alpha*pot(1):g} "
                           f"--v2-potential {alpha*pot(np.sqrt(2)):g} --v3-potential {alpha*pot(2):g} "
                           f"--v4-potential {alpha*pot(np.sqrt(5)):g} > /dev/null;\n")
            file.write(f"cd ../..\n")

    file.close()
