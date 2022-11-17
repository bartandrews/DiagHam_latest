import numpy as np
import math


def pot(distance):
    return np.exp(-np.abs(distance)**4)


if __name__ == "__main__":

    stats = "fermions"  # "fermions" or "bosons"
    FCIHofstadterModel = "~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel"

    file = open(f"02_tune_script.sh", "w")
    file.write("#!/bin/bash\n\n")

    file.write("runs() {\n")
    for alpha in np.linspace(0, 1, 11):
        if stats == "fermions":
            file.write(f"echo {FCIHofstadterModel} -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 "
                       f"--lanczos-precision 1e-10 "
                       f"--u-potential {pot(1.0)} --v-potential {alpha*math.e*pot(np.sqrt(2))} "
                       f"--v2-potential {alpha*math.e*pot(2.0)} --v3-potential {alpha*math.e*pot(np.sqrt(5))}\n")
        else:
            file.write(f"echo {FCIHofstadterModel} --boson -p 6 -X 16 -Y 12 -x 3 -y 4 -m 8000 -S --processors 4 -n 20 "
                       f"--lanczos-precision 1e-10 "
                       f"--u-potential {pot(0.0)} --v-potential {alpha * pot(1.0)} "
                       f"--v2-potential {alpha * pot(np.sqrt(2))} --v3-potential {alpha * pot(2.0)} "
                       f"--v4-potential {alpha * pot(np.sqrt(5.0))}\n")

    file.write("}\n")
    file.write("export -f runs\n")
    file.write("runs | nohup nice parallel -j 4 &\n")

    file.close()
