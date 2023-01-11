import numpy as np


def pot(distance):
    return np.exp(1-np.abs(distance)**4)


if __name__ == "__main__":

    program = "tune"  # "tune" or "turn"
    stats = "fermions"  # "fermions" or "bosons"
    FCIHofstadterModel = "~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel"

    file = open(f"02_script.sh", "w")
    file.write("#!/bin/bash\n\n")

    file.write("runs() {\n")
    if program == "tune":
        for alpha in np.linspace(0, 1, 11):
            if stats == "fermions":
                file.write(f"echo {FCIHofstadterModel} -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 "
                           f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                           f"--u-potential {pot(1)} --v-potential {alpha*pot(np.sqrt(2))} "
                           f"--v2-potential {alpha*pot(2)} --v3-potential {alpha*pot(np.sqrt(5))}\n")
            else:
                file.write(f"echo {FCIHofstadterModel} --boson -p 6 -X 16 -Y 12 -x 3 -y 4 -m 8000 "
                           f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                           f"--u-potential {pot(0)} --v-potential {pot(1)} "
                           f"--v2-potential {alpha*pot(np.sqrt(2))} --v3-potential {alpha*pot(2)} "
                           f"--v4-potential {alpha*pot(np.sqrt(5))}\n")
    else:
        if stats == "fermions":
            N = 6
            X, Y, x, y = 18, 9, 3, 6
            file.write(f"echo {FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--atan_0-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--atan_1_4-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--atan_1_3-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--atan_1_2-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--atan_2_3-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--atan_3_4-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--atan_1-potential {1}\n")
        else:
            N = 6
            X, Y, x, y = 4, 2, 2, 6
            file.write(f"echo {FCIHofstadterModel} --boson -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--u-potential 2.72 --atan_0-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} --boson -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--u-potential 2.72 --atan_1_4-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} --boson -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--u-potential 2.72 --atan_1_3-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} --boson -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--u-potential 2.72 --atan_1_2-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} --boson -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--u-potential 2.72 --atan_2_3-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} --boson -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--u-potential 2.72 --atan_3_4-potential {1}\n")
            file.write(f"echo {FCIHofstadterModel} --boson -p {N} -X {X} -Y {Y} -x {x} -y {y} -m 8000 "
                       f"-S --processors 4 -n 20 --lanczos-precision 1e-10 "
                       f"--u-potential 2.72 --atan_1-potential {1}\n")

    file.write("}\n")
    file.write("export -f runs\n")
    file.write("runs | nohup nice parallel -j 4 > /dev/null &\n")

    file.close()
