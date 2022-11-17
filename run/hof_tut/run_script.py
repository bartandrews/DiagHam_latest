import os

if __name__ == "__main__":

    executable = "/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel"

    # fermions with filling fraction 1/3
    os.makedirs("fermions", exist_ok=True)  # make the parent directory
    os.chdir("fermions")  # enter the parent directory
    for q in range(3, 7):  # loop over integer q from 3 to 6
        x, y = 3, 6  # number of MUCs
        os.makedirs(f"x_{x}_y_{y}_X_{q}_Y_3", exist_ok=True)  # make the directory
        os.chdir(f"x_{x}_y_{y}_X_{q}_Y_3")  # enter the directory
        flags = f" -p 6 -x {x} -y {y} -X {q} -Y 3 -m 10000 -S --processors 10 -n 1 --lanczos-precision 1e-10"
        command = executable+flags  # define the command
        os.system(command)  # execute the command
        os.chdir("..")  # leave the directory
