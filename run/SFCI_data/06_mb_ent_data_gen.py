import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import subprocess


if __name__ == "__main__":

    stats = "bosons"  # "fermions" or "bosons"
    q = 81

    for alpha in [1]:
        os.chdir(f"/home/bart/DiagHam_latest/run/SFCI_data/{stats}_alpha_{alpha}/q_{q:g}/ent")
        out_file = open(f"mb_ent_q_{q:g}.txt", "w")

        for t3hop in np.linspace(-0.25, 0, 26):

            # extract the ent spec image file
            for file in os.listdir("."):
                if f"t3_{t3hop:.2f}_" in file and file.endswith(".parentspec.ps"):
                    spec_image_file = file

            # plot the ent spec
            os.system(f"gmt psconvert {spec_image_file} -Tg")
            plt.imshow(mpimg.imread(spec_image_file.replace(".ps", ".png")))
            plt.axis('off')
            plt.show()

            # determine the threshold
            threshold = float(input("Enter an entanglement threshold: "))

            # find the entanglement gap
            FindEntanglementGap = "/home/bart/DiagHam_latest/scripts_bart/entanglement/FindEntanglementGap.sh "
            command = FindEntanglementGap+spec_image_file.replace(".ps", "")+f" {threshold}"
            ent_gap = float(subprocess.check_output(command, shell=True))
            print("ent_gap = ", float(ent_gap))

            out_file.write(f"{t3hop:.2f}\t{ent_gap}\n")

        out_file.close()
