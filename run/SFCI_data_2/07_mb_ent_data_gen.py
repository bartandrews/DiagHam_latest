import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import subprocess


def can_convert_to_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


if __name__ == "__main__":

    stats = "bosons"  # "fermions" or "bosons"
    alpha = 0
    q = 144
    ts = np.linspace(-0.25, 0.25, 11)
    if stats == "bosons":
        exp_count = 660
    else:
        exp_count = 2730

    os.chdir(f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats}_alpha_{alpha}/q_{q:g}/ent")
    out_file = open(f"mb_ent_q_{q:g}.txt", "w")

    threshold, threshold_count = 123, 0  # initial ansatz
    for t6hop in ts:
        for t9hop in ts:
            t3hop = (-9 * t6hop - 16 * t9hop - 1) / 4  # quartic plane

            # extract the ent spec image file
            for file in os.listdir("."):
                if f"t6_{t6hop:.2f}_t9_{t9hop:.2f}_" in file and file.endswith(".parentspec.ps"):
                    spec_image_file = file

            # extract the ent spec file
            for file in os.listdir("."):
                if f"t6_{t6hop:.2f}_t9_{t9hop:.2f}_" in file and file.endswith(".parentspec"):
                    spec_file = file

            # compute the spectral gap
            ents = []
            with open(spec_file, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter=' ')
                for j, row in enumerate(plots):
                    if can_convert_to_float(row[0]):  # if value is a number
                        ents.append(float(row[5]))

            # # compute the new threshold count
            # new_threshold_count = 0
            # with open(spec_file, 'r') as csvfile:
            #     plots = csv.reader(csvfile, delimiter=' ')
            #     for j, row in enumerate(plots):
            #         if can_convert_to_float(row[0]):  # if value is a number
            #             ent = float(row[5])
            #             if ent < threshold:
            #                 new_threshold_count += 1
            # print(f"threshold count = {new_threshold_count}")
            #
            # if new_threshold_count != threshold_count or new_threshold_count == 0:
            #     # plot the ent spec
            #     os.system(f"gmt psconvert {spec_image_file} -Tg")
            #     plt.imshow(mpimg.imread(spec_image_file.replace(".ps", ".png")))
            #     plt.axis('off')
            #     plt.show()
            #
            #     # determine the threshold
            #     threshold = float(input("Enter an entanglement threshold (0 if no gap present): "))
            #
            #     # compute the actual new threshold count
            #     new_threshold_count = 0
            #     with open(spec_file, 'r') as csvfile:
            #         plots = csv.reader(csvfile, delimiter=' ')
            #         for j, row in enumerate(plots):
            #             if can_convert_to_float(row[0]):  # if value is a number
            #                 ent = float(row[5])
            #                 if ent < threshold:
            #                     new_threshold_count += 1
            #     print(f"threshold count = {new_threshold_count}")
            #
            #     # reset the threshold count
            #     threshold_count = new_threshold_count

            if threshold == 0:
                ent_gap = np.nan
                print("ent_gap = ", ent_gap)
            elif threshold == 123:
                ents = np.sort(ents)
                ent_gap = ents[exp_count] - ents[exp_count-1]
            else:
                # find the entanglement gap
                FindEntanglementGap = "/home/bart/DiagHam_latest/scripts_bart/entanglement/FindEntanglementGap.sh "
                command = FindEntanglementGap+spec_image_file.replace(".ps", "")+f" {threshold}"
                ent_gap = float(subprocess.check_output(command, shell=True))
                print("ent_gap = ", float(ent_gap))

            out_file.write(f"{t6hop:.2f}\t{t9hop:.2f}\t{ent_gap}\n")

    out_file.close()
