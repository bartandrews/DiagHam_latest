import numpy as np
import csv
import matplotlib.pyplot as plt

if __name__ == "__main__":

    fig = plt.figure()
    ax1 = plt.subplot(111)

    # sp_data = f"q_4_2.txt"
    # t6 = []
    # TISM = []
    # with open(sp_data, 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     for i, row in enumerate(plots):
    #         t6.append(float(row[0]))
    #         TISM.append(float(row[2]))
    # ax1.plot(t6, TISM, '.-', label="1/4")
    #
    # sp_data = f"q_16_2.txt"
    # t6 = []
    # TISM = []
    # with open(sp_data, 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     for i, row in enumerate(plots):
    #         t6.append(float(row[0]))
    #         TISM.append(float(row[2]))
    # ax1.plot(t6, TISM, '.-', label="1/16")

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_16.txt"
    t6 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            if float(row[0]) == float(row[1]):
                t6.append(float(row[0]))
                TISM.append(float(row[2]))
    ax1.plot(t6, TISM, '.-', label="1/16")

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_24.txt"
    t6 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            if float(row[0]) == float(row[1]):
                t6.append(float(row[0]))
                TISM.append(float(row[2]))
    ax1.plot(t6, TISM, '.-', label="1/24")

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_49.txt"
    t6 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            if float(row[0]) == float(row[1]):
                t6.append(float(row[0]))
                TISM.append(float(row[2]))
    ax1.plot(t6, TISM, '.-', label="1/49")

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_54.txt"
    t6 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            if float(row[0]) == float(row[1]):
                t6.append(float(row[0]))
                TISM.append(float(row[2]))
    ax1.plot(t6, TISM, '.-', label="1/54")

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_81.txt"
    t6 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            if float(row[0]) == float(row[1]):
                t6.append(float(row[0]))
                TISM.append(float(row[2]))
    ax1.plot(t6, TISM, '.-', label="1/81")

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_96.txt"
    t6 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            if float(row[0]) == float(row[1]):
                t6.append(float(row[0]))
                TISM.append(float(row[2]))
    ax1.plot(t6, TISM, '.-', label="1/96")

    ax1.legend(title="$n_\phi$", loc="best")
    ax1.set_xlabel('$t_6=t_9$')
    ax1.set_ylabel('$\\langle \\mathcal{T} \\rangle$')
    ax1.set_xlim([-0.25, 0.05])
    ax1.set_ylim([0.01, 0.015])
    ax1.axvline(1/79, c='g', ls='--')

    plt.savefig(f"diag_cross_section_new_Fig_3.png", bbox_inches='tight', dpi=300)
    plt.show()
