import numpy as np
import csv
import matplotlib.pyplot as plt

if __name__ == "__main__":

    fig = plt.figure()
    ax1 = plt.subplot(111)

    sp_data = f"q_4.txt"
    t3 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            t3.append(float(row[0]))
            TISM.append(float(row[3]))
    ax1.plot(t3, TISM, '.-', label="1/4")

    sp_data = f"q_6.txt"
    t3 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            t3.append(float(row[0]))
            TISM.append(float(row[3]))
    ax1.plot(t3, TISM, '.-', label="1/6")

    sp_data = f"q_24.txt"
    t3 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            t3.append(float(row[0]))
            TISM.append(float(row[3]))
    ax1.plot(t3, TISM, '.-', label="1/24")

    sp_data = f"q_54.txt"
    t3 = []
    TISM = []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            t3.append(float(row[0]))
            TISM.append(float(row[3]))
    ax1.plot(t3, TISM, '.-', label="1/54")

    ax1.legend(title="$n_\phi$", loc="upper center")
    ax1.set_xlabel('$t_3$')
    ax1.set_ylabel('$\\langle \\mathcal{T} \\rangle$')

    plt.savefig(f"reproduction_Fig_3.png", bbox_inches='tight', dpi=300)
    plt.show()
