import csv
import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

if __name__ == "__main__":

    alpha = 1
    stats = "fermions"
    q = 96
    scale_factor = q**2 / 2 if stats == "fermions" else q
    file_name = "Fig11c"

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data/old/sp_data/q_{q:g}.txt"
    mb_data = f"/home/bart/DiagHam_latest/run/SFCI_data/{stats}_alpha_{alpha:g}/q_{q:g}/mb_ener_q_{q:g}.txt"

    t3hop, tism = [], []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            t3hop.append(float(row[0]))
            tism.append(float(row[1]))

    gap = []
    with open(mb_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            gap.append(float(row[1])*scale_factor)

    fig = plt.figure()
    ax1 = plt.subplot(111)

    sc = ax1.scatter(tism, gap, c=t3hop, zorder=5)
    ax1.set_xlim([-0.002, 0.016])
    # ax1.set_ylim([0.605, 0.618])  # bosons_alpha_0
    # ax1.set_ylim([1.520, 1.555])  # bosons_alpha_1
    ax1.set_ylim([0.482, 0.522])

    cbar = plt.colorbar(sc)
    cbar.set_label('$t_3$')
    plt.grid()

    # leg = ax1.legend(loc='best', title='$p$', ncol=3, handletextpad=0.5, handlelength=0, labelspacing=0,
    #                  borderpad=0.35, framealpha=1, markerscale=0.8, fontsize=10, columnspacing=0.5)

    ax1.set_xlabel('$\\langle \\mathcal{T} \\rangle$')
    scale_label = "q^2" if stats == "fermions" else "q"
    ax1.set_ylabel(f'${scale_label}\\Delta_\\mathrm{{m.b.}}$')
    # ax1.set_title(title)
    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data/plots/{file_name}.png", bbox_inches='tight', dpi=300)

    plt.show()
