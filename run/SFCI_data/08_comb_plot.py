import csv
import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

if __name__ == "__main__":

    alpha = 0
    stats = "bosons"
    q = 16
    scale_factor = q**2 if stats == "fermions" else q
    file_name = f"{stats}_alpha_{alpha}_q_{q}_comb"

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data/sp_data/q_{q:g}.txt"
    sp_data_bart = f"/home/bart/DiagHam_latest/run/SFCI_data/sp_data/q_{q:g}_bart.txt"
    mb_data = f"/home/bart/DiagHam_latest/run/SFCI_data/{stats}_alpha_{alpha:g}/q_{q:g}/mb_ener_q_{q:g}.txt"

    t3hop, tism, gap_width, berry_fluc = [], [], [], []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            t3hop.append(float(row[0]))
            tism.append(float(row[3]))
    with open(sp_data_bart, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            gap_width.append(float(row[1]))
            berry_fluc.append(float(row[2]))

    gap = []
    with open(mb_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            gap.append(float(row[1])*scale_factor)

    fig = plt.figure()

    color = 'C0'
    ax1 = plt.subplot(111)
    ax1.plot(t3hop, gap, '.-', zorder=5, color=color, marker="x")
    ax1.set_xlabel('$t_3$')
    scale_label = "q^2" if stats == "fermions" else "q"
    ax1.set_ylabel(f'${scale_label}\\Delta_\\mathrm{{m.b.}}$', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.axvline(t3hop[np.argmax(gap)], c='k', ls='--')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'C2'
    # ax2.plot(t3hop, gap_width, '.-', zorder=5, color=color, label="$\\Delta / W$", marker="o")
    ax2.plot(t3hop, berry_fluc, '.-', zorder=5, color=color, label="$\\sigma_\\mathcal{B} / \\langle \\mathcal{B} \\rangle$", marker="^")
    # ax2.plot(t3hop, tism, '.-', zorder=5, color=color, label="$\\langle \\mathcal{T} \\rangle$", marker="s")
    ax2.tick_params(axis='y', labelcolor=color)
    # ax1.set_xlim([-0.002, 0.016])
    # ax1.set_ylim([0.605, 0.618])  # bosons_alpha_0
    # ax1.set_ylim([1.520, 1.555])  # bosons_alpha_1

    # cbar = plt.colorbar(sc)
    # cbar.set_label('$t_3$')

    leg = ax2.legend(loc='best', ncol=3, handletextpad=0.5, handlelength=0, labelspacing=0,
                     borderpad=0.35, framealpha=1, markerscale=0.8, fontsize=10, columnspacing=0.5)

    # ax1.set_title(title)
    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data/plots/{file_name}.png", bbox_inches='tight', dpi=300)

    plt.show()
