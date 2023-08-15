import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

if __name__ == "__main__":

    alpha = 0
    stats = "bosons"
    q = 16
    scale_factor = q**2 if stats == "fermions" else q
    file_name = f"{stats}_alpha_{alpha}_q_{q}_comb"

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data/sp_data/q_{q:g}.txt"
    mb_data = f"/home/bart/DiagHam_latest/run/SFCI_data/{stats}_alpha_{alpha:g}/q_{q:g}/mb_ener_q_{q:g}.txt"
    mb_data_ent = f"/home/bart/DiagHam_latest/run/SFCI_data/{stats}_alpha_{alpha:g}/q_{q:g}/ent/mb_ent_q_{q:g}.txt"

    # single-particle quantities
    t3hop, tism, dism, berry_fluc, fs_fluc, gap_width = [], [], [], [], [], []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            t3hop.append(float(row[0]))
            tism.append(float(row[1]))
            dism.append(float(row[2]))
            berry_fluc.append(float(row[3]))
            fs_fluc.append(float(row[4]))
            gap_width.append(float(row[5]))

    # many-body quantities
    gap, ent = [], []
    with open(mb_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            gap.append(float(row[1])*scale_factor)
    with open(mb_data_ent, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            ent.append(float(row[1]))

    fig = plt.figure(figsize=(6, 9))
    gs = gridspec.GridSpec(6, 1, hspace=0, height_ratios=[3, 1, 1, 1, 1, 1])

    ax1 = plt.subplot(gs[0])
    ax1.set_title(f"{stats}, $n_\\phi=1/{q}$, $\\alpha={alpha}$")
    #
    ax1_color = 'C0'
    ax1.tick_params('x', direction='in', bottom=True)
    ax1.plot(t3hop, gap, '.-', zorder=5, color=ax1_color)
    scale_label = "q^2" if stats == "fermions" else "q"
    ax1.set_ylabel(f'${scale_label}\\Delta_\\mathrm{{m.b.}}$', color=ax1_color)
    ax1.tick_params(axis='y', labelcolor=ax1_color)
    ax1.axvline(t3hop[np.argmax(gap)], ls='--', color=ax1_color)
    #
    ax1_2_color = 'C3'
    ax1_2 = ax1.twinx()
    ax1_2.plot(t3hop, ent, '.-', zorder=5, color=ax1_2_color)
    ax1_2.set_ylabel(f'$\\Delta_\\xi$', color=ax1_2_color)
    ax1_2.tick_params(axis='y', labelcolor=ax1_2_color)
    ax1.axvline(t3hop[np.argmax(ent)], ls='--', color=ax1_2_color)
    ax1.set_title(f"{stats}, $n_\\phi=1/{q}$, $\\alpha={alpha}$")

    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.tick_params('x', direction='in', bottom=True)
    ax2.plot(t3hop, tism, '.-', zorder=5, c='k')
    ax2.set_ylabel('$\\langle \\mathcal{T} \\rangle$')
    ax2.axvline(t3hop[np.argmin(tism)], c='k', ls='--')

    ax3 = plt.subplot(gs[2], sharex=ax1)
    ax3.tick_params('x', direction='in', bottom=True)
    ax3.plot(t3hop, dism, '.-', zorder=5, c='k')
    ax3.set_ylabel('$\\langle \\mathcal{D} \\rangle$')
    ax3.axvline(t3hop[np.argmin(dism)], c='k', ls='--')

    ax4 = plt.subplot(gs[3], sharex=ax1)
    ax4.plot(t3hop, gap_width, '.-', zorder=5, c='k')
    ax4.set_ylabel('$\\log(\\Delta / W)$')
    ax4.axvline(t3hop[np.argmax(gap_width)], c='k', ls='--')

    ax5 = plt.subplot(gs[4], sharex=ax1)
    ax5.tick_params('x', direction='in', bottom=True)
    ax5.plot(t3hop, berry_fluc, '.-', zorder=5, c='k')
    ax5.set_ylabel('$\\log(\\sigma_\\mathcal{B})$')
    # ax5.axvline(t3hop[np.argmin(berry_fluc)], c='k', ls='--')

    ax6 = plt.subplot(gs[5], sharex=ax1)
    ax6.tick_params('x', direction='in', bottom=True)
    ax6.plot(t3hop, fs_fluc, '.-', zorder=5, c='k')
    ax6.set_xlabel('$t_3$')
    ax6.set_ylabel('$\\log(\\sigma_g)$')
    # ax6.axvline(t3hop[np.argmin(fs_fluc)], c='k', ls='--')

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data/plots/{file_name}.png", bbox_inches='tight', dpi=300)
    plt.show()
