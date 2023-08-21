import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

if __name__ == "__main__":

    stats = "bosons"
    alpha = 0
    q = 16
    file_name = f"laughlin_2d"

    if stats == "bosons":
        scale_factor = q
        title_str = "bosons"
    else:  # fermions
        scale_factor = q**2
        title_str = "fermions"

    sp_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_{q:g}.txt"
    mb_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats}_alpha_{alpha:g}/q_{q:g}/mb_ener_q_{q:g}.txt"
    mb_data_ent = f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats}_alpha_{alpha:g}/q_{q:g}/ent/mb_ent_q_{q:g}.txt"

    # single-particle quantities
    t6hop, tism, dism, berry_fluc, fs_fluc, gap_width = [], [], [], [], [], []
    with open(sp_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[1]) == 0:
                t6hop.append(float(row[0]))
                tism.append(float(row[2]))
                dism.append(float(row[3]))
                berry_fluc.append(float(row[4]))
                fs_fluc.append(float(row[5]))
                gap_width.append(float(row[6]))

    # many-body quantities
    gap, ent = [], []
    with open(mb_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[1]) == 0:
                gap.append(float(row[2])*scale_factor)
    with open(mb_data_ent, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[1]) == 0:
                ent.append(float(row[2]))

    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(6, 2, hspace=0, wspace=0.5, height_ratios=[3, 1, 1, 1, 1, 1])

    ##########
    # bosons #
    ##########

    ax0 = plt.subplot(gs[0])
    ax0.set_title(f"{title_str}, $n_\\phi=1/{q}$")
    #
    ax0_color = 'C0'
    ax0.tick_params('x', direction='in', bottom=True)
    ax0.plot(t6hop, gap, '.-', zorder=5, color=ax0_color)
    scale_label = "q^2" if stats == "fermions" else "q"
    ax0.set_ylabel(f'${scale_label}\\Delta_\\mathrm{{m.b.}}$', color=ax0_color)
    ax0.tick_params(axis='y', labelcolor=ax0_color)
    ax0.axvline(t6hop[np.argmax(gap)], ls='--', color=ax0_color)
    #
    ax0_2_color = 'C3'
    ax0_2 = ax0.twinx()
    ax0_2.plot(t6hop, ent, '.-', zorder=5, color=ax0_2_color)
    ax0_2.set_ylabel(f'$\\Delta_\\xi$', color=ax0_2_color)
    ax0_2.tick_params(axis='y', labelcolor=ax0_2_color)
    if np.max(ent) < 1000:
        ax0.axvline(t6hop[np.argmax(ent)], ls='--', color=ax0_2_color)

    ax2 = plt.subplot(gs[2], sharex=ax0)
    ax2.tick_params('x', direction='in', bottom=True)
    ax2.plot(t6hop, tism, '.-', zorder=5, c='k')
    ax2.set_ylabel('$\\langle \\mathcal{T} \\rangle$')
    ax2.axvline(t6hop[np.argmin(tism)], c='k', ls='--')

    ax4 = plt.subplot(gs[4], sharex=ax0)
    ax4.tick_params('x', direction='in', bottom=True)
    ax4.plot(t6hop, dism, '.-', zorder=5, c='k')
    ax4.set_ylabel('$\\langle \\mathcal{D} \\rangle$')
    ax4.axvline(t6hop[np.argmin(dism)], c='k', ls='--')

    ax6 = plt.subplot(gs[6], sharex=ax0)
    ax6.tick_params('x', direction='in', bottom=True)
    ax6.plot(t6hop, fs_fluc, '.-', zorder=5, c='k')
    ax6.set_ylabel('$\\log(\\sigma_g)$')
    # ax6.axvline(t6hop[np.argmin(fs_fluc)], c='k', ls='--')

    ax8 = plt.subplot(gs[8], sharex=ax0)
    ax8.tick_params('x', direction='in', bottom=True)
    ax8.plot(t6hop, berry_fluc, '.-', zorder=5, c='k')
    ax8.set_ylabel('$\\log(\\sigma_\\mathcal{B})$')
    # ax8.axvline(t6hop[np.argmin(berry_fluc)], c='k', ls='--')

    ax10 = plt.subplot(gs[10], sharex=ax0)
    ax10.plot(t6hop, gap_width, '.-', zorder=5, c='k')
    ax10.set_ylabel('$\\log(\\Delta / W)$')
    # ax10.axvline(t6hop[np.argmax(gap_width)], c='k', ls='--')

    ax10.set_xlabel('$t_6$')

    ############
    # fermions #
    ############

    ax1 = plt.subplot(gs[1])
    ax1.set_title(f"{title_str}, $n_\\phi=1/{q}$")
    #
    ax1_color = 'C0'
    ax1.tick_params('x', direction='in', bottom=True)
    ax1.plot(t6hop, gap, '.-', zorder=5, color=ax1_color)
    scale_label = "q^2" if stats == "fermions" else "q"
    ax1.set_ylabel(f'${scale_label}\\Delta_\\mathrm{{m.b.}}$', color=ax1_color)
    ax1.tick_params(axis='y', labelcolor=ax1_color)
    ax1.axvline(t6hop[np.argmax(gap)], ls='--', color=ax1_color)
    #
    ax1_2_color = 'C3'
    ax1_2 = ax1.twinx()
    ax1_2.plot(t6hop, ent, '.-', zorder=5, color=ax1_2_color)
    ax1_2.set_ylabel(f'$\\Delta_\\xi$', color=ax1_2_color)
    ax1_2.tick_params(axis='y', labelcolor=ax1_2_color)
    if np.max(ent) < 1000:
        ax1.axvline(t6hop[np.argmax(ent)], ls='--', color=ax1_2_color)

    ax3 = plt.subplot(gs[3], sharex=ax1)
    ax3.tick_params('x', direction='in', bottom=True)
    ax3.plot(t6hop, tism, '.-', zorder=5, c='k')
    ax3.set_ylabel('$\\langle \\mathcal{T} \\rangle$')
    ax3.axvline(t6hop[np.argmin(tism)], c='k', ls='--')

    ax5 = plt.subplot(gs[5], sharex=ax1)
    ax5.tick_params('x', direction='in', bottom=True)
    ax5.plot(t6hop, dism, '.-', zorder=5, c='k')
    ax5.set_ylabel('$\\langle \\mathcal{D} \\rangle$')
    ax5.axvline(t6hop[np.argmin(dism)], c='k', ls='--')

    ax7 = plt.subplot(gs[7], sharex=ax1)
    ax7.tick_params('x', direction='in', bottom=True)
    ax7.plot(t6hop, fs_fluc, '.-', zorder=5, c='k')
    ax7.set_ylabel('$\\log(\\sigma_g)$')
    # ax7.axvline(t6hop[np.argmin(fs_fluc)], c='k', ls='--')

    ax9 = plt.subplot(gs[9], sharex=ax1)
    ax9.tick_params('x', direction='in', bottom=True)
    ax9.plot(t6hop, berry_fluc, '.-', zorder=5, c='k')
    ax9.set_ylabel('$\\log(\\sigma_\\mathcal{B})$')
    # ax9.axvline(t6hop[np.argmin(berry_fluc)], c='k', ls='--')

    ax11 = plt.subplot(gs[11], sharex=ax1)
    ax11.plot(t6hop, gap_width, '.-', zorder=5, c='k')
    ax11.set_ylabel('$\\log(\\Delta / W)$')
    # ax11.axvline(t6hop[np.argmax(gap_width)], c='k', ls='--')

    ax11.set_xlabel('$t_6$')

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/{file_name}.png", bbox_inches='tight', dpi=300)
    plt.show()
