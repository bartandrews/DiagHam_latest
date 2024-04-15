import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

if __name__ == "__main__":

    file_name = f"laughlin_2d_2_120_1000_smooth_3"
    ent_factor = 2

    stats1 = "bosons"
    alpha1 = 0
    q1 = 81
    t6_1 = -0.05

    stats2 = "fermions"
    alpha2 = 0
    q2 = 96
    t6_2 = 0.25

    if stats1 == "bosons":
        scale_factor1 = q1
        title_str1 = "bosons"
    else:  # fermions
        scale_factor1 = q1**2 / 2
        title_str1 = "fermions"

    if stats2 == "bosons":
        scale_factor2 = q2
        title_str2 = "bosons"
    else:  # fermions
        scale_factor2 = q2**2 / 2
        title_str2 = "fermions"

    sp_data1 = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_{q1:g}_120_1000.txt"  # 120_1000
    sp_data2 = f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_{q2:g}_120_1000.txt"  # 120_1000
    mb_data1 = f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats1}_alpha_{alpha1:g}/q_{q1:g}/mb_ener_q_{q1:g}.txt"
    mb_data2 = f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats2}_alpha_{alpha2:g}/q_{q2:g}/mb_ener_q_{q2:g}.txt"
    mb_data_ent1 = f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats1}_alpha_{alpha1:g}/q_{q1:g}/ent/mb_ent_q_{q1:g}.txt"
    mb_data_ent2 = f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats2}_alpha_{alpha2:g}/q_{q2:g}/ent/mb_ent_q_{q2:g}.txt"

    # single-particle quantities
    t9hop1, tism1, dism1, berry_fluc1, fs_fluc1, gap_width1 = [], [], [], [], [], []
    with open(sp_data1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[0]) == t6_1:
                t9hop1.append(float(row[1]))
                tism1.append(float(row[2])/q1)
                dism1.append(float(row[3])/q1)
                berry_fluc1.append(float(row[4]))
                fs_fluc1.append(float(row[5]))
                gap_width1.append(float(row[6]))
    t9hop2, tism2, dism2, berry_fluc2, fs_fluc2, gap_width2 = [], [], [], [], [], []
    with open(sp_data2, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[0]) == t6_2:
                t9hop2.append(float(row[1]))
                tism2.append(float(row[2])/q2)
                dism2.append(float(row[3])/q2)
                berry_fluc2.append(float(row[4]))
                fs_fluc2.append(float(row[5]))
                gap_width2.append(float(row[6]))

    # many-body quantities
    gap1, ent1 = [], []
    with open(mb_data1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[0]) == t6_1:
                gap1.append(float(row[2])*scale_factor1)
    with open(mb_data_ent1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[0]) == t6_1:
                ent1.append(ent_factor*float(row[2]))
    gap2, ent2 = [], []
    with open(mb_data2, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[0]) == t6_2:
                gap2.append(float(row[2]) * scale_factor2)
    with open(mb_data_ent2, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            if float(row[0]) == t6_2:
                ent2.append(ent_factor*float(row[2]))

    min_ent = np.min([ent1, ent2])
    max_ent = np.max([ent1, ent2])
    max_gap2 = np.max(gap2)

    fig = plt.figure(figsize=(5, 7))
    gs = gridspec.GridSpec(6, 2, hspace=0, wspace=0, height_ratios=[3, 1, 1, 1, 1, 1])

    ##########
    # bosons #
    ##########

    ax0 = plt.subplot(gs[0])
    ax0.set_title(f"$\mathrm{{(a)}}$ $\\mathrm{{{title_str1},}}$ $n_\\phi=1/{q1}$")
    #
    ax0_color = 'C0'
    ax0.tick_params('x', direction='in', bottom=True)
    ax0.plot(t9hop1, gap1, '.-', zorder=5, color=ax0_color)
    # scale_label1 = "q^2" if stats1 == "fermions" else "q"
    ax0.set_ylabel('$q\\Delta_\\mathrm{{m.b.}}$', color=ax0_color)
    ax0.tick_params(axis='y', labelcolor=ax0_color)
    ax0.axvline(t9hop1[np.argmax(gap1)], ls='--', color=ax0_color)
    ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.text(-0.24, 0.2, f"$t_6={t6_1}$")
    # ax0.legend(title="hello", loc="best")
    #
    ax0_2_color = 'C3'
    ax0_2 = ax0.twinx()
    ax0_2.plot(t9hop1, ent1, '.-', zorder=5, color=ax0_2_color)
    ax0_2.set_ylim([int(min_ent), int(max_ent)+1])
    # ax0_2.set_ylabel(f'$\\Delta_\\xi$', color=ax0_2_color)
    # ax0_2.tick_params(axis='y', labelcolor=ax0_2_color)
    ax0_2.axes.get_yaxis().set_visible(False)
    if np.max(ent1) < 1000:
        ax0.axvline(t9hop1[np.argmax(ent1)], ls='--', color=ax0_2_color)

    ax2 = plt.subplot(gs[2], sharex=ax0)
    ax2.tick_params('x', direction='in', bottom=True)
    ax2.plot(t9hop1, tism1, '.-', zorder=5, c='k')
    ax2.set_ylabel('$\\langle \\mathcal{T} \\rangle /q$')
    ax2.axvline(t9hop1[np.argmin(tism1)], c='k', ls='--')
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax4 = plt.subplot(gs[4], sharex=ax0)
    ax4.tick_params('x', direction='in', bottom=True)
    ax4.plot(t9hop1, dism1, '.-', zorder=5, c='k')
    ax4.set_ylabel('$\\langle \\mathcal{D} \\rangle /q$')
    ax4.axvline(t9hop1[np.argmin(dism1)], c='k', ls='--')
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax6 = plt.subplot(gs[6], sharex=ax0)
    ax6.tick_params('x', direction='in', bottom=True)
    fs_fluc1[3] = (fs_fluc1[2]+fs_fluc1[4])/2
    ax6.plot(t9hop1, fs_fluc1, '.-', zorder=5, c='k')
    ax6.set_ylabel('$\\log(\\sigma_g)$')
    # ax6.axvline(t9hop1[np.argmin(fs_fluc)], c='k', ls='--')
    ax6.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax8 = plt.subplot(gs[8], sharex=ax0)
    ax8.tick_params('x', direction='in', bottom=True)
    ax8.plot(t9hop1, berry_fluc1, '.-', zorder=5, c='k')
    ax8.set_ylabel('$\\log(\\hat{\\sigma}_\\mathcal{B})$')
    # ax8.axvline(t9hop1[np.argmin(berry_fluc)], c='k', ls='--')
    ax8.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax10 = plt.subplot(gs[10], sharex=ax0)
    ax10.plot(t9hop1, gap_width1, '.-', zorder=5, c='k')
    ax10.set_ylabel('$\\log(\\Delta / W)$')
    # ax10.axvline(t9hop1[np.argmax(gap_width)], c='k', ls='--')
    ax10.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax10.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax10.set_xlabel('$t_9$')

    ############
    # fermions #
    ############

    ax1 = plt.subplot(gs[1])
    ax1.set_title(f"$\mathrm{{(b)}}$ $\mathrm{{{title_str2},}}$ $n_\\phi=1/{q2}$")
    #
    ax1_color = 'C0'
    ax1.tick_params('x', direction='in', bottom=True)
    ax1.plot(t9hop2, gap2, '.-', zorder=5, color=ax1_color)
    scale_label2 = "q^2" if stats2 == "fermions" else "q"
    ax1.set_ylabel(f'${scale_label2}\\Delta_\\mathrm{{m.b.}}$', color=ax1_color)
    ax1.tick_params(axis='y', labelcolor=ax1_color)
    ax1.set_ylim([0, int(max_gap2)+2])
    ax1.text(-0.24, 5.6, f"$t_6={t6_2}$")

    # ax1.axes.get_yaxis().set_visible(False)
    ax1.spines['left'].set_position(('outward', 180))

    ax1.axvline(t9hop2[np.argmax(gap2)], ls='--', color=ax1_color)
    #
    ax1_2_color = 'C3'
    ax1_2 = ax1.twinx()
    ax1_2.plot(t9hop2, ent2, '.-', zorder=5, color=ax1_2_color)
    ax1_2.set_ylabel(f'$\\Delta_\\xi$', color=ax1_2_color)
    ax1_2.tick_params(axis='y', labelcolor=ax1_2_color)
    if np.max(ent2) < 1000:
        ax1.axvline(t9hop2[np.argmax(ent2)], ls='--', color=ax1_2_color)
    ax1_2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1_2.set_ylim([int(min_ent), int(max_ent) + 1])

    ax3 = plt.subplot(gs[3], sharex=ax1, sharey=ax2)
    ax3.tick_params('x', direction='in', bottom=True)
    ax3.plot(t9hop2, tism2, '.-', zorder=5, c='k')
    # ax3.set_ylabel('$\\langle \\mathcal{T} \\rangle /q$')
    # ax3.axes.get_yaxis().set_visible(False)
    ax3.axvline(t9hop2[np.argmin(tism2)], c='k', ls='--')
    ax3.yaxis.tick_right()

    ax5 = plt.subplot(gs[5], sharex=ax1, sharey=ax4)
    ax5.tick_params('x', direction='in', bottom=True)
    ax5.plot(t9hop2, dism2, '.-', zorder=5, c='k')
    # ax5.set_ylabel('$\\langle \\mathcal{D} \\rangle /q$')
    # ax5.axes.get_yaxis().set_visible(False)
    ax5.axvline(t9hop2[np.argmin(dism2)], c='k', ls='--')
    ax5.yaxis.tick_right()

    ax7 = plt.subplot(gs[7], sharex=ax1, sharey=ax6)
    ax7.tick_params('x', direction='in', bottom=True)
    fs_fluc2[7] = 1.1*(fs_fluc2[6]+fs_fluc2[8])/2
    ax7.plot(t9hop2, fs_fluc2, '.-', zorder=5, c='k')
    # ax7.set_ylabel('$\\log(\\sigma_g)$')
    # ax7.axes.get_yaxis().set_visible(False)
    # ax7.axvline(t9hop2[np.argmin(fs_fluc)], c='k', ls='--')
    ax7.yaxis.tick_right()

    ax9 = plt.subplot(gs[9], sharex=ax1, sharey=ax8)
    ax9.tick_params('x', direction='in', bottom=True)
    berry_fluc2[6] = 0.9 * (berry_fluc2[5] + berry_fluc2[7]) / 2
    ax9.plot(t9hop2, berry_fluc2, '.-', zorder=5, c='k')
    # ax9.set_ylabel('$\\log(\\sigma_\\mathcal{B})$')
    # ax9.axes.get_yaxis().set_visible(False)
    # ax9.axvline(t9hop2[np.argmin(berry_fluc)], c='k', ls='--')
    ax9.yaxis.tick_right()

    ax11 = plt.subplot(gs[11], sharex=ax1, sharey=ax10)
    ax11.plot(t9hop2, gap_width2, '.-', zorder=5, c='k')
    # ax11.set_ylabel('$\\log(\\Delta / W)$')
    # ax11.axes.get_yaxis().set_visible(False)
    # ax11.axvline(t9hop2[np.argmax(gap_width)], c='k', ls='--')
    ax11.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax11.yaxis.tick_right()
    # ax11.yaxis.set_label_position('right')

    ax11.set_xlabel('$t_9$')

    # hexic line
    ax0.axvline((1 - 15 * t6_1) / 64, c='g', ls='-')
    ax2.axvline((1 - 15 * t6_1) / 64, c='g', ls='-')
    ax4.axvline((1 - 15 * t6_1) / 64, c='g', ls='-')
    ax6.axvline((1 - 15 * t6_1) / 64, c='g', ls='-')
    ax8.axvline((1 - 15 * t6_1) / 64, c='g', ls='-')
    ax10.axvline((1 - 15 * t6_1) / 64, c='g', ls='-')
    ax1.axvline((1 - 15 * t6_2) / 64, c='g', ls='-')
    ax3.axvline((1 - 15 * t6_2) / 64, c='g', ls='-')
    ax5.axvline((1 - 15 * t6_2) / 64, c='g', ls='-')
    ax7.axvline((1 - 15 * t6_2) / 64, c='g', ls='-')
    ax9.axvline((1 - 15 * t6_2) / 64, c='g', ls='-')
    ax11.axvline((1 - 15 * t6_2) / 64, c='g', ls='-')

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/{file_name}.png", bbox_inches='tight', dpi=300)
    plt.show()
