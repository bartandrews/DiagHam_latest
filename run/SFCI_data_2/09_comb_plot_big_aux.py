import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def hexic_line(t6_val):
    return (1/64) * (1 - 15*t6_val)


if __name__ == "__main__":

    stats = "fermions"
    alpha = 0
    ts = np.linspace(-0.25, 0.25, 11)
    file_name = f"{stats}_alpha_{alpha}_big_aux"

    if stats == "bosons":
        qs = [100, 121, 144, 169, 196, 225, 256]
        scale_power = 1  # power of q for many-body gap scaling
        scale_divisor = 1  # divisor for many-body gap
        cbtitle = "$q \\Delta_\\mathrm{m.b.}$"
        cbtitle2 = "$q \\delta / 10^{-1}$"
        title_str = f"Bosons, $\\nu=1/2$, $V_{{ij}} = (1-\\alpha)\\delta_{{ij}} + \\alpha e^{{-|r_{{ij}}|^4}}$, $\\alpha={alpha}$"
        ent_factor = 2
    else:  # fermions
        qs = [126, 140, 150, 160, 168, 176, 187]
        scale_power = 2  # power of q for many-body gap scaling
        scale_divisor = 2  # divisor for many-body gap
        cbtitle = "$q^2 \\Delta_\\mathrm{m.b.}$"
        cbtitle2 = "$q^2 \\delta $"
        title_str = f"Fermions, $\\nu=1/3$, $V_{{ij}} = (1-\\alpha)\\delta_{{\\langle ij \\rangle}} + \\alpha e^{{-|r_{{ij}}|^4}}$, $\\alpha={alpha}$"
        ent_factor = 2

    sp_data, mb_data, mb_ent_data = [], [], []
    sp_data.append(f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_96.txt")
    for i, _ in enumerate(qs):
        # sp_data.append(f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_{qs[i]:g}.txt")
        mb_data.append(f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats}_alpha_{alpha:g}/q_{qs[i]:g}/mb_ener_q_{qs[i]:g}.txt")
        mb_ent_data.append(f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats}_alpha_{alpha:g}/q_{qs[i]:g}/ent/mb_ent_q_{qs[i]:g}.txt")

    ts_min = np.min(ts)
    ts_max = np.max(ts)
    ts_len = len(ts)
    ts_step = (ts_max - ts_min) / (ts_len - 1)

    plot_min = ts_min - ts_step / 2
    plot_max = ts_max + ts_step / 2
    plot_ts = np.linspace(plot_min, plot_max, 100)

    tisms = np.zeros((len(qs), len(ts), len(ts)))
    gaps = np.zeros((len(qs), len(ts), len(ts)))
    spread = np.zeros((len(qs), len(ts), len(ts)))
    ent_gaps = np.zeros((len(qs), len(ts), len(ts)))

    # single-particle quantities
    for iq, q in enumerate(qs):
        with open(sp_data[0], 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for i, row in enumerate(plots):
                it6 = int(i / len(ts))
                it9 = i % len(ts)
                tisms[iq, it6, it9] = float(row[2])/96
    min_tism, max_tism = np.nanmin(tisms), np.nanmax(tisms)
    # many-body quantities
    for iq, q in enumerate(qs):
        with open(mb_data[iq], 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for i, row in enumerate(plots):
                it6 = int(i / len(ts))
                it9 = i % len(ts)
                gaps[iq, it6, it9] = q**scale_power * float(row[2]) / scale_divisor
                if stats == "bosons":
                    spread[iq, it6, it9] = (q**scale_power * float(row[3]) / scale_divisor)*10
                else:
                    spread[iq, it6, it9] = (q ** scale_power * float(row[3]) / scale_divisor)
                    spread[iq, it6, it9] = (q ** scale_power * float(row[3]) / scale_divisor)
    if stats == "bosons":
        min_gap, max_gap = np.nanmin(gaps), np.nanmax(gaps)
    else:
        min_gap, max_gap = np.nanmin(gaps), 2  # np.nanmax(gaps[:1])  # [:2]
    min_spread, max_spread = np.nanmin(spread), np.nanmax(spread)
    for iq, _ in enumerate(qs):
        with open(mb_ent_data[iq], 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for i, row in enumerate(plots):
                it6 = int(i / len(ts))
                it9 = i % len(ts)
                ent_gaps[iq, it6, it9] = ent_factor*float(row[2])
    min_ent_gap, max_ent_gap = np.nanmin(ent_gaps), np.nanmax(ent_gaps)

    fig = plt.figure(figsize=(14, 8))
    # fig.suptitle(f'{title_str}')
    gs = gridspec.GridSpec(4, 7, wspace=0, hspace=0.2)
    gs.update(top=0.92)
    axlabel_fontsize = 12
    title_fontsize = 13

    ########
    # TISM #
    ########

    ax0 = plt.subplot(gs[0])
    iq = 0
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if tisms[iq, i, j] > 1000:
                tisms[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax0.imshow(tisms[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_tism, vmax=max_tism)
    ax0.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax0.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\langle \\mathcal{T} \\rangle$', fontsize=axlabel_fontsize)
    ax0.grid()
    # ax0.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    ax0.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    # ax0.tick_params(labelleft=False)
    ax0.tick_params(labelbottom=False)
    ax0.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax1 = plt.subplot(gs[1])
    iq = 1
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if tisms[iq, i, j] > 1000:
                tisms[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax1.imshow(tisms[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_tism, vmax=max_tism)
    ax1.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax1.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\langle \\mathcal{T} \\rangle$', fontsize=axlabel_fontsize)
    ax1.grid()
    # ax1.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax1.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax1.tick_params(labelleft=False)
    ax1.tick_params(labelbottom=False)
    ax1.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax2 = plt.subplot(gs[2])
    iq = 2
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if tisms[iq, i, j] > 1000:
                tisms[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax2.imshow(tisms[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_tism, vmax=max_tism)
    ax2.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax2.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\langle \\mathcal{T} \\rangle$', fontsize=axlabel_fontsize)
    ax2.grid()
    # ax2.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax2.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax2.tick_params(labelleft=False)
    ax2.tick_params(labelbottom=False)
    ax2.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax3 = plt.subplot(gs[3])
    iq = 3
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if tisms[iq, i, j] > 1000:
                tisms[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax3.imshow(tisms[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_tism, vmax=max_tism)
    ax3.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax3.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\langle \\mathcal{T} \\rangle$', fontsize=axlabel_fontsize)
    ax3.grid()
    # ax3.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax3.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax3.tick_params(labelleft=False)
    ax3.tick_params(labelbottom=False)
    ax3.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax4 = plt.subplot(gs[4])
    iq = 4
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if tisms[iq, i, j] > 1000:
                tisms[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax4.imshow(tisms[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_tism, vmax=max_tism)
    ax4.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax4.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\langle \\mathcal{T} \\rangle$', fontsize=axlabel_fontsize)
    ax4.grid()
    # ax4.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax4.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax4.tick_params(labelleft=False)
    ax4.tick_params(labelbottom=False)
    ax4.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax5 = plt.subplot(gs[5])
    iq = 5
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if tisms[iq, i, j] > 1000:
                tisms[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax5.imshow(tisms[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_tism, vmax=max_tism)
    ax5.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax5.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\langle \\mathcal{T} \\rangle$', fontsize=axlabel_fontsize)
    ax5.grid()
    # ax5.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax5.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax5.tick_params(labelleft=False)
    ax5.tick_params(labelbottom=False)
    ax5.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax6 = plt.subplot(gs[6])
    iq = 6
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if tisms[iq, i, j] > 1000:
                tisms[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax6.imshow(tisms[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_tism, vmax=max_tism)
    ax6.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax6.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # #cbar = plt.colorbar(sc)
    cb_ax = fig.add_axes([.9, .744, .00429, .176])
    fmt = lambda x, pos: '${:g}$'.format(x)
    cbar = fig.colorbar(sc, orientation='vertical', cax=cb_ax, format=FuncFormatter(fmt))
    cbar.set_label('$\\langle \\mathcal{T} \\rangle / q$', fontsize=axlabel_fontsize)

    ax6.grid()
    # ax6.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax6.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax6.tick_params(labelleft=False)
    ax6.tick_params(labelbottom=False)
    ax6.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax6.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax6.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # if stats == "fermions":
    #     ax6.plot(-0.15, 0.05, 'x', c='r')

    #################
    # Many-body gap #
    #################

    ax7 = plt.subplot(gs[7])
    iq = 0
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax7.imshow(gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_gap, vmax=max_gap)
    ax7.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax7.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax7.grid()
    # ax7.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    ax7.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    # ax7.tick_params(labelleft=False)
    ax7.tick_params(labelbottom=False)
    # ax7.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax7.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax7.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax8 = plt.subplot(gs[8])
    iq = 1
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax8.imshow(gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_gap, vmax=max_gap)
    ax8.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax8.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax8.grid()
    # ax8.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax8.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax8.tick_params(labelleft=False)
    ax8.tick_params(labelbottom=False)
    # ax8.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax8.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax8.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax9 = plt.subplot(gs[9])
    iq = 2
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax9.imshow(gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_gap, vmax=max_gap)
    ax9.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax9.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax9.grid()
    # ax9.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax9.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax9.tick_params(labelleft=False)
    ax9.tick_params(labelbottom=False)
    # ax9.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax9.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax9.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax10 = plt.subplot(gs[10])
    iq = 3
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax10.imshow(gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_gap, vmax=max_gap)
    ax10.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax10.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax10.grid()
    # ax10.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax10.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax10.tick_params(labelleft=False)
    ax10.tick_params(labelbottom=False)
    # ax10.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax10.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax10.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax11 = plt.subplot(gs[11])
    iq = 4
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax11.imshow(gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_gap, vmax=max_gap)
    ax11.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax11.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax11.grid()
    # ax11.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax11.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax11.tick_params(labelleft=False)
    ax11.tick_params(labelbottom=False)
    # ax11.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax11.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax11.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax12 = plt.subplot(gs[12])
    iq = 5
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax12.imshow(gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_gap, vmax=max_gap)
    ax12.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax12.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax12.grid()
    # ax12.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax12.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax12.tick_params(labelleft=False)
    ax12.tick_params(labelbottom=False)
    # ax12.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax12.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax12.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax13 = plt.subplot(gs[13])
    iq = 6
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax13.imshow(gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_gap, vmax=max_gap)
    ax13.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax13.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    cb_ax = fig.add_axes([.9, .532, .00429, .176])
    if stats == "bosons":
        cbar = fig.colorbar(sc, orientation='vertical', cax=cb_ax, format=FuncFormatter(fmt))
    else:
        cbar = fig.colorbar(sc, orientation='vertical', cax=cb_ax, extend='max', format=FuncFormatter(fmt))
        cbar.cmap.set_over('yellow')
    cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax13.grid()
    # ax13.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax13.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax13.tick_params(labelleft=False)
    ax13.tick_params(labelbottom=False)
    # ax13.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax13.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax13.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # if stats == "fermions":
    #     ax13.plot(-0.15, 0.05, 'x', c='r')

    #####################
    # degeneracy spread #
    #####################

    ax14 = plt.subplot(gs[14])
    iq = 0
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax14.imshow(spread[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_spread, vmax=max_spread)
    ax14.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax14.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax14.grid()
    # ax14.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    ax14.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    # ax14.tick_params(labelleft=False)
    ax14.tick_params(labelbottom=False)
    # ax14.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax14.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax14.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax15 = plt.subplot(gs[15])
    iq = 1
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax15.imshow(spread[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_spread, vmax=max_spread)
    ax15.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax15.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax15.grid()
    # ax15.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax15.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    ax15.tick_params(labelleft=False)
    ax15.tick_params(labelbottom=False)
    # ax15.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax15.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax15.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax16 = plt.subplot(gs[16])
    iq = 2
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax16.imshow(spread[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_spread, vmax=max_spread)
    ax16.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax16.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax16.grid()
    # ax16.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax16.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    ax16.tick_params(labelleft=False)
    ax16.tick_params(labelbottom=False)
    # ax16.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax16.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax16.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax17 = plt.subplot(gs[17])
    iq = 3
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax17.imshow(spread[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_spread, vmax=max_spread)
    ax17.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax17.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax17.grid()
    # ax17.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax17.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    ax17.tick_params(labelleft=False)
    ax17.tick_params(labelbottom=False)
    # ax17.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax17.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax17.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax18 = plt.subplot(gs[18])
    iq = 4
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax18.imshow(spread[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_spread, vmax=max_spread)
    ax18.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax18.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax18.grid()
    # ax18.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax18.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    ax18.tick_params(labelleft=False)
    ax18.tick_params(labelbottom=False)
    # ax18.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax18.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax18.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax19 = plt.subplot(gs[19])
    iq = 5
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax19.imshow(spread[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_spread, vmax=max_spread)
    ax19.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax19.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label(f'{cbtitle}', fontsize=axlabel_fontsize)
    ax19.grid()
    # ax19.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax19.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    ax19.tick_params(labelleft=False)
    ax19.tick_params(labelbottom=False)
    # ax19.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax19.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax19.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax20 = plt.subplot(gs[20])
    iq = 6
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax20.imshow(spread[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_spread, vmax=max_spread)
    ax20.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax20.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    cb_ax = fig.add_axes([.9, .321, .00429, .176])
    cbar = fig.colorbar(sc, orientation='vertical', cax=cb_ax, format=FuncFormatter(fmt))
    cbar.set_label(f'{cbtitle2}', fontsize=axlabel_fontsize)
    ax20.grid()
    # ax20.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax20.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    ax20.tick_params(labelleft=False)
    ax20.tick_params(labelbottom=False)
    # ax20.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax20.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax20.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # if stats == "fermions":
    #     ax20.plot(-0.15, 0.05, 'x', c='r')

    ####################
    # Entanglement gap #
    ####################

    ax21 = plt.subplot(gs[21])
    iq = 0
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if ent_gaps[iq, i, j] > 1000:
                ent_gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax21.imshow(ent_gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_ent_gap, vmax=max_ent_gap)
    ax21.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax21.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\Delta_\\xi$', fontsize=axlabel_fontsize)
    ax21.grid()
    ax21.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    ax21.set_ylabel('$t_9$', fontsize=axlabel_fontsize, labelpad=-8)
    # ax21.tick_params(labelleft=False)
    # ax21.tick_params(labelbottom=False)
    # ax21.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax21.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax21.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax22 = plt.subplot(gs[22])
    iq = 1
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if ent_gaps[iq, i, j] > 1000:
                ent_gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax22.imshow(ent_gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_ent_gap, vmax=max_ent_gap)
    ax22.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax22.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\Delta_\\xi$', fontsize=axlabel_fontsize)
    ax22.grid()
    ax22.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax22.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax22.tick_params(labelleft=False)
    # ax22.tick_params(labelbottom=False)
    # ax22.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax22.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax22.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax23 = plt.subplot(gs[23])
    iq = 2
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if ent_gaps[iq, i, j] > 1000:
                ent_gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax23.imshow(ent_gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_ent_gap, vmax=max_ent_gap)
    ax23.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax23.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\Delta_\\xi$', fontsize=axlabel_fontsize)
    ax23.grid()
    ax23.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax23.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax23.tick_params(labelleft=False)
    # ax23.tick_params(labelbottom=False)
    # ax23.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax23.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax23.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax24 = plt.subplot(gs[24])
    iq = 3
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if ent_gaps[iq, i, j] > 1000:
                ent_gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax24.imshow(ent_gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_ent_gap, vmax=max_ent_gap)
    ax24.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax24.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\Delta_\\xi$', fontsize=axlabel_fontsize)
    ax24.grid()
    ax24.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax24.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax24.tick_params(labelleft=False)
    # ax24.tick_params(labelbottom=False)
    # ax24.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax24.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax24.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax25 = plt.subplot(gs[25])
    iq = 4
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if ent_gaps[iq, i, j] > 1000:
                ent_gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax25.imshow(ent_gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_ent_gap, vmax=max_ent_gap)
    ax25.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax25.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\Delta_\\xi$', fontsize=axlabel_fontsize)
    ax25.grid()
    ax25.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax25.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax25.tick_params(labelleft=False)
    # ax25.tick_params(labelbottom=False)
    # ax25.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax25.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax25.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax26 = plt.subplot(gs[26])
    iq = 5
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if ent_gaps[iq, i, j] > 1000:
                ent_gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax26.imshow(ent_gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                     extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_ent_gap, vmax=max_ent_gap)
    ax26.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax26.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    # cbar.set_label('$\\Delta_\\xi$', fontsize=axlabel_fontsize)
    ax26.grid()
    ax26.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax26.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax26.tick_params(labelleft=False)
    # ax26.tick_params(labelbottom=False)
    # ax26.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax26.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax26.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax27 = plt.subplot(gs[27])
    iq = 6
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if ent_gaps[iq, i, j] > 1000:
                ent_gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax27.imshow(ent_gaps[iq].T, cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max], vmin=min_ent_gap, vmax=max_ent_gap)
    ax27.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax27.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    # cbar = plt.colorbar(sc)
    cb_ax = fig.add_axes([.9, .11, .00429, .176])
    cbar = fig.colorbar(sc, orientation='vertical', cax=cb_ax, format=FuncFormatter(fmt))
    cbar.set_label('$\\Delta_\\xi$', fontsize=axlabel_fontsize)
    ax27.grid()
    ax27.set_xlabel('$t_6$', fontsize=axlabel_fontsize)
    # ax27.set_ylabel('$t_9$', fontsize=axlabel_fontsize)
    ax27.tick_params(labelleft=False)
    # ax27.tick_params(labelbottom=False)
    # ax27.set_title(f'$n_\\phi = 1/{qs[iq]}$', fontsize=title_fontsize)
    ax27.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax27.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # if stats == "fermions":
    #     ax27.plot(-0.15, 0.05, 'x', c='r')

    label_fontsize = 13
    # fig.text(0.095, 0.923, "$\mathrm{(a)}$", fontsize=label_fontsize)
    # fig.text(0.095, 0.71, "$\mathrm{(b)}$", fontsize=label_fontsize)
    # fig.text(0.095, 0.5, "$\mathrm{(c)}$", fontsize=label_fontsize)
    # fig.text(0.095, 0.285, "$\mathrm{(d)}$", fontsize=label_fontsize)

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/{file_name}.png", bbox_inches='tight', dpi=300)
    plt.show()
