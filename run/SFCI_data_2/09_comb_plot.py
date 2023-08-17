import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def hexic_line(t6_val):
    return (1/64) * (1 - 15*t6_val)


if __name__ == "__main__":

    stats = "bosons"
    alpha = 0
    qs = [16, 49, 81]
    ts = np.linspace(-0.25, 0.25, 11)
    file_name = f"{stats}_alpha_{alpha}"

    if stats == "bosons":
        cbtitle = "$q \\Delta_\\mathrm{m.b.}$"
        title_str = f"Bosons, $\\nu=1/2$, $V_{{ij}} = (1-\\alpha)\\delta_{{ij}} + \\alpha e^{{-|r_{{ij}}|^4}}$, $\\alpha={alpha}$"
    else:  # fermions
        cbtitle = "$q^2 \\Delta_\\mathrm{m.b.}$"
        title_str = f"Fermions, $\\nu=1/3$, $V_{{ij}} = (1-\\alpha)\\delta_{{\\langle ij \\rangle}} + \\alpha e^{{-|r_{{ij}}|^4}}$, $\\alpha={alpha}$"

    sp_data, mb_data, mb_ent_data = [], [], []
    for i, _ in enumerate(qs):
        sp_data.append(f"/home/bart/DiagHam_latest/run/SFCI_data_2/sp_data/q_{qs[i]:g}.txt")
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
    ent_gaps = np.zeros((len(qs), len(ts), len(ts)))

    # single-particle quantities
    for iq, _ in enumerate(qs):
        with open(sp_data[iq], 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for i, row in enumerate(plots):
                it6 = int(i / len(ts))
                it9 = i % len(ts)
                tisms[iq, it6, it9] = float(row[2])

    # many-body quantities
    for iq in [0]:
        with open(mb_data[iq], 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for i, row in enumerate(plots):
                it6 = int(i / len(ts))
                it9 = i % len(ts)
                gaps[iq, it6, it9] = float(row[2])
    for iq in [0]:
        with open(mb_ent_data[iq], 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for i, row in enumerate(plots):
                it6 = int(i / len(ts))
                it9 = i % len(ts)
                ent_gaps[iq, it6, it9] = float(row[2])

    fig = plt.figure(figsize=(9, 9))
    fig.suptitle(f'{title_str}')
    gs = gridspec.GridSpec(3, 3)
    gs.update(top=0.9)

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
    sc = ax0.imshow(tisms[iq], cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max])
    ax0.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax0.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    cbar = plt.colorbar(sc)
    # cbar.set_label('$\\langle \\mathcal{T} \\rangle$')
    ax0.grid()
    # ax0.set_xlabel('$t_6$')
    ax0.set_ylabel('$t_9$')
    # ax0.tick_params(labelleft=False)
    ax0.tick_params(labelbottom=False)
    ax0.set_title(f'$n_\\phi = 1/{qs[iq]}$')
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
    sc = ax1.imshow(tisms[iq], cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max])
    ax1.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax1.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    cbar = plt.colorbar(sc)
    # cbar.set_label('$\\langle \\mathcal{T} \\rangle$')
    ax1.grid()
    # ax1.set_xlabel('$t_6$')
    # ax1.set_ylabel('$t_9$')
    ax1.tick_params(labelleft=False)
    ax1.tick_params(labelbottom=False)
    ax1.set_title(f'$n_\\phi = 1/{qs[iq]}$')
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
    sc = ax2.imshow(tisms[iq], cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max])
    ax2.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax2.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    cbar = plt.colorbar(sc)
    cbar.set_label('$\\langle \\mathcal{T} \\rangle$')
    ax2.grid()
    # ax2.set_xlabel('$t_6$')
    # ax2.set_ylabel('$t_9$')
    ax2.tick_params(labelleft=False)
    ax2.tick_params(labelbottom=False)
    ax2.set_title(f'$n_\\phi = 1/{qs[iq]}$')
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    #################
    # Many-body gap #
    #################

    ax3 = plt.subplot(gs[3])
    iq = 0
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if gaps[iq, i, j] > 1000:
                gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax3.imshow(gaps[iq], cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max])
    ax3.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax3.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    cbar = plt.colorbar(sc)
    cbar.set_label(f'{cbtitle}')
    ax3.grid()
    # ax3.set_xlabel('$t_6$')
    ax3.set_ylabel('$t_9$')
    # ax3.tick_params(labelleft=False)
    ax3.tick_params(labelbottom=False)
    # ax3.set_title(f'$n_\\phi = 1/{qs[iq]}$')
    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ####################
    # Entanglement gap #
    ####################

    ax6 = plt.subplot(gs[6])
    iq = 0
    # remove outliers
    for i in range(ts_len):
        for j in range(ts_len):
            if ent_gaps[iq, i, j] > 1000:
                ent_gaps[iq, i, j] = np.nan
    # smooth plot
    itrpol_method = 'None'  # 'bicubic'
    # plot
    sc = ax6.imshow(ent_gaps[iq], cmap='magma', origin='lower', interpolation=itrpol_method,
                    extent=[plot_min, plot_max, plot_min, plot_max])
    ax6.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax6.scatter(1 / 7, -1 / 56, c='g', label="octic point")
    cbar = plt.colorbar(sc)
    cbar.set_label('$\\Delta_\\xi$')
    ax6.grid()
    ax6.set_xlabel('$t_6$')
    ax6.set_ylabel('$t_9$')
    # ax6.tick_params(labelleft=False)
    # ax6.tick_params(labelbottom=False)
    # ax6.set_title(f'$n_\\phi = 1/{qs[iq]}$')
    ax6.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax6.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/{file_name}.png", bbox_inches='tight', dpi=300)
    plt.show()
