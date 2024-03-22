import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def hexic_line(t6_val):
    return (1/64) * (1 - 15*t6_val)


if __name__ == "__main__":

    stats = "bosons"
    alpha = 0
    q = 144
    ts = np.linspace(-0.25, 0.25, 11)
    file_name = f"{stats}_alpha_{alpha}_q_{q}_ent"

    mb_data = f"/home/bart/DiagHam_latest/run/SFCI_data_2/{stats}_alpha_{alpha:g}/q_{q:g}/ent/mb_ent_q_{q:g}.txt"

    ts_min = np.min(ts)
    ts_max = np.max(ts)
    ts_len = len(ts)
    ts_step = (ts_max - ts_min) / (ts_len-1)

    plot_min = ts_min - ts_step/2
    plot_max = ts_max + ts_step/2
    plot_ts = np.linspace(plot_min, plot_max, 100)

    gaps = np.zeros((len(ts), len(ts)))

    with open(mb_data, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            it6 = int(i / len(ts))
            it9 = i % len(ts)
            gaps[it6, it9] = 2*float(row[2])

    fig = plt.figure()
    ax1 = plt.subplot(111)

    sc = ax1.imshow(gaps.T, cmap='magma', origin='lower',
                    extent=[plot_min, plot_max, plot_min, plot_max])
    ax1.plot(plot_ts, [hexic_line(i) for i in plot_ts], c='g', label="hexic line")
    ax1.scatter(1/7, -1/56, c='g', label="octic point")

    cbar = plt.colorbar(sc)
    cbar.set_label('$\\Delta_\\xi$')

    plt.grid()

    leg = ax1.legend(loc='best', ncol=3, handletextpad=0.5, handlelength=1, labelspacing=0,
                     borderpad=0.35, framealpha=1, markerscale=0.8, fontsize=10, columnspacing=0.5)

    ax1.set_xlabel('$t_6$')
    ax1.set_ylabel('$t_9$')
    ax1.set_title(f'quartic plane, $n_\\phi = 1/{q}$')

    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    
    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/{file_name}.png", bbox_inches='tight', dpi=300)
    plt.show()
