import csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

if __name__ == "__main__":

    stats = "bosons"  # "fermions" or "bosons"

    if stats == "fermions":
        scale_exp = 2
        N_list = [6]
        p_list = [9, 19, 33, 51, 73, 99, 129, 163]
    else:
        scale_exp = 1
        N_list = [6]
        p_list = [13, 28, 49, 76, 109, 148, 193]

    alpha = np.empty((len(N_list), len(p_list)), dtype=object)
    gap = np.empty((len(N_list), len(p_list)), dtype=object)

    for N_idx, _ in enumerate(N_list):
        for p_idx, _ in enumerate(p_list):
            alpha[N_idx, p_idx] = []
            gap[N_idx, p_idx] = []

    with open(f"{stats}/{stats}_alpha_gap.dat", 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            for N_idx, N in enumerate(N_list):
                for p_idx, p in enumerate(p_list):
                    if int(row[0]) == N and int(row[1]) == p:
                        alpha[N_idx, p_idx].append(float(row[2]))
                        gap[N_idx, p_idx].append((p-1)**scale_exp*float(row[3]))

    fig = plt.figure()
    ax1 = plt.subplot(111)

    for N_idx, N in enumerate(N_list):
        for p_idx, p in enumerate(p_list):
            ax1.plot(alpha[N_idx, p_idx], gap[N_idx, p_idx], '.-', label=f"${p}$")

    leg = ax1.legend(loc='best', title='$p$', ncol=3, handletextpad=0.5, handlelength=0, labelspacing=0,
                     borderpad=0.35, framealpha=1, markerscale=0.8, fontsize=10, columnspacing=0.5)

    if stats == "fermions":
        ylabel = "$q^2\\Delta_{\mathrm{m.b.}}$"
        title = "$N=6$ fermions with interaction strength " \
                "$V_{ij} = (1-\\alpha) \\delta_{\\langle ij \\rangle} + \\alpha e^{1-|\\mathbf{r}_{ij}|^4}$"
    else:
        ylabel = "$q\\Delta_{\mathrm{m.b.}}$"
        title = "$N=6$ bosons with interaction strength " \
                "$V_{ij} = (1-\\alpha) \\delta_{ij} + \\alpha e^{-|\\mathbf{r}_{ij}|^4}$"

    ax1.set_xlabel('$\\alpha$')
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    plt.savefig(f"/home/bart/DiagHam_latest/run/hof_int/tune/{stats}/{stats}_alpha_gap.png",
                bbox_inches='tight', dpi=300)
    plt.show()
