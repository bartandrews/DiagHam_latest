import csv
import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

if __name__ == "__main__":

    program = "tune"  # "tune" or "turn"
    stats = "fermions"  # "fermions" or "bosons"

    if program == "tune":
        if stats == "fermions":
            scale_exp = 2
            N_list = [6]
            p_list = [9, 19, 33, 51, 73, 99, 129, 163]
            outliers = [(6, 163, 0.4)]  # [(N, p, alpha)]
        else:
            scale_exp = 1
            N_list = [6]
            p_list = [13, 28, 49, 76, 109, 148, 193]
            outliers = []  # [(N, p, alpha)]

        alpha = np.empty((len(N_list), len(p_list)), dtype=object)
        gap = np.empty((len(N_list), len(p_list)), dtype=object)

        for N_idx, _ in enumerate(N_list):
            for p_idx, _ in enumerate(p_list):
                alpha[N_idx, p_idx] = []
                gap[N_idx, p_idx] = []

        with open(f"tune/{stats}/{stats}_alpha_gap.dat", 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for row in plots:
                for N_idx, N in enumerate(N_list):
                    for p_idx, p in enumerate(p_list):
                        if int(row[0]) == N and int(row[1]) == p:
                            alpha_val = float(row[2])
                            gap_val = float(row[3])
                            if (N, p, alpha_val) in outliers:
                                continue
                            else:
                                alpha[N_idx, p_idx].append(alpha_val)
                                gap[N_idx, p_idx].append((p-1)**scale_exp*gap_val)

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
                    "$V_{ij} = (1-\\alpha) (e\\delta_{ij} + \\delta_{\\langle ij \\rangle}) " \
                    "+ \\alpha e^{1-|\\mathbf{r}_{ij}|^4}$"

        ax1.set_xlabel('$\\alpha$')
        ax1.set_ylabel(ylabel)
        ax1.set_title(title)
        plt.savefig(f"/home/bart/DiagHam_latest/run/hof_int/tune/{stats}/{stats}_alpha_gap.png",
                    bbox_inches='tight', dpi=300)
    else:
        if stats == "fermions":
            scale_exp = 2
            N_list = [6]
            p_list = [9]
            outliers = [(6, 163, 0.4)]  # [(N, p, theta)]
        else:
            scale_exp = 1
            N_list = [6]
            p_list = [13, 28, 49, 76, 109, 148, 193]
            outliers = []  # [(N, p, theta)]

        theta = np.empty((len(N_list), len(p_list)), dtype=object)
        gap = np.empty((len(N_list), len(p_list)), dtype=object)

        for N_idx, _ in enumerate(N_list):
            for p_idx, _ in enumerate(p_list):
                theta[N_idx, p_idx] = []
                gap[N_idx, p_idx] = []

        with open(f"turn/{stats}/{stats}_theta_gap.dat", 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for row in plots:
                for N_idx, N in enumerate(N_list):
                    for p_idx, p in enumerate(p_list):
                        if int(row[0]) == N and int(row[1]) == p:
                            theta_val = float(row[2])
                            gap_val = float(row[3])
                            if (N, p, theta_val) in outliers:
                                continue
                            else:
                                theta[N_idx, p_idx].append(theta_val)
                                gap[N_idx, p_idx].append((p - 1) ** scale_exp * gap_val)

        fig = plt.figure()
        ax1 = plt.subplot(111)

        for N_idx, N in enumerate(N_list):
            for p_idx, p in enumerate(p_list):
                ax1.plot(theta[N_idx, p_idx], gap[N_idx, p_idx], '.-', label=f"${p}$")

        leg = ax1.legend(loc='best', title='$p$', ncol=3, handletextpad=0.5, handlelength=0, labelspacing=0,
                         borderpad=0.35, framealpha=1, markerscale=0.8, fontsize=10, columnspacing=0.5)

        if stats == "fermions":
            ylabel = "$q^2\\Delta_{\mathrm{m.b.}}$"
            title = "$N=6$ fermions with " \
                    "$V_{ij} = \\sum_{ij} e^{1-|\\mathbf{r}_{ij}|^2} \[ \\delta(y_j -\\tan(\\theta) x_i) " \
                    "+ \\delta(y_j + \\tan^{-1}(\\theta) x_i) \]$"
        else:
            ylabel = "$q\\Delta_{\mathrm{m.b.}}$"
            title = "$N=6$ bosons with " \
                    "$V_{ij} = \\sum_{ij} e^{1-|\\mathbf{r}_{ij}|^2} \[ \\delta(y_j -\\tan(\\theta) x_i) " \
                    "+ \\delta(y_j + \\tan^{-1}(\\theta) x_i) \]$"

        ax1.set_xlabel('$\\theta$')
        ax1.set_ylabel(ylabel)
        ax1.set_title(title)
        plt.savefig(f"/home/bart/DiagHam_latest/run/hof_int/turn/{stats}/{stats}_theta_gap.png",
                    bbox_inches='tight', dpi=300)
    plt.show()
