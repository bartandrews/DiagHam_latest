import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import csv


def can_convert_to_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def H(p, q, k, t3=0.0, t6=0.0, t9=0.0):

    kx = k[0]
    ky = k[1]

    nphi = p/q
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    def A(m, ky):
        A_term = (- 2*np.cos(2*np.pi*nphi*m + ky) - 2*t3*np.cos(4*np.pi*nphi*m + 2*ky)
                  - 2*t6*np.cos(6*np.pi*nphi*m + 3*ky) - 2*t9*np.cos(8*np.pi*nphi*m + 4*ky))
        return A_term

    def B(kx):
        return -np.exp(1j * kx)

    def C(kx):
        return -t3*np.exp(1j*2*kx)

    def D(kx):
        return -t6*np.exp(1j*3*kx)

    def F(kx):
        return -t9*np.exp(1j*4*kx)

    # A terms
    for i in range(q):
        Hamiltonian[i][i] = A(i, ky)

    # B terms
    for i in range(q-1):
        Hamiltonian[i][i+1] = B(kx)
        Hamiltonian[i+1][i] = np.conj(B(kx))

    # C terms
    if t3 != 0 and q >= 5:
        for i in range(q-2):
            Hamiltonian[i][i+2] = C(kx)
            Hamiltonian[i + 2][i] = np.conj(C(kx))
    elif t3 != 0:
        raise ValueError("q is too small")

    # D terms
    if t6 != 0 and q >= 7:
        for i in range(q-3):
            Hamiltonian[i][i+3] = D(kx)
            Hamiltonian[i+3][i] = np.conj(D(kx))
    elif t6 != 0:
        raise ValueError("q is too small")

    # F terms
    if t9 != 0 and q >= 9:
        for i in range(q-4):
            Hamiltonian[i][i+4] = F(kx)
            Hamiltonian[i+4][i] = np.conj(F(kx))
    elif t9 != 0:
        raise ValueError("q is too small")

    # B boundary terms
    Hamiltonian[0][q-1] = np.conj(B(kx))
    Hamiltonian[q-1][0] = B(kx)

    # C boundary terms
    if t3 != 0 and q >= 5:
        Hamiltonian[0][q-2] = np.conj(C(kx))
        Hamiltonian[1][q-1] = np.conj(C(kx))
        Hamiltonian[q-2][0] = C(kx)
        Hamiltonian[q-1][1] = C(kx)
    elif t3 != 0:
        raise ValueError("q is too small")

    # D boundary terms
    if t6 != 0 and q >= 7:
        Hamiltonian[0][q-3] = np.conj(D(kx))
        Hamiltonian[1][q-2] = np.conj(D(kx))
        Hamiltonian[2][q-1] = np.conj(D(kx))
        Hamiltonian[q-3][0] = D(kx)
        Hamiltonian[q-2][1] = D(kx)
        Hamiltonian[q-1][2] = D(kx)
    elif t6 != 0:
        raise ValueError("q is too small")

    # F boundary terms
    if t9 != 0 and q >= 9:
        Hamiltonian[0][q-4] = np.conj(F(kx))
        Hamiltonian[1][q-3] = np.conj(F(kx))
        Hamiltonian[2][q-2] = np.conj(F(kx))
        Hamiltonian[3][q-1] = np.conj(F(kx))
        Hamiltonian[q-4][0] = F(kx)
        Hamiltonian[q-3][1] = F(kx)
        Hamiltonian[q-2][2] = F(kx)
        Hamiltonian[q-1][3] = F(kx)
    elif t9 != 0:
        raise ValueError("q is too small")

    return Hamiltonian


if __name__ == "__main__":


    GA = np.array([0, 0])
    Y = np.array([0, 0.5])
    S = np.array([0.5, 0.5])
    X = np.array([0.5, 0])
    sym_points = [GA, Y, S, X]
    num_samples = 100
    num_bands = 9
    b1 = (2. * np.pi) * np.array([1 / num_bands, 0])
    b2 = (2. * np.pi) * np.array([0, 1])
    bvec = np.vstack((b1, b2))

    # construct bands
    num_paths = len(sym_points)
    points_per_path = int(num_samples / num_paths)
    num_points = num_paths * points_per_path
    eigenvalues = np.zeros((num_bands, num_points))  # real
    count = 0
    for i in range(num_paths):
        for j in range(points_per_path):
            k = sym_points[i] + (sym_points[(i + 1) % num_paths] - sym_points[i]) * float(j) / float(
                points_per_path - 1)
            k = np.matmul(k, bvec)
            eigvals = np.linalg.eigvals(H(1, num_bands, k, -0.1125, -0.15, 0.05))
            idx = np.argsort(eigvals)
            for band in range(num_bands):
                eigenvalues[band, count] = np.real(eigvals[idx[band]])
            count += 1

    # construct butterfly
    q_val = 99
    nphi_list, E_list = [], []
    for p in range(1, q_val):
        nphi_val = p/q_val
        nphi_list.append([nphi_val]*q_val)
        E_list.append(np.linalg.eigvals(H(p, q_val, np.array([0, 0]), -0.1125, -0.15, 0.05)))

    # construct many-body spectrum
    mb_file = "fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_-0.11_t6_-0.15_t9_0.05_u_1_gx_0_gy_0_ext.dat"
    lin_K = []
    mb_E = []
    with open(mb_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K.append(float(row[0])*6+float(row[1]))
                mb_E.append(float(row[2])/2)

    # construct entanglement spectrum
    ent_file = "fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_-0.11_t6_-0.15_t9_0.05_u_1_gx_0_gy_0.na_4.parentspec"
    lin_K_ent = []
    ent = []
    with open(ent_file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i, row in enumerate(plots):
            if can_convert_to_float(row[0]):
                lin_K_ent.append(float(row[3]))
                ent.append(float(row[5]))

    # construct spectral flow
    gammas = []
    gamma_E_0, gamma_E_2, gamma_E_4 = [], [], []
    for gamma in np.linspace(0, 1, 21):
        gammas_E_0, gammas_E_2, gammas_E_4 = [], [], []
        gammas.append(gamma)
        flow_file = f"fermions_hofstadter_X_6_Y_4_q_1_n_8_x_4_y_6_t3_-0.11_t6_-0.15_t9_0.05_u_1_gx_{gamma:g}_gy_0.dat"
        with open(flow_file, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=' ')
            for i, row in enumerate(plots):
                if can_convert_to_float(row[0]):
                    if float(row[0]) == 0 and float(row[1]) == 0:
                        gammas_E_0.append(float(row[2])/2)
                    elif float(row[0]) == 0 and float(row[1]) == 2:
                        gammas_E_2.append(float(row[2]) / 2)
                    elif float(row[0]) == 0 and float(row[1]) == 4:
                        gammas_E_4.append(float(row[2]) / 2)
        gamma_E_0.append(np.min(gammas_E_0))
        gamma_E_2.append(np.min(gammas_E_0))
        gamma_E_4.append(np.min(gammas_E_0))

    # construct entanglement flow
    Nas = []
    ents = []
    for realNa in [1, 2, 3, 4, 5, 6, 7]:
        if realNa > 4:
            Na = 8 - realNa
        else:
            Na = realNa
        ent_flow_file = f"fermions_hofstadter_X_12_Y_8_q_1_n_8_x_4_y_6_t3_-0.11_t6_-0.15_t9_0.05_u_1_gx_0_gy_0.na_{Na}.parentspec"
        with open(ent_flow_file, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=' ')
            for i, row in enumerate(plots):
                if can_convert_to_float(row[0]):
                    if float(row[3]) == 0:
                        if realNa <= 4:
                            Nas.append(int(row[0]))
                        else:
                            Nas.append(realNa)
                        ents.append(float(row[5]))

    ##########
    # figure #
    ##########

    fig = plt.figure(figsize=(6, 9))
    gs = gridspec.GridSpec(3, 2, hspace=0.5, wspace=0.5)

    ############################
    # single-particle spectrum #
    ############################

    ax0 = plt.subplot(gs[0])
    ax0.set_title(f"$n_\phi = 1/{num_bands}$")
    for i in range(num_bands):
        if i == 5:
            color = 'b'
        else:
            color = 'r'
        ax0.plot(eigenvalues[i], c=color)
    xtick_vals = [0]
    for i in range(1, num_paths):
        ax0.axvline(i * points_per_path, color='k', linewidth=0.05, ls='--')
        xtick_vals.append(i * points_per_path)
    xtick_vals.append(num_paths*points_per_path)

    ax0.set_xticks(xtick_vals)
    ax0.set_xticklabels(["$\Gamma$", "$Y$", "$S$", "$X$", "$\Gamma$"])

    ax0.set_xlim([0, num_points])
    ax0.set_ylabel('$E$')
    # ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    Chern_list = [1, 1, 1, 1, 1, -8, 1, 1, 1]
    for i in range(num_bands):
        if i < 5:
            ax0.text(1, eigenvalues[i][10], f"$C_{i+1}={Chern_list[i]}$")
        elif i == 5:
            ax0.text(37, 1.5*eigenvalues[i][50], f"$C_{i + 1}={Chern_list[i]}$")
        else:
            ax0.text(77, 0.95*eigenvalues[i][90], f"$C_{i + 1}={Chern_list[i]}$")

    #############
    # butterfly #
    #############

    ax1 = plt.subplot(gs[1])
    ax1.set_title(f"$n_\phi = p/{q_val}$")
    ax1.plot(nphi_list, E_list, '.', c='k', markersize=1)
    ax1.set_ylabel('$E$')
    ax1.set_xlabel('$n_\phi$')
    ax1.set_xlim([0, 1])
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ######################
    # many-body spectrum #
    ######################

    ax2 = plt.subplot(gs[2])
    ax2.set_title(f"$n_\phi = 1/96$")
    ax2.plot(lin_K, np.subtract(mb_E, min(mb_E))*1e5, '+', c='k', markersize=5)
    ax2.set_ylabel('$(E_\\mathrm{m.b.} - E_{\\mathrm{m.b.},0})/10^{-5}$')
    ax2.set_xlabel('$k_x L_y + k_y$')
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.annotate(s='', xy=(3, 0.000128*1e5), xytext=(3, -0.000003*1e5), arrowprops=dict(arrowstyle='<->'))
    ax2.text(3.5, 0.00006*1e5, "$\\Delta_\\mathrm{m.b.}$")

    #########################
    # entanglement spectrum #
    #########################

    ax3 = plt.subplot(gs[3])
    ax3.set_title(f"$n_\phi = 1/96$")
    ax3.plot(lin_K_ent, ent, '+', c='k', markersize=2)
    ax3.set_ylabel('$\\xi$')
    ax3.set_xlabel('$k_x L_y + k_y$')
    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.annotate(s='', xy=(10, 17.5), xytext=(10, 9.5), arrowprops=dict(arrowstyle='<->'))
    ax3.text(10.5, 13, "$\\Delta_\\xi$")

    #################
    # spectral flow #
    #################

    offset = np.min(np.array([gamma_E_0, gamma_E_2, gamma_E_4]))

    ax4 = plt.subplot(gs[4])
    ax4.set_title(f"$n_\phi = 1/24$")
    ax4.plot(gammas, np.subtract(gamma_E_0, offset)*1e10, 's', c='b', markersize=5, label="$0$")
    ax4.plot(gammas, np.subtract(gamma_E_2, offset)*1e10, 'o', c='g', markersize=5, label="$2$")
    ax4.plot(gammas, np.subtract(gamma_E_4, offset)*1e10, '+', c='r', markersize=5, label="$4$")
    ax4.set_ylabel('$(E_\\mathrm{m.b.}-\\min_{\Phi}(E_{\\mathrm{m.b.},0}))/10^{-10}$')
    ax4.set_xlabel('$\\Phi/2\pi$')
    ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    ax4.legend(loc='best', handletextpad=0.3, handlelength=0.5, labelspacing=0.1, borderpad=0.3,
               framealpha=1, edgecolor='k', markerscale=0.8, fontsize=10, ncol=6, columnspacing=0.5, title='$k_y$')

    #####################
    # entanglement flow #
    #####################

    ax5 = plt.subplot(gs[5])
    ax5.set_title(f"$n_\phi = 1/96$")
    ax5.plot(Nas, ents, '+', c='k', markersize=5)
    ax5.set_ylabel('$\\xi$')
    ax5.set_xlabel('$N_A$')
    ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    fig.text(0.05, 0.9, "(a)", fontsize=12)
    fig.text(0.5, 0.9, "(b)", fontsize=12)
    fig.text(0.05, 0.61, "(c)", fontsize=12)
    fig.text(0.5, 0.61, "(d)", fontsize=12)
    fig.text(0.05, 0.32, "(e)", fontsize=12)
    fig.text(0.5, 0.32, "(f)", fontsize=12)

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/spmb/spmb.png", bbox_inches='tight', dpi=300)
    plt.show()
