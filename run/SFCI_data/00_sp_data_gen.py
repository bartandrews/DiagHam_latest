import numpy as np
from joblib import Parallel, delayed
import time


def H(t_3, q, kx, ky):

    n_phi = 1/q

    def A(m, ky):
        return -2*(np.cos(2*np.pi*n_phi*m + ky) + t_3*np.cos(4*np.pi*n_phi*m + 2*ky))

    def B(kx):
        return -np.exp(1j*kx)

    def C(kx):
        return -t_3*np.exp(1j*2*kx)

    B_arr_temp = np.diag(np.array([B(kx) for _ in range(q)]))
    C_arr_temp = np.diag(np.array([C(kx) for _ in range(q)]))

    A_arr = np.diag(np.array([A(i, ky) for i in range(q)]))
    B_arr = np.roll(B_arr_temp, 1, axis=1)
    Bc_arr = np.conjugate(np.roll(B_arr_temp, 1, axis=0))
    C_arr = np.roll(C_arr_temp, 2, axis=1)
    Cc_arr = np.conjugate(np.roll(C_arr_temp, 2, axis=0))

    return A_arr + B_arr + Bc_arr + C_arr + Cc_arr


def H_eigbasis(t_3, q, k_x, k_y):
    evals_arr = np.empty((q, len(k_x), len(k_y)))
    evecs_arr = np.empty((q, len(k_x), len(k_y), q), dtype="complex_")
    for i, kx in enumerate(k_x):
        for j, ky in enumerate(k_y):
            evals, evecs = np.linalg.eigh(H(t_3, q, kx, ky))
            for a in range(q):
                evals_arr[a, i, j] = evals[a]
                evecs_arr[a, i, j, :] = evecs[:, a]
    return evals_arr, evecs_arr


def qgt(H_evecs_arr, i, ikx, iky):
    tdx, tdy = (k_x[1] - k_x[0]) / 1000, (k_y[1] - k_y[0]) / 1000
    eig_0 = H_evecs_arr[0][i, ikx % len(k_x), iky % len(k_y), :]
    eig_x = H_evecs_arr[1][i, ikx % len(k_x), iky % len(k_y), :]
    eig_y = H_evecs_arr[2][i, ikx % len(k_x), iky % len(k_y), :]
    eigs = {}
    eigs.update({"x": eig_x})
    eigs.update({"y": eig_y})
    grad = {}
    grad.update({"x": (eig_x - eig_0) / tdx})
    grad.update({"y": (eig_y - eig_0) / tdy})
    return np.array(
        [[np.vdot(grad[u], grad[v])
          - np.vdot(grad[u], eig_0) * np.vdot(eig_0, grad[v])
          for u in ["x", "y"]] for v in ["x", "y"]])


def fs_metric(H_evecs_arr, i, ikx, iky):
    return np.real(qgt(H_evecs_arr, i, ikx, iky))


# compute Berry curvature using qgt
def berry_curv(H_evecs_arr, i, ikx, iky):
    return -2*np.imag(qgt(H_evecs_arr, i, ikx, iky))


# compute Berry curvature using link variables (faster)
def berry_curv_link(_eigenvectors, _band, _idx_x, _idx_y):

    def _U(var_num, __eigenvectors, __band, __idx_x, __idx_y):

        vec1 = __eigenvectors[__band, __idx_x, __idx_y, :]
        if var_num == 1:
            vec2 = __eigenvectors[__band, __idx_x + 1, __idx_y, :]
        elif var_num == 2:
            vec2 = __eigenvectors[__band, __idx_x, __idx_y + 1, :]
        else:
            raise ValueError("link variable number must be in [1, 2].")

        return np.conj(vec1).dot(vec2)

    Berry_curv = - np.imag(np.log(_U(1, _eigenvectors, _band, _idx_x, _idx_y)
                                  * _U(2, _eigenvectors, _band, _idx_x+1, _idx_y)
                                  * _U(1, _eigenvectors, _band, _idx_x, _idx_y+1)**-1
                                  * _U(2, _eigenvectors, _band, _idx_x, _idx_y)**-1))

    return Berry_curv


def TISM(H_evecs_arr, i, xi, yi):
    return np.trace(fs_metric(H_evecs_arr, i, xi, yi)) - np.abs(berry_curv(H_evecs_arr, i, xi, yi)[0, 1])


def DISM(H_evecs_arr, i, xi, yi):
    return np.linalg.det(fs_metric(H_evecs_arr, i, xi, yi)) - (berry_curv(H_evecs_arr, i, xi, yi)[0, 1])**2/4


def band_geom(t_3, q):
    print(f"t3 = {t_3}")

    k_x, k_y = np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi, grain)
    dx, dy = (k_x[1] - k_x[0]), (k_y[1] - k_y[0])
    tdx, tdy = (k_x[1] - k_x[0])/1000, (k_y[1] - k_y[0])/1000

    H_evals_arr_0, H_evecs_arr_0 = H_eigbasis(t_3, q, k_x, k_y)
    _, H_evecs_arr_1 = H_eigbasis(t_3, q, k_x+tdx, k_y)
    _, H_evecs_arr_2 = H_eigbasis(t_3, q, k_x, k_y+tdy)
    H_temp = np.array([H_evecs_arr_0, H_evecs_arr_1, H_evecs_arr_2])

    tisms = np.empty((len(k_x)-1, len(k_y)-1))
    disms = np.empty((len(k_x)-1, len(k_y)-1))
    berry_fluxes = np.zeros((len(k_x)-1, len(k_y)-1))
    g_array = np.zeros((2, 2))

    for u in range(2):
        for v in range(2):
            double_g = np.zeros((len(k_x) - 1, len(k_y) - 1))
            single_g1 = np.zeros((len(k_x) - 1, len(k_y) - 1))
            single_g2 = np.zeros((len(k_x) - 1, len(k_y) - 1))
            for xi in range(len(k_x)-1):
                for yi in range(len(k_y)-1):
                    if u == 0 and v == 0:
                        tisms[xi][yi] = TISM(H_temp, 0, xi, yi)
                        disms[xi][yi] = DISM(H_temp, 0, xi, yi)
                        berry_fluxes[xi, yi] = berry_curv_link(H_evecs_arr_0, 0, xi, yi)
                    double_g[xi][yi] = fs_metric(H_temp, 0, xi, yi)[u, v] * fs_metric(H_temp, 0, xi, yi)[v, u]
                    single_g1[xi][yi] = fs_metric(H_temp, 0, xi, yi)[u, v]
                    single_g2[xi][yi] = fs_metric(H_temp, 0, xi, yi)[v, u]
            double_g_av = np.average(double_g)
            single_g1_av = np.average(single_g1)
            single_g2_av = np.average(single_g2)
            g_array[u, v] = double_g_av - single_g1_av * single_g2_av

    tism_int = np.sum(tisms) * dx * dy / (2 * np.pi)
    dism_int = np.sum(disms) * dx * dy / (2 * np.pi)
    fs_fluc = np.log10(np.sqrt(0.5 * np.sum(g_array)))

    # berry fluctuations
    A_BZ = np.sum(berry_fluxes) / np.average(berry_fluxes)
    berry_fluc = np.log10((A_BZ / (2 * np.pi)) * np.std(berry_fluxes))

    # gap-to-width
    gap = np.min(H_evals_arr_0[1]) - np.max(H_evals_arr_0[0])
    width = np.max(H_evals_arr_0[0]) - np.min(H_evals_arr_0[0])
    gap_width = np.log10(gap) - np.log10(width)

    return [t_3, tism_int, dism_int, berry_fluc, fs_fluc, gap_width]


if __name__ == "__main__":

    t0 = time.time()

    q = 54
    grain = 100
    ts = np.linspace(-0.25, 0, 26)

    k_x, k_y = np.linspace(0, 2*np.pi/q, grain), np.linspace(0, 2*np.pi, grain)
    results = np.array(Parallel(n_jobs=6)(delayed(band_geom)(t3, q) for t3 in ts))

    file = open(f"sp_data/q_{q}.txt", "w")
    for i, t3hop in enumerate(ts):
        file.write(f"{results[:, 0][i]:.2f}\t{results[:, 1][i]}\t{results[:, 2][i]}\t{results[:, 3][i]}\t"
                   f"{results[:, 4][i]}\t{results[:, 5][i]}\n")
    file.close()

    print("Total time taken (seconds) = ", time.time() - t0)
