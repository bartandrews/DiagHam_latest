import numpy as np
from joblib import Parallel, delayed
# from mpmath import mpf


def H(t_3, q_val, k):

    kx, ky = k[0], k[1]
    n_phi = 1/q_val

    def A(m, ky_val):
        return -2*(1*np.cos(2*np.pi*n_phi*m + ky_val) + t_3*np.cos(4*np.pi*n_phi*m + 2*ky_val))

    def B(kx_val):
        return -1*np.exp(1j*kx_val)

    def C(kx_val):
        return -t_3*np.exp(1j*2*kx_val)

    B_arr_temp = np.diag(np.array([B(kx) for _ in range(q_val)]))
    C_arr_temp = np.diag(np.array([C(kx) for _ in range(q_val)]))

    A_arr = np.diag(np.array([A(i, ky) for i in range(q_val)]))
    B_arr = np.roll(B_arr_temp, 1, axis=1)
    Bc_arr = np.conj(np.roll(B_arr_temp, 1, axis=0))
    C_arr = np.roll(C_arr_temp, 2, axis=1)
    Cc_arr = np.conj(np.roll(C_arr_temp, 2, axis=0))

    return A_arr + B_arr + Bc_arr + C_arr + Cc_arr


def berry_curv(_eigenvectors, _band, _idx_x, _idx_y):

    def _U(var_num, __eigenvectors, __band, __idx_x, __idx_y):

        vec1 = __eigenvectors[:, __band, __idx_x, __idx_y]
        if var_num == 1:
            vec2 = __eigenvectors[:, __band, __idx_x + 1, __idx_y]
        elif var_num == 2:
            vec2 = __eigenvectors[:, __band, __idx_x, __idx_y + 1]
        else:
            raise ValueError("link variable number must be in [1, 2].")

        return np.conj(vec1).dot(vec2)

    Berry_curv = - np.imag(np.log(_U(1, _eigenvectors, _band, _idx_x, _idx_y)
                                  * _U(2, _eigenvectors, _band, _idx_x+1, _idx_y)
                                  * _U(1, _eigenvectors, _band, _idx_x, _idx_y+1)**-1
                                  * _U(2, _eigenvectors, _band, _idx_x, _idx_y)**-1))

    return Berry_curv



def H_eigenvectors2(t_3, q, k_x, k_y):
    eigenvectors_arr = np.empty((q, len(k_x), len(k_y), q), dtype="complex_")
    for i, kx in enumerate(k_x):
        for j, ky in enumerate(k_y):
            matrix = H(t_3, q, kx, ky)
            assert(np.all(0 == (matrix - np.conj(matrix.T))))  # Hermitian
            eigenvalues, eigenvectors = np.linalg.eigh(matrix)
            for a in range(q):
                eigenvectors_arr[a, i, j, :] = eigenvectors[:, a]
    return eigenvectors_arr


def qgt(k_x, k_y, H_eigenvectors_arr, i, ikx, iky):
    dx, dy = (k_x[1] - k_x[0]) / 1000, (k_y[1] - k_y[0]) / 1000
    eig_0 = H_eigenvectors_arr[0][i, ikx % len(k_x), iky % len(k_y), :]
    eig_x = H_eigenvectors_arr[1][i, ikx % len(k_x), iky % len(k_y), :]
    eig_y = H_eigenvectors_arr[2][i, ikx % len(k_x), iky % len(k_y), :]
    eigs = {}
    eigs.update({"x": eig_x})
    eigs.update({"y": eig_y})
    grad = {}
    grad.update({"x": (eig_x - eig_0) / dx})
    grad.update({"y": (eig_y - eig_0) / dy})
    return np.array(
        [[np.vdot(grad[u], grad[v]) - np.vdot(grad[u], eig_0) * np.vdot(eig_0, grad[v]) for u in ["x", "y"]] for v in
         ["x", "y"]])


def fs_metric(k_x, k_y, H_eigenvectors_arr, i, ikx, iky):
    return np.real(qgt(k_x, k_y, H_eigenvectors_arr, i, ikx, iky))


def berry_curv(k_x, k_y, H_eigenvectors_arr, i, ikx, iky):
    return -2 * np.imag(qgt(k_x, k_y, H_eigenvectors_arr, i, ikx, iky))


def TISM(k_x, k_y, H_eigenvectors_arr, i, xi, yi):
    return (np.trace(fs_metric(k_x, k_y, H_eigenvectors_arr, i, xi, yi))
            - np.abs(berry_curv(k_x, k_y, H_eigenvectors_arr, i ,xi ,yi)[0,1]))


def TISM_int(H_eigenvectors_arr, t_3, q, grain, i):
    k_x, k_y = np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi, grain)
    dx, dy = (k_x[1] - k_x[0]), (k_y[1] - k_y[0])
    tdx, tdy = (k_x[1] - k_x[0])/1000, (k_y[1] - k_y[0])/1000
    H_eigenvectors_arr_0 = H_eigenvectors_arr(t_3, q, np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi/q, grain))
    H_eigenvectors_arr_1 = H_eigenvectors_arr(t_3, q, np.linspace(0, 2 * np.pi/q, grain)+tdx, np.linspace(0, 2 * np.pi/q, grain))
    H_eigenvectors_arr_2 = H_eigenvectors_arr(t_3, q, np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi/q, grain)+tdy)
    H_temp = np.array([H_eigenvectors_arr_0, H_eigenvectors_arr_1, H_eigenvectors_arr_2])
    n_phi = 1/q
    tisms = np.empty((len(k_x)-1, len(k_y)-1))


    for xi in range(len(k_x)-1):
            for yi in range(len(k_y)-1):
                tisms[xi][yi] = TISM(k_x, k_y, H_temp, 0, xi, yi)

    return np.sum(tisms)*dx*dy/(2*np.pi)






def my_task(q_val, t3):
    print(f"t3 = {t3}")

    # construct bands
    eigenvalues = np.zeros((q_val, num_samples, num_samples))  # real
    eigenvectors = np.zeros((q_val, q_val, num_samples, num_samples), dtype=np.complex128)  # complex
    for band in range(q_val):
        for idx_x in range(num_samples):
            frac_kx = idx_x / (num_samples - 1)
            for idx_y in range(num_samples):
                frac_ky = idx_y / (num_samples - 1)
                k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                eigvals, eigvecs = np.linalg.eig(H(t3, q_val, k))
                idx = np.argsort(eigvals)
                eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    # compute Berry fluxes
    berry_fluxes = np.zeros((q_val, num_samples-1, num_samples-1))  # real
    for band in range(q_val):
        for idx_x in range(num_samples - 1):
            for idx_y in range(num_samples - 1):
                berry_fluxes[band, idx_x, idx_y] = berry_curv(eigenvectors, band, idx_x, idx_y)

    # gap-to-width
    gap = np.min(eigenvalues[1])-np.max(eigenvalues[0])
    width = np.max(eigenvalues[0])-np.min(eigenvalues[0])
    gap_width = np.log(gap) - np.log(width)

    # berry fluctuations
    A_BZ = np.sum(berry_fluxes[0])/np.average(berry_fluxes[0])
    berry_fluc = np.log((A_BZ/(2*np.pi)) * np.std(berry_fluxes[0]))

    # TISM
    tism = TISM_int(t3, q, 0)

    return [t3, gap_width, berry_fluc]


if __name__ == "__main__":

    num_samples = 100
    q = 16

    b1 = (2.*np.pi) / 1 * np.array([1/q, 0])
    b2 = (2.*np.pi) / 1 * np.array([0, 1])
    bvec = np.vstack((b1, b2))

    file = open(f"sp_data/q_{q}_aux.txt", "w")

    ts = list(np.linspace(-0.25, 0, 26))
    results = np.array(Parallel(n_jobs=6)(delayed(my_task)(q, t_temp) for t_temp in ts))
    results = results[results[:, 0].argsort()]  # sort by first column

    for i in range(np.shape(results)[0]):
        file.write(f"{results[i, 0]:.2f}\t{results[i, 1]}\t{results[i, 2]}\n")
    file.close()
