import numpy as np
from joblib import Parallel, delayed


def H(t_3, q, kx, ky):

    n_phi = p/q

    def A(m, ky):
        return -2 * (t_1 * np.cos(2 * np.pi * n_phi * m + ky) + t_3 * np.cos(4 * np.pi * n_phi * m + 2 * ky))

    def B(kx):
        return -t_1 * np.exp(1j * kx)

    def C(kx):
        return -t_3 * np.exp(1j * 2 * kx)

    B_arr_temp = np.diag(np.array([B(kx) for i in range(q)]))
    C_arr_temp = np.diag(np.array([C(kx) for i in range(q)]))

    A_arr = np.diag(np.array([A(i, ky) for i in range(q)]))
    B_arr = np.roll(B_arr_temp, 1, axis=1)
    Bc_arr = np.conjugate(np.roll(B_arr_temp, 1, axis=0))
    C_arr = np.roll(C_arr_temp, 2, axis=1)
    Cc_arr = np.conjugate(np.roll(C_arr_temp, 2, axis=0))

    return A_arr + B_arr + Bc_arr + C_arr + Cc_arr


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


def qgt(H_eigenvectors_arr, i, ikx, iky):
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


def fs_metric(H_eigenvectors_arr, i, ikx, iky):
    return np.real(qgt(H_eigenvectors_arr, i, ikx, iky))


def berry_curv(H_eigenvectors_arr, i, ikx, iky):
    return -2 * np.imag(qgt(H_eigenvectors_arr, i, ikx, iky))


def TISM(H_eigenvectors_arr, i, xi, yi):
    return np.trace(fs_metric(H_eigenvectors_arr, i, xi, yi)) - np.abs(berry_curv(H_eigenvectors_arr, i ,xi ,yi)[0,1])


def TISM_int(t_3, q, i):
    k_x, k_y = np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi, grain)
    dx, dy = (k_x[1] - k_x[0]), (k_y[1] - k_y[0])
    tdx, tdy = (k_x[1] - k_x[0])/1000, (k_y[1] - k_y[0])/1000
    H_eigenvectors_arr_0 = H_eigenvectors2(t_3, q, np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi/q, grain))
    H_eigenvectors_arr_1 = H_eigenvectors2(t_3, q, np.linspace(0, 2 * np.pi/q, grain)+tdx, np.linspace(0, 2 * np.pi/q, grain))
    H_eigenvectors_arr_2 = H_eigenvectors2(t_3, q, np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi/q, grain)+tdy)
    H_temp = np.array([H_eigenvectors_arr_0, H_eigenvectors_arr_1, H_eigenvectors_arr_2])
    n_phi = p/q
    tisms = np.empty((len(k_x)-1, len(k_y)-1))
    for xi in range(len(k_x)-1):
            for yi in range(len(k_y)-1):
                tisms[xi][yi] = (np.trace(fs_metric(H_temp, 0, xi, yi)) - np.abs(berry_curv(H_temp, 0 ,xi ,yi)[0,1]))
    return np.sum(tisms)*dx*dy/(2*np.pi)


if __name__ == "__main__":

    print('running')

    t_1 = 1
    grain = 100
    p = 1
    qs = np.array([96])

    def my_task(q, t_temp):
        n_phi = 1/q
        points = np.array([None, None])
        t_3 = t_temp
        k_x, k_y = np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi, grain)
        tism = np.array([t_3, TISM_int(t_3, q, 0)])
        return tism

    for q in qs:
        k_x, k_y = np.linspace(0, 2 * np.pi/q, grain), np.linspace(0, 2 * np.pi, grain)
        results = np.array(Parallel(n_jobs=3)(delayed(my_task)(q, t_temp) for t_temp in np.linspace(-0.25, 0, 26)))

    file = open(f"sp_data/q_96.txt", "w")
    for i, t3hop in enumerate(np.linspace(-0.25, 0, 26)):
        file.write(f"{results[:, 0][i]:.2f}\t{results[:, 1][i]}\n")
    file.close()
