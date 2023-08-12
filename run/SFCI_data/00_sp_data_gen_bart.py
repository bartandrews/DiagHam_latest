import numpy as np
from joblib import Parallel, delayed


def berry_curv(_eigenvectors, _band, _idx_x, _idx_y, _group_size=1):

    def _U(var_num, __eigenvectors, __band, __idx_x, __idx_y, __group_size):

        link_matrix = np.zeros((__group_size, __group_size), dtype=complex)
        for i in range(__group_size):
            for j in range(__group_size):
                vec1 = __eigenvectors[:, __band+i, __idx_x, __idx_y]
                if var_num == 1:
                    vec2 = __eigenvectors[:, __band+j, __idx_x + 1, __idx_y]
                elif var_num == 2:
                    vec2 = __eigenvectors[:, __band + j, __idx_x, __idx_y + 1]
                else:
                    raise ValueError("link variable number must be in [1, 2].")
                link_matrix[i, j] = np.conj(vec1).dot(vec2)
        link_var = np.linalg.det(link_matrix)
        return link_var

    Berry_curv = - np.imag(np.log(_U(1, _eigenvectors, _band, _idx_x, _idx_y, _group_size)
                                  * _U(2, _eigenvectors, _band, _idx_x+1, _idx_y, _group_size)
                                  * _U(1, _eigenvectors, _band, _idx_x, _idx_y+1, _group_size)**-1
                                  * _U(2, _eigenvectors, _band, _idx_x, _idx_y, _group_size)**-1))

    return Berry_curv


def H(t_3, q, k):

    kx = k[0]
    ky = k[1]

    n_phi = 1/q

    def A(m, ky):
        return -2 * (1 * np.cos(2 * np.pi * n_phi * m + ky) + t_3 * np.cos(4 * np.pi * n_phi * m + 2 * ky))

    def B(kx):
        return -1 * np.exp(1j * kx)

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


def my_task(q, t3):

    # i = ts.index(t3)
    print(f"t3 = {t3}")

    # construct bands
    eigenvalues = np.zeros((q, num_samples, num_samples))  # real
    eigenvectors = np.zeros((q, q, num_samples, num_samples), dtype=np.complex128)  # complex
    for band in range(q):
        for idx_x in range(num_samples):
            frac_kx = idx_x / (num_samples - 1)
            for idx_y in range(num_samples):
                frac_ky = idx_y / (num_samples - 1)
                k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                eigvals, eigvecs = np.linalg.eig(H(t3, q, k))
                idx = np.argsort(eigvals)
                eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    # compute Berry fluxes
    berry_fluxes = np.zeros((q, num_samples - 1, num_samples - 1))  # real
    for band in range(q):
        for idx_x in range(num_samples - 1):
            for idx_y in range(num_samples - 1):
                berry_fluxes[band, idx_x, idx_y] = berry_curv(eigenvectors, band, idx_x, idx_y)

    gap_width = (np.min(eigenvalues[1])-np.max(eigenvalues[0]))/(np.max(eigenvalues[0])-np.min(eigenvalues[0]))
    berry_fluc = np.std(berry_fluxes[0])/np.average(berry_fluxes[0])

    return [t3, gap_width, berry_fluc]


if __name__ == "__main__":

    num_samples = 100
    q = 16

    b1 = (2. * np.pi) / 1 * np.array([1 / q, 0])
    b2 = (2. * np.pi) / 1 * np.array([0, 1])
    bvec = np.vstack((b1, b2))

    file = open(f"sp_data/q_{q}_bart.txt", "w")

    ts = list(np.linspace(-0.25, 0, 26))
    results = np.array(Parallel(n_jobs=6)(delayed(my_task)(q, t_temp) for t_temp in ts))
    results = results[results[:, 0].argsort()]  # sort by first column

    for i in range(np.shape(results)[0]):
        file.write(f"{results[i, 0]:.2f}\t{results[i, 1]}\t{results[i, 2]}\n")
    file.close()
