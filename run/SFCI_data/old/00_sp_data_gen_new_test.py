import numpy as np
from joblib import Parallel, delayed
from time import perf_counter


def H(t3, q, kx, ky, t6=0, t9=0):

    nphi = 1/q
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


def H_eigbasis(t3, q, kx_list, ky_list, t6=0, t9=0):
    evals_arr = np.empty((q, len(kx_list), len(ky_list)))
    evecs_arr = np.empty((q, len(kx_list), len(ky_list), q), dtype="complex_")
    for ikx, kx in enumerate(kx_list):
        for iky, ky in enumerate(ky_list):
            evals, evecs = np.linalg.eigh(H(t3, q, kx, ky, t6, t9))
            for a in range(q):
                evals_arr[a, ikx, iky] = evals[a]
                evecs_arr[a, ikx, iky, :] = evecs[:, a]
    return evals_arr, evecs_arr


def qgt(evecs_comb, band, ikx, iky, dkx_r, dky_r):

    len_kx_list = np.shape(evecs_comb[0])[1]
    len_ky_list = np.shape(evecs_comb[0])[2]

    evec_0 = evecs_comb[0][band, ikx % len_kx_list, iky % len_ky_list, :]
    evec_x = evecs_comb[1][band, ikx % len_kx_list, iky % len_ky_list, :]
    evec_y = evecs_comb[2][band, ikx % len_kx_list, iky % len_ky_list, :]
    grad = {}
    grad.update({"x": (evec_x - evec_0) / dkx_r})
    grad.update({"y": (evec_y - evec_0) / dky_r})

    return np.array(
        [[np.vdot(grad[u], grad[v])
          - np.vdot(grad[u], evec_0) * np.vdot(evec_0, grad[v])
          for u in ["x", "y"]] for v in ["x", "y"]])


def fs_metric(evecs_comb, band, ikx, iky, dkx_r, dky_r):
    return np.real(qgt(evecs_comb, band, ikx, iky, dkx_r, dky_r))


# compute Berry curvature using qgt
def berry_curv(evecs_comb, band, ikx, iky, dkx_r, dky_r):
    return -2*np.imag(qgt(evecs_comb, band, ikx, iky, dkx_r, dky_r))


def TISM(evecs_comb, band, ikx, iky, dkx_r, dky_r):
    return (np.trace(fs_metric(evecs_comb, band, ikx, iky, dkx_r, dky_r))
            - np.abs(berry_curv(evecs_comb, band, ikx, iky, dkx_r, dky_r)[0, 1]))


def band_geom(t3, q, grain, grain_r):
    print(f"t3 = {t3}")

    # define the mesh
    kx_list, ky_list = np.linspace(0, 2*np.pi/q, grain), np.linspace(0, 2*np.pi, grain)
    dkx, dky = (kx_list[1] - kx_list[0]), (ky_list[1] - ky_list[0])  # mesh
    dkx_r, dky_r = (kx_list[1] - kx_list[0])/grain_r, (ky_list[1] - ky_list[0])/grain_r  # reduced mesh (for qgt)

    # compute the eigenbasis
    evals, evecs = H_eigbasis(t3, q, kx_list, ky_list)
    _, evecs_x = H_eigbasis(t3, q, kx_list+dkx_r, ky_list)
    _, evecs_y = H_eigbasis(t3, q, kx_list, ky_list+dky_r)
    evecs_comb = np.array([evecs, evecs_x, evecs_y])

    # initialize the arrays
    tisms = np.empty((len(kx_list)-1, len(ky_list)-1))

    # populate the arrays
    for ikx in range(len(kx_list)-1):
        for iky in range(len(ky_list)-1):
            tisms[ikx][iky] = TISM(evecs_comb, 0, ikx, iky, dkx_r, dky_r)

    # quantum geometry
    tism_int = np.sum(tisms) * dkx * dky / (2 * np.pi)

    return [t3, tism_int]


if __name__ == "__main__":

    t0 = perf_counter()

    q_val = 24
    grain_val = 100
    grain_r_val = 1000
    ts = np.linspace(-0.15, -0.15, 1)

    results = np.array(Parallel(n_jobs=1)(delayed(band_geom)(t3, q_val, grain_val, grain_r_val) for t3 in ts))

    print(results)

    print(f"Total time taken (seconds) = {perf_counter() - t0:.1f}")
