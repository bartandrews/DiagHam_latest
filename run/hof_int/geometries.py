import numpy as np

if __name__ == "__main__":

    # Chern number
    C = 1
    # filling factor
    r = 1
    s = 3

    p_list = np.arange(4, 200)
    q_geometries = np.zeros(len(p_list), dtype=object)

    final_geometries = np.zeros(len(p_list), dtype=object)

    for p_idx, p in enumerate(p_list):  # cycle through p

        q_geometries[p_idx] = []
        final_geometries[p_idx] = []

        q = np.abs(C)*p - np.sign(C)  # compute MUC area q
        print(f"(p, q) = ({p}, {q})")

        # factorize q (MUC area) for each p
        for i in range(q):
            if q % (i+1) == 0:
                # print(f"q={q}={i+1}x{int(q/(i+1))}")
                q_geometries[p_idx].append((i+1, int(q/(i+1))))  # create a list of all q-factor tuples
        temp = [tuple(sorted(sub)) for sub in q_geometries[p_idx]]  # sort the tuples in list
        q_geometries[p_idx] = list(set(temp))  # eliminate duplicates in list

        # list of possible, x,y,epsilon,N
        for j_idx, j in enumerate(q_geometries[p_idx]):
            Y = j[0]  # shorter length in Y direction (gauge convention)
            X = j[1]
            for x in range(1, 100):
                for y in range(1, 100):
                    if x*y % s == 0 and np.abs((Y*y)/(X*x) - 1) <= 0.01:  # int N and approx square system
                        N_MUC = x*y
                        N = int(N_MUC/s)
                        if N == 6:  # tractable N
                            epsilon = np.abs((Y*y)/(X*x) - 1)
                            final_geometries[p_idx].append((N, X, Y, x, y, epsilon))

        final_geometries[p_idx] = sorted(final_geometries[p_idx], key=lambda tup: tup[0])
        if len(final_geometries[p_idx]) != 0:
            print(final_geometries[p_idx])
