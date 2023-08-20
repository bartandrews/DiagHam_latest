import numpy as np
from fractions import Fraction
from math import gcd
import matplotlib.pyplot as plt


def chern(pval, qval):

    nphi = Fraction(pval, qval)
    p = nphi.numerator
    q = nphi.denominator

    # determine r and s
    sr_list, tr_list = [], []

    for r in range(q+1):
        if q % 2 == 0 and r == q/2:
            continue
        for tr in range(-int(q/2), int(q/2)+1):
            for sr in range(-q, q+1):
                if r == q*sr + p*tr:
                    print(sr, tr)
                    sr_list.append(sr)
                    tr_list.append(tr)
                    break
            else:
                continue  # only executed if the inner loop did NOT break
            break  # only executed if the inner loop DID break

    # print(tr_list)

    Chern_list = []
    if q % 2 != 0:
        numb_band_groups = q
    else:
        numb_band_groups = q-1

    for i in range(numb_band_groups):
        Chern_list.append(tr_list[i+1] - tr_list[i])

    if q % 2 == 0:
        Chern_list.insert(q//2-1, Chern_list[q//2-1])

    print("sum cherns = ", np.sum(Chern_list))

    return Chern_list


if __name__ == "__main__":

    print(chern(2, 99))

    # fig = plt.figure(figsize=(6, 9))
    # ax = fig.add_subplot(111)
    # sc = ax.scatter([[1, 1, 1], [2, 2, 2]], [[1, 2, 3], [1, 2, 3]], c=[[1, -40, 0], [1, -20, 0]], cmap='viridis')
    # cbar = plt.colorbar(sc)
    #
    # plt.show()
