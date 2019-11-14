#!/usr/bin/env python3

import numpy as np
import scipy.linalg


def lj_potential(r_ij):
    return 4 * (1/np.linalg.norm(r_ij)**12 - 1/np.linalg.norm(r_ij)**6)

def lj_force(r_ij):
    return 24.0 * (2.0/np.linalg.norm(r_ij)**14 - 1/np.linalg.norm(r_ij)**8) * r_ij

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.rc('text', usetex=True)

    d = np.zeros((2, 1000))
    d[0,:] = np.linspace(0.85, 2.5, 1000)

    plt.plot(d[0,:], np.apply_along_axis(lj_potential, 0, d), 'b-')
    plt.xlabel(r'$d$')
    plt.ylabel(r'$V_\mathrm{LJ}(d)$')
    plt.show()

    plt.plot(d[0,:], np.apply_along_axis(lj_force, 0, d)[0], 'r-')
    plt.xlabel(r'$d$')
    plt.ylabel(r'$F_{\mathrm{LJ},x}(d)$')
    plt.show()
