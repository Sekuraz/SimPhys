#!/usr/bin/env python3

import numpy as np
import scipy.constants


def force(r_ij, m_i, m_j, g):
    return - g * m_i * m_j * r_ij / np.linalg.norm(r_ij) ** 3


def step_euler(x, v, dt, mass, g):
    f = forces(x, mass, g)
    x += v * dt
    v += f / mass * dt
    return x, v


def forces(x, masses, g):
    ret = np.zeros(x.shape)

    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            f = force(x[j] - x[i], masses[i], masses[j], g)
            ret[i] += f
            ret[j] -= f

    return ret


if __name__ == "__main__":
    import pathlib
    import matplotlib.pyplot as plt

    current_dir = pathlib.Path(__file__).resolve()
    data = np.load(current_dir.parent.parent.joinpath(
        'files/solar_system.npz'))
    names = data['names']
    x_init = data['x_init']
    v_init = data['v_init']
    masses = data['m']
    g = data['g']

    x = np.array([np.array(p) for p in zip(x_init[0], x_init[1])])
    v = np.array([np.array(p) for p in zip(v_init[0], v_init[1])])


    dt = 0.0001
    trajectories = []
    for s in range(10000):




    print(names)
    print(x)
    print(v)