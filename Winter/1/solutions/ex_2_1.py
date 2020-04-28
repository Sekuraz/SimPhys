#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)


def force(mass, gravity):
    return np.array([0.0, -mass * gravity])


def step_euler(x, v, dt, mass, gravity, f):
    x += v * dt
    v += f / mass * dt
    return x, v


if __name__ == "__main__":
    gravity = 9.81
    mass = 2.0
    x = np.zeros(2)
    v = np.array([50.0, 50.0])
    dt = 0.1

    f = force(mass, gravity)  # constant

    pos_x = []
    pos_y = []
    while x[1] >= 0:
        pos_x.append(x[0])
        pos_y.append(x[1])

        step_euler(x, v, dt, mass, gravity, f)
    
    plt.plot(pos_x, pos_y, "b-")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y(x)$')
    plt.show()
