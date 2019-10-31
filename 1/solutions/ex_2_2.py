#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


import ex_2_1


def force(mass, gravity, v, gamma, v_0):
    ret = np.array([0.0, -mass * gravity])
    ret -= gamma * (v - v_0)
    return ret


def step(x, v, dt, mass, gravity, gamma, v_0):
    f = force(mass, gravity, v, gamma, v_0)
    x += v * dt
    v += f / mass * dt

    return x, v


def trajectory(v_w, friction):
    pos_x = []
    pos_y = []

    gravity = 9.81
    mass = 2.0
    x = np.zeros(2)
    v = np.array([50.0, 50.0])
    dt = 0.1

    gamma = 0.1 if friction else 0
    v_0 = np.array([v_w, 0.0])

    while x[1] >= 0:
        pos_x.append(x[0])
        pos_y.append(x[1])

        step(x, v, dt, mass, gravity, gamma, v_0)
    return pos_x, pos_y


if __name__ == "__main__":
    # x, y = trajectory(0, False)
    # plt.plot(x, y, "-", label="y = 0")
    # x, y = trajectory(0, True)
    # plt.plot(x, y, "-", label="y = 0.1")
    # x, y = trajectory(-50, True)
    # plt.plot(x, y, "-", label="y = 0.1, v_w = -50")
    # plt.legend()
    # plt.show()

    for v_w in range(0, 201, 25):
        x, y = trajectory(-v_w, True)
        plt.plot(x, y, "-", label="v_w = -%d" % v_w)

    plt.legend()
    plt.show()
