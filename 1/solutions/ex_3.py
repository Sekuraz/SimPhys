#!/usr/bin/env python3

import numpy as np
import scipy.constants


def force(r_ij, m_i, m_j, g):
    return - g * m_i * m_j * r_ij / np.linalg.norm(r_ij) ** 3


def step_euler(x, v, dt, mass, g):
    x += v * dt
    f = forces(x, mass, g)
    # calculate acceleration per coordinate dimension
    f[:, 0] /= mass
    f[:, 1] /= mass
    v += f * dt
    return x, v


def step_euler_symplectic(x, v, dt, mass, g):
    f = forces(x, mass, g)
    f[:, 0] /= mass
    f[:, 1] /= mass
    v += f * dt
    x += v * dt
    return x, v


def step_velocity_verlet(x, v, dt, mass, g):
    f = forces(x, mass, g)
    f[:, 0] /= mass
    f[:, 1] /= mass  # now it's acceleration

    x += v * dt + f / 2 * dt * dt

    # a( t + dt )
    f_t = forces(x, mass, g)
    f_t[:, 0] /= mass
    f_t[:, 1] /= mass  # now it's acceleration

    v += (f + f_t) / 2 * dt

    return x, v


def forces(x, masses, g):
    ret = np.zeros(x.shape)

    for i in range(len(x)):
        # Forces are opposite equal, so we only need this shortened loop
        for j in range(i + 1, len(x)):
            f = force(x[j] - x[i], masses[i], masses[j], g)
            ret[i] -= f
            ret[j] += f

    return ret


def generate_trajectories(steps, dt, x, v, masses, g, integrator):
    # in order not to mess with the original data
    x = x.copy()
    v = v.copy()

    trajectories = [{'x': [x[i][0]], 'y': [x[i][1]]} for i in range(len(names))]
    for s in range(steps):
        integrator(x, v, dt, masses, g)
        for i in range(len(names)):
            trajectories[i]['x'].append(x[i][0])
            trajectories[i]['y'].append(x[i][1])

    return trajectories


def plot_all(trajectories, names):
    for i in range(len(names)):
        plt.plot(trajectories[i]['x'], trajectories[i]['y'], "-", label=names[i].decode("utf-8"))
    plt.legend(loc="upper left", ncol=2)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.show()


def plot_ref(trajectories, reference, moving, label, do_plot=True):
    """
    Plot the moving object in the frame of reference
    :param trajectories: Trajectories object
    :param reference: index of the reference object
    :param moving: index of the moving object
    :param label: The label for the resulting curve
    :param do_plot: show the plot or just add the generated curve to the current plot
    """
    mov_x = [ref - mov for ref, mov in zip(trajectories[reference]['x'], trajectories[moving]['x'])]
    mov_y = [ref - mov for ref, mov in zip(trajectories[reference]['y'], trajectories[moving]['y'])]

    plt.plot(mov_x, mov_y, label=label)
    if do_plot:
        plt.legend()
        plt.show()


def plot_dist(trajectories, b1, b2, dt, label, do_plot=True):
    """
    Plot the distance between b1 and b2 over time
    :param trajectories: Trajectories object
    :param b1: Index of first body
    :param b2: Index of second body
    :param dt: The timestep, this is used for scaling of the axis
    :param label: The label for the resulting curve
    :param do_plot: show the plot or just add the generated curve to the current plot
    """
    dist_x = [ref - mov for ref, mov in zip(trajectories[b1]['x'], trajectories[b2]['x'])]
    dist_y = [ref - mov for ref, mov in zip(trajectories[b1]['y'], trajectories[b2]['y'])]

    time = np.array(range(len(dist_x)), dtype=np.float) * dt
    dist = [np.sqrt(x * x + y * y) for x, y in zip(dist_x, dist_y)]
    plt.semilogy(time, dist, label=label)
    if do_plot:
        plt.legend()
        plt.show()


if __name__ == "__main__":
    import pathlib
    import matplotlib.pyplot as plt

    plt.rc('text', usetex=True)

    current_dir = pathlib.Path(__file__).resolve()
    data = np.load(current_dir.parent.parent.joinpath(
        'files/solar_system.npz'))
    names = data['names']
    x_init = data['x_init']
    v_init = data['v_init']
    masses = data['m']
    g = data['g']

    # convert to list of position and velocity vectors
    x = np.array([np.array(p) for p in zip(x_init[0], x_init[1])])
    v = np.array([np.array(p) for p in zip(v_init[0], v_init[1])])


    # # 3.1
    # plt.xlabel(r"$x$ (AU)")
    # plt.ylabel(r"$y$ (AU)")
    # integrator = step_euler
    # trajectories = generate_trajectories(10000, 0.0001, x, v, masses, g, integrator)
    # # 3.1 full
    # plot_all(trajectories, names)
    # # 3.1 small step moon
    # plot_ref(trajectories, 1, 2, r"$\Delta t = 1e^{-4}$", False)
    #
    # # large step euler
    # trajectories = generate_trajectories(1000, 0.001, x, v, masses, g, integrator)
    # # 3.1 large step moon
    # plot_ref(trajectories, 1, 2, r"$\Delta t = 1e^{-3}$")


    # # 3.2
    # plt.xlabel(r"$x$ (AU)")
    # plt.ylabel(r"$y$ (AU)")
    # trajectories = generate_trajectories(100, 0.01, x, v, masses, g, step_euler)
    # plot_ref(trajectories, 1, 2, "Euler", False)
    #
    # trajectories = generate_trajectories(100, 0.01, x, v, masses, g, step_euler_symplectic)
    # plot_ref(trajectories, 1, 2, "Symplectic euler", False)
    #
    # trajectories = generate_trajectories(100, 0.01, x, v, masses, g, step_velocity_verlet)
    # plot_ref(trajectories, 1, 2, "Velocity verlet")


    # # 3.3
    # plt.xlabel(r"time (years)")
    # plt.ylabel(r"distance (AU)")
    # trajectories = generate_trajectories(1000, 0.01, x, v, masses, g, step_euler)
    # plot_dist(trajectories, 1, 2, 0.01, "Euler", False)
    #
    # trajectories = generate_trajectories(1000, 0.01, x, v, masses, g, step_euler_symplectic)
    # plot_dist(trajectories, 1, 2, 0.01, "Symplectic Euler", False)
    #
    # trajectories = generate_trajectories(1000, 0.01, x, v, masses, g, step_velocity_verlet)
    # plot_dist(trajectories, 1, 2, 0.01, "Velocity Verlet")





