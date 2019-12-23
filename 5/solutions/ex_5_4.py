import numpy as np
import matplotlib.pyplot as plt
import itertools

N_TEMPERATURE = 50

def calculate_energy(spin_config, Nx, Ny):
    ret = 0.0
    for nx in range(Nx):
        for ny in range(Ny):
            ret += -0.5 * spin_config[nx * Ny + ny] * (spin_config[((nx - 1) % Nx) * Ny + ny]\
                + spin_config[((nx + 1) % Nx) * Ny + ny] + spin_config[nx * Ny + (ny - 1) % Ny]\
                + spin_config[nx * Ny + (ny + 1) % Ny])
    return ret

def calculate_magnetization(spin_config):
    return np.sum(spin_config)

def ising_mean(Nx, Ny, T):
    Z = 0.0
    mean_energy = 0.0
    mean_magnetization = 0.0
    iterator = itertools.product([-1.0, +1.0], repeat = Nx*Ny)
    for i in iterator:
        m = calculate_magnetization(list(i))
        E = calculate_energy(list(i), Nx, Ny)
        Z += np.exp( -E / T)
        mean_energy += E * np.exp( -E / T)
        mean_magnetization += np.abs(calculate_magnetization(list(i))) * np.exp( -E / T)
    mean_energy /= (Z * Nx * Ny)
    mean_magnetization /= (Z * Nx * Ny)
    return mean_energy, mean_magnetization


E = np.zeros(N_TEMPERATURE)
M = np.zeros(N_TEMPERATURE)

for i in range(N_TEMPERATURE):
    T = 1.0 + i*0.1
    E[i], M[i] = ising_mean(4, 4, T)
    print(E[i])
    print(M[i])

plt.plot(E)
plt.show()

plt.plot(M)
plt.show()
