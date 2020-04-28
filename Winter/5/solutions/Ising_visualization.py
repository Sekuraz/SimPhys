import numpy as np
import matplotlib.pyplot as plt
from numba import njit
#from matplotlib.animation import FuncAnimation

NX = 1000
NY = 1000
N_STEPS = 1000000000
T = 1.0

@njit
def calculate_energy(spin_config, Nx, Ny):
    ret = 0.0
    for nx in range(Nx):
        for ny in range(Ny):
            ret += -0.5 * spin_config[nx * Ny + ny] * (spin_config[((nx - 1) % Nx) * Ny + ny]\
                + spin_config[((nx + 1) % Nx) * Ny + ny] + spin_config[nx * Ny + (ny - 1) % Ny]\
                + spin_config[nx * Ny + (ny + 1) % Ny])
    return ret

@njit
def ising_monte_carlo(Nx, Ny, T, N_STEPS):
    N_SPINS = Nx * Ny
    spin_config = [(-1) ** np.random.randint(0, high=2) for i in range(N_SPINS)]
    energy = calculate_energy(spin_config, Nx, Ny)
    #yield spin_config
    while True:
        for i in range(N_STEPS):
            #Choose a random spin to flip
            random_spin = np.random.randint(0, high=N_SPINS)
        
            #Calculate the coordinates of the spin
            ny = random_spin % Ny
            nx = int((random_spin - ny) / Ny)
        
            #Calculate the new energy for this configuration
            energy_trial = energy + 2 * spin_config[nx * Ny + ny] * (spin_config[((nx - 1) % Nx) * Ny + ny]\
            + spin_config[((nx + 1) % Nx) * Ny + ny] + spin_config[nx * Ny + (ny - 1) % Ny]\
            + spin_config[nx * Ny + (ny + 1) % Ny])
        
            #Decide if new spin configuration is accepted
            if np.random.random_sample() < min(1, np.exp( (energy - energy_trial)/ T)):
                spin_config[random_spin] *= -1
                energy = energy_trial
        
        yield spin_config
        
#Main Program
monte_carlo_step = ising_monte_carlo(NX, NY, T, N_STEPS)

#next(monte_carlo_step)
#for i in range(10):
    #next(monte_carlo_step)
data = next(monte_carlo_step)
   # print(data)
plt.imshow(np.reshape(data, (NX, NY)), cmap = 'binary')
plt.show()
#    plt.savefig('/home/dbeyer/Ising/ising{}'.format(i), format = 'pdf')
