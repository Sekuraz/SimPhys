import argparse

import matplotlib.pyplot as plt
import numpy as np

# Number of particles
N = 10000

# Number of spatial dimensions 
NDIM = 1

# Time step
DT = 0.01

width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))
plt.rc('text', usetex=True)

parser = argparse.ArgumentParser()
parser.add_argument('gamma', type=float, help='Friction coefficient')
parser.add_argument('T', type=float, help='Temperature')
args = parser.parse_args()

# Rounds and integration steps per round
ROUNDS = 3
STEPS = 300

def gaussian_distribution(x, mean, sigma):
    ret = np.exp(-(x - mean)**2/(2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
    return ret

# Velocity Verlet and Langevin step
def step_vv_langevin(x, v, f, dt, T, gamma):

    # update positions
    x += dt * v * (1 - dt * 0.5 * gamma) + 0.5 * dt * dt * f
    
    # half upate velocity
    v = (v * (1 - 0.5 * gamma * dt) + 0.5 * dt * f) / (1 + 0.5 * dt * gamma)
    
    #calculate new random force
    f = np.random.random_sample(np.shape(x))
    f -= 0.5
    f *= np.sqrt(12 * 2 * T * gamma / dt)
    
    # second half update of the velocity
    v += 0.5 * dt * f / (1 + 0.5 * dt * gamma)
    
    return x, v, f


def plot(pos, time, color, T, gamma):
    # Boundaries of the histogram (make them symmetric)
    hist_range = max(-np.amin(pos), np.amax(pos))

    # Sample positions into a histogram
    H = np.histogram(pos, bins=200, density=True)

    # Calculate bin centers
    bin_centers = (H[1][:-1] + H[1][1:]) / 2
    plt.plot(bin_centers, H[0], label='Numerical solution')

    # Plot the analytical solution
    D = T / gamma
    x = np.linspace(-10.0, 10.0, 1000)
    plt.plot(x, gaussian_distribution(x, 0.0, np.sqrt(2 * D * time)), linestyle='-.', label='Analytical solution')


# Initial positions (NDIM coordinate per particle)
x = np.zeros((N,NDIM))

# Initial velocities
v = np.random.normal(0.0, np.sqrt(args.T), (N,NDIM))

# Initial forces
f = np.zeros((N,NDIM))


plt.figure()
colors = ['k', 'r', 'b']

for i in range(ROUNDS):
    for j in range(STEPS):
        x, v, f = step_vv_langevin(x, v, f, DT, args.T, args.gamma)
    plot(x, (i * STEPS + j) * DT, colors[i], args.T, args.gamma)

plt.xlim((-6,6))
plt.xlabel(r'$x/x_0$')
plt.ylabel(r'$P(x, t)\cdot x_0$')
plt.legend()
plt.tight_layout()
plt.show()
