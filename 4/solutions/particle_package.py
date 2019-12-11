import argparse

import matplotlib.pyplot as plt
import numpy as np

# Number of particles
N = 10000

# Time step
DT = 0.01

parser = argparse.ArgumentParser()
parser.add_argument('gamma', type=float, help='Friction coefficient')
parser.add_argument('T', type=float, help='Temperature')
args = parser.parse_args()

# Rounds and integration steps per round
ROUNDS = 3
STEPS = 300


# Velocity Verlet and Langevin step
def step_vv_langevin(x, v, f, dt, T, gamma):
    # INSERT CODE HERE
    return x, v, f


def plot(pos, time, color):
    # Boundaries of the histogram (make them symmetric)
    hist_range = max(-np.amin(pos), np.amax(pos))

    # Sample positions into a histogram
    # INSERT CODE HERE
    # H = np.histogram(...)

    # Calculate bin centers
    bin_centers = (H[1][:-1] + H[1][1:]) / 2

    # INSERT CODE HERE TO PLOT THE ANALYTICAL SOLUTION


# Initial positions (1 coordinate per particle)
# INSERT CODE HERE
# x = np....

# Initial velocities
# INSERT CODE HERE
# v = np...

# Initial forces
# INSERT CODE HERE
# f = np...


plt.figure()
colors = ['k', 'r', 'b']

for i in range(ROUNDS):
    for j in range(STEPS):
        x, v, f = step_vv_langevin(x, v, f, DT, args.T, args.gamma)
    plot(x, (i * STEPS + j) * DT, colors[i])

plt.xlabel("x")
plt.ylabel("P(x, t)")
plt.legend()

plt.show()
