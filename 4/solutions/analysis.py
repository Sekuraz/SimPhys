import argparse
import gzip
import pickle
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('id', type=int, help='Simulation id')
args = parser.parse_args()

# Read data from file
datafilename = f'{args.id}.dat.gz'
print(f"Reading simulation data from {datafilename}.")

datafile = gzip.open(datafilename, 'rb')
N, T, GAMMA_LANGEVIN, x, v, ts, Es, Tms, vels, traj = pickle.load(datafile)
datafile.close()

timesteps = len(ts)

def MSD(x):
    dx = np.empty(timesteps)
    for i in range(2, timesteps):
        for j in range((timesteps - 1) // i):
            start =  i * j
            end = i * (j + 1) - 1
            print(f"{i} {j} {start} {end}")
            d = x[end] - x[start]
            dx[i] += np.sum(np.power(d, 2)) / N # mean value over all particles
        dx[i] /= timesteps // i
    return dx

dx = MSD(traj)

plt.plot(ts, dx)
plt.show()
