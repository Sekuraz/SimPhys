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
    dx = np.zeros(timesteps)
    for i in range(2, timesteps):
        for j in range((timesteps - 1) // i):
            start = i * j
            end = i * (j + 1) - 1
            print(f"{i} {j} {start} {end}")
            d = x[end] - x[start]
            dx[i] += np.sum(np.power(d, 2)) / N  # mean value over all particles
        dx[i] /= timesteps // i
    return dx

dx = MSD(traj)
plt.plot(ts, dx)
plt.show()


def VACF(v):
    v = vels.T
    r = np.zeros(timesteps * 2 - 1)
    for d in range(v.shape[0]):  # number of dimensions
        for part in range(v.shape[1]):  # for each particle
            r += np.convolve(v[d][part], np.flip(v[d][part]), 'full')
    r /= (4000.0*v.shape[0]*v.shape[1])
    plt.plot(r[4000:4050])
    trange = np.linspace(0, 50, 1000)
    plt.plot(trange, r[4000]*np.exp(-0.4*trange))
    plt.xlim((0.0,50.0))
    plt.show()

    D = np.trapz(r[4000:],dx=0.5)
    print(D)

VACF(vels)