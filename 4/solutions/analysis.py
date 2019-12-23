import argparse
import gzip
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))
plt.rc('text', usetex=True)

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

def f(t, D):
    return 6.0*D*t

def vacf_exp(t, a, b):
    return a * np.exp( -t * b)

def standard_error(x):
    N = len(x)
    ret = np.dot(x,x) - np.sum(x)**2 / N
    ret /= N * (N - 1)
    ret = np.sqrt(ret)
    return ret

def MSD(traj):
    dx = np.zeros((timesteps, N))
    error = np.zeros(timesteps)
    for l in range(N):
        for i in range(1, timesteps):
            for j in range(int(timesteps/(i + 1))):
                start = (i + 1) * j
                end = (i + 1) * (j + 1) - 1
                for k in range(traj.shape[2]):
                    dx[i, l] += np.power(traj[end, l, k] - traj[start, l, k], 2)
            dx[i, l] /= int(timesteps/(i + 1))
            error[i] = standard_error(dx[i,:])
        
        dx_average = np.average(dx, axis=1)
    return dx_average, error

def VACF(v):
    v = v.T
    r = np.zeros(timesteps * 2 - 1)
    for d in range(v.shape[0]):  # number of dimensions
        for part in range(v.shape[1]):  # for each particle
            r += np.convolve(v[d][part], np.flip(v[d][part]), 'full')
    r /= (len(ts)*v.shape[0]*v.shape[1])
    plt.plot(np.linspace(-0.5*timesteps, 0.5*timesteps, len(r)), r/0.3)
    plt.xlabel(r'$t/t_0$')
    plt.ylabel(r'VACF $\langle v(0)v(t)\rangle \cdot \frac{m}{k_\mathrm{B} T}$')
    plt.tight_layout()
    plt.xlim((-30.0, 30.0))
    plt.show()

    D = 0.5 * np.trapz(r,dx=0.5)
    print(D)
    
    param, cov = curve_fit(vacf_exp, ts, r[timesteps - 1:])
    print(*param)
    

dx, error = MSD(traj)
#popt, pcov = curve_fit(f, ts[0:1000], dx[0:1000])
#plt.errorbar(ts, dx, yerr = error,  ecolor='r', label = r'$\langle\Delta x^2\rangle\left(\Delta t\right)$')
plt.plot(ts[0:20], dx[0:20], label = r'$\langle\Delta x^2\rangle\left(\Delta t\right)$')
#plt.plot(ts, f(ts, popt), color = 'k', label=r'Linear fit')
plt.xlabel(r'$\Delta t/t_0$')
plt.legend()
plt.tight_layout()
plt.show()

print(*popt)

VACF(vels)
