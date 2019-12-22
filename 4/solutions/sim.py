import argparse
import gzip
import pickle

import numpy as np

# SYSTEM CONSTANTS
# timestep
DT = 0.01
# length of run
TIME_MAX = 2000.0
# desired temperature
T = 0.3
# total number of particles
N = 50
# Number of spatial dimensions
NDIM = 3
# friction coefficient
GAMMA_LANGEVIN = 0.8
# number of steps to do before the next measurement
MEASUREMENT_STRIDE = 50


parser = argparse.ArgumentParser()
parser.add_argument('id', type=int, help='Simulation id')
args = parser.parse_args()


def compute_temperature(v):
    ret = compute_energy(v) * 2 / (NDIM * N)
    return ret


def compute_energy(v):
    return (v * v).sum() / 2.


def step_vv(x, v, f, dt):
    # update positions
    x += v * dt + 0.5 * f * dt * dt

    # half update of the velocity
    v += 0.5 * f * dt

    # for this excercise no forces from other particles
    f = np.zeros_like(x)

    # second half update of the velocity
    v += 0.5 * f * dt

    return x, v, f


def step_vv_langevin(x, v, f, dt, gamma):

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


# SET UP SYSTEM OR LOAD IT
print("Starting simulation...")
t = 0.0
step = 0

# random particle positions
x = np.random.random((N, NDIM))
v = np.zeros((N, NDIM))

# variables to cumulate data
ts = []
Es = []
Tms = []
vels = []
traj = []


# main loop
f = np.zeros_like(x)


print(f"Simulating until tmax={TIME_MAX}...")

while t < TIME_MAX:
    x, v, f = step_vv_langevin(x, v, f, DT, GAMMA_LANGEVIN)

    t += DT

    if step % MEASUREMENT_STRIDE == 0:
        E = compute_energy(v)
        Tm = compute_temperature(v)
        vels.append(v.copy())
        traj.append(x.copy())
        print(f"t={t}, E={E}, T_m={Tm}")

        ts.append(t)
        Es.append(E)
        Tms.append(Tm)
    step += 1


# at the end of the simulation, write out the final state
datafilename = f'{args.id}.dat.gz'
print(f"Writing simulation data to {datafilename}.")
vels = np.array(vels)
traj = np.array(traj)

datafile = gzip.open(datafilename, 'wb')
pickle.dump([N, T, GAMMA_LANGEVIN, x, v, ts, Es, Tms, vels, traj], datafile)
datafile.close()


print("Finished simulation.")
