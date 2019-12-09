#!/usr/bin/env python3

import pickle
import argparse
import warnings

import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
width = 2*2.893
height = width
plt.rc('figure', figsize=(width,height))

parser = argparse.ArgumentParser()
parser.add_argument('file', help="Path to pickle file.")
parser.add_argument('-t', '--t-eq', type=int, help="Path to pickle file.")

args = parser.parse_args()

with open(args.file, 'rb') as fp:
    data = pickle.load(fp)


def running_average(O, M):
    ret = np.empty_like(O)
    N = len(O)
    for i in range(0, N):
        if (i < M) or (i >= N - M):
            ret[i] = np.nan
        else:
            ret[i] = np.sum(O[i-M:i+M+1]) / (2*M + 1)
    return ret


rdfs = data['rdfs']
#bins = ['rdfs'][0][1]

DENSITY = 0.316
dr = 0.042


def average_rdf(rdfs, t_eq):
    rdfs = rdfs[t_eq:]
    ret = np.zeros_like(rdfs[0])
    for i in range(len(rdfs)):
        ret += rdfs[i]
    ret = ret/len(rdfs)
    for i in range(0, 100):
        ret[i] /= 100*np.pi*DENSITY*(2*dr*(0.8+i*dr) + dr*dr)
    return ret


observables = ['potential_energies', 'kinetic_energies', 'pressures', 'temperatures']
if len(data[observables[0]]) > 1000:
    warnings.warn("Data might be inaccurate due to numerical instabilities!")


def print_teq():
    if not hasattr(args, "t_eq"):
        return
    for o in observables:
        d = data[o]
        d = d[args.t_eq:]
        print("Average equilibrium %s: %s" % (o, sum(d) / len(d)))


print_teq()


def plot_observables():
    fig, axes = plt.subplots(nrows=len(observables), sharex=True)
    for i, o in enumerate(observables):
        d = data[o][args.t_eq:]
        y = list(range(len(d)))
        p = axes[i]
        #p.title.set_text(o)
        p.plot(np.linspace(0, len(d)*0.01*3, len(d)), d, label="Measured time series")
        p.plot(np.linspace(0, len(d)*0.01*3, len(d)),running_average(d, 10), label="Running average for $M=10$")
        p.plot(np.linspace(0, len(d)*0.01*3, len(d)),running_average(d, 100), label="Running average for $M=100$")
    plt.legend(loc='lower center', bbox_to_anchor=(0.5, -1.5), ncol=2)
    axes[0].title.set_text(r'Potential Energy')
    axes[0].set_ylabel(r'$\frac{E_\mathrm{pot}}{\epsilon}$')
    axes[1].title.set_text(r'Kinetic Energy')
    axes[1].set_ylabel(r'$\frac{E_\mathrm{kin}}{\epsilon}$')
    axes[2].title.set_text(r'Pressure')
    axes[2].set_ylabel(r'$\frac{P\cdot\sigma^2}{\epsilon}$')
    axes[3].title.set_text(r'Temperature')
    axes[3].set_xlabel(r'$\frac{t}{\sigma}\cdot\sqrt{\frac{\epsilon}{m}}$')
    axes[3].set_ylabel(r'$\frac{k_\mathrm{B} T}{\epsilon}$')
    plt.tight_layout()
    plt.show()
    
#def plot_rdf():
#    plt.plot(average_rdf(rdfs, args.t_eq))
#    plt.show()

plot_observables()
plt.show()
#plot_rdf()
plt.plot(np.linspace(0.8,5.0,100),average_rdf(rdfs, args.t_eq))
plt.xlabel(r'$\frac{r}{\sigma}$')
plt.ylabel(r'RDF $g(r)$')
#plt.plot(rdfs[0][1])
plt.show()
