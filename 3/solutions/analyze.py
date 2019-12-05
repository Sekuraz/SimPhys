#!/usr/bin/env python3

import pickle
import argparse
import warnings

import matplotlib.pyplot as plt
import numpy as np

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


rdfs = ['rdfs'][:][0]
bins = ['rdfs'][0][1]


def average_rdf(rdfs, t_eq):
    rdfs = rdfs[t_eq:]
    ret = np.zeros_like(rdfs[0])
    for i in range(len(rdfs)):
        ret += rdfs[i]
    ret /= len(rdfs)
    return ret


observables = ['energies', 'pressures', 'temperatures']
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
        p.title.set_text(o)
        p.plot(d, label="raw")
        p.plot(running_average(d, 10), label="running 10")
        p.plot(running_average(d, 100), label="running 100")
    plt.legend()

    plt.show()


plot_observables()
