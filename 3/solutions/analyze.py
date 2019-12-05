#!/usr/bin/env python3

import pickle
import argparse

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
            ret[i] =  np.sum(O[i-M:i+M+1]) / (2*M + 1)           
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

def print_teq():
    if not hasattr(args, "t_eq"):
        return
    for o in observables:
        d = data[o]
        if len(d) < args.t_eq:
            return
        d = d[args.t_eq:]
        print("Average equilibrium %s: %s" % (o, sum(d) / len(d)))


print_teq()
