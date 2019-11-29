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
    pass


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
