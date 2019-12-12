#!/usr/bin/env python3

import time
import numpy as np
import matplotlib.pyplot as plt

n_steps = 1000
n_walks = 100
x0 = 0.0

plt.rc('text', usetex=True)


def linear_congruental_generator():
    a = 1103515245
    c = 13245
    m = 2**32
    #seed = 123647
    seed = int(time.time()*10**7)
    X = seed
    while True:
        X = (a*X + c) % m
        yield X/m
        
def random_walk(N, x0):
    ret = np.empty(N+1)
    ret[0] = x0
    k = linear_congruental_generator()
    for i in range(N):
        ret[i+1] = ret[i] + next(k) - 0.5
    return ret


walks = []

for l in range(n_walks):
    walks.append(random_walk(n_steps, x0))
    
for i in range(len(walks)):
    plt.plot(walks[i])

plt.xlabel(r'$n$')
plt.ylabel(r'$x$')
plt.show()
