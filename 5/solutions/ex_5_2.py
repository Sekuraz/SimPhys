import numpy as np
import matplotlib.pyplot as plt
from sympy import *

def f(x):
    ret = (-2 * x**2 * np.sin(x) * np.cos(x) - 2 * x * np.sin(x)**2)\
        * np.exp(-x**2 * np.sin(x)**2)
    return ret

def exact_integral(a, b):
    x = Symbol('x')
    ret = integrate((-2 * x**2 * sin(x) * cos(x) - 2 * x * sin(x)**2)\
        * exp(-x**2 * sin(x)**2), (x, a, b))
    return ret

def simple_sampling(f, a, b, N):
    dx = (b - a) / N
    ret = 0.0
    for i in range(N):
        ret += f((b - a) * np.random.random_sample() + a)
    ret *= dx
    return ret

x = np.linspace(0.1, 10.0, 1000)
plt.plot(x, f(x))
plt.show()

#print(exact_integral(0.1, 10.0))
for i in range(2, 25):
    print(simple_sampling(f, 0.1, 10.0, 2**i))
