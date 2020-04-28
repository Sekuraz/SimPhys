import numpy as np
import matplotlib.pyplot as plt
from sympy import *

width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))
plt.rc('text', usetex=True)

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
    sample = 0.0
    squared_sum  = 0.0
    for i in range(N):
        sample = f((b - a) * np.random.random_sample() + a)
        ret += sample
        squared_sum += sample ** 2
    error = np.sqrt(squared_sum - (ret ** 2) / N)
    error *= (b - a) / np.sqrt(N * (N-1))
    ret *= dx
    return ret, error

#Create a plot of the function f(x)
x = np.linspace(0.1, 50.0, 1000)
plt.plot(x, f(x))
plt.xlabel(r'$x$')
plt.ylabel(r'$f\left(x\right)$')
plt.tight_layout()
plt.show()

#Integrate the function f(x) analytically
exact_solution = exact_integral(0.1, 50.0)

#Perform the Monte Carlo integration
monte_carlo_solution = [simple_sampling(f, 0.1, 50.0, 2**i) for i in range(2, 21)]
monte_carlo_solution = np.array(monte_carlo_solution)

#Calculate the actual error
error = np.abs(exact_solution - monte_carlo_solution[:,0])

N = np.arange(2, 21)
plt.plot(N, monte_carlo_solution[:,0])
plt.xlabel(r'$\mathrm{log}_2\left(N\right)$')
plt.ylabel(r'Monte Carlo estimate of $\int_{0.1}^{50.0}\mathrm{d}x\,f(x)$')
plt.tight_layout()
plt.show()

plt.plot(N, monte_carlo_solution[:,1], label=r'statistical error')
plt.plot(N, error, label=r'actual error')
plt.xlabel(r'$\mathrm{log}_2\left(N\right)$')
plt.ylabel(r'$\mathrm{Error}$')
plt.legend()
plt.tight_layout()
plt.show()
