#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

N_random_numbers = 10000
mean_random_numbers = 1.0
sigma_random_numbers = 4.0

N_velocities = 10000
mean_velocities = 0.0
sigma_velocities = 1.0

width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))
plt.rc('text', usetex=True)

# Generator for normally distributed numbers with mean value mean and standard deviation sigma
def box_muller(mean, sigma):
    while True:
        u1 = np.random.random()
        u2 = np.random.random()
        n1 = mean + sigma*np.sqrt(-2*np.log(u1))*np.cos(2*np.pi*u2)
        n2 = mean + sigma*np.sqrt(-2*np.log(u1))*np.sin(2*np.pi*u2)
        yield n1
        yield n2
    
# Function which returns the Gaussian distribution with mean value mean and standard deviation sigma
def gaussian_distribution(x, mean, sigma):
    ret = np.exp(-(x - mean)**2/(2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
    return ret

# Function which returns the threedimensional Maxwell-Boltzmann-distribution with mean value mean and standard deviation sigma
def maxwell_boltzmann_distribution(x, mean, sigma):
    ret = 4*np.pi*x**2*np.exp(-(x - mean)**2/(2*sigma**2)) / (np.power((2*np.pi*sigma**2),1.5))
    return ret
        
random_numbers = box_muller(mean_random_numbers, sigma_random_numbers)
 
x = [next(random_numbers) for i in range(N_random_numbers)]
 
random_velocites = box_muller(0.0, 1.0)
v = [np.array([next(random_velocites),next(random_velocites),next(random_velocites)]) for i in range(N_velocities)] 
v = np.array(v)

y = np.linspace(-15.0, 15.0, 1000)

width = 5.787
height = width*0.8
plt.rc('figure', figsize=(width,height))

plt.hist(x, bins='auto', density = True, label='normalized histogram')
plt.plot(y, gaussian_distribution(y, mean_random_numbers, sigma_random_numbers), label=r'$p(x) = \frac{1}{\sqrt{2\pi\sigma^2}}\cdot\exp\left(-\frac{\left(x-\mu\right)^2}{2\sigma^2}                                                                                                                                       \right)$')
plt.xlabel(r'$x$')
#plt.xlim((-15.0,20.0))
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4))
plt.tight_layout()
plt.show()

z = np.linspace(-5, 5.0, 1000)
plt.hist(v[:,0], bins='auto', density = True, label='normalized histogram for $v_x$')
plt.hist(v[:,1], bins='auto', density = True, label='normalized histogram for $v_y$')
plt.hist(v[:,2], bins='auto', density = True, label='normalized histogram for $v_z$')
plt.plot(z, gaussian_distribution(z, mean_velocities, sigma_velocities), label=r'$p(x) =$')
plt.legend()
plt.show()

width = 5.787
height = width*0.8
plt.rc('figure', figsize=(width,height))
z = np.linspace(0.0, 5.0, 1000)
plt.hist(np.linalg.norm(v/sigma_velocities, axis=1), bins='auto', density = True, label=r'normalized histogram for $v\cdot\sqrt{m\beta}$')
plt.plot(z/sigma_velocities, maxwell_boltzmann_distribution(z, mean_velocities, sigma_velocities)*sigma_velocities, label=r'$p(v)\cdot\sqrt{m\beta}^{-1} = 4\pi\sqrt{\frac{m\beta}{2\pi}}^3v^2\cdot\exp\left(-\beta\frac{mv^2}{2}\right) \cdot\sqrt{m\beta}^{-1}$')
plt.xlabel(r'$v\cdot \sqrt{m\beta}$')
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4))
plt.tight_layout()
plt.show()
