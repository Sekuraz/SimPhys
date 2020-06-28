import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import fsolve


dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/positions_1.0.dat'
data = np.loadtxt(dir_file,unpack=False) - np.sqrt(np.pi) * 28.2/2.0

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

r = np.linalg.norm(data, axis=1)
hist, bin_edges = np.histogram(data[10000:], bins = np.linspace(1.0, 28.2, endpoint=True, num=1000))

real_hist = np.cumsum(hist)
real_hist = real_hist / np.amax(real_hist)

#Parameters
R = 28.2
l_B = 1.0
r_0 = 1.0
r = np.linspace(1.0, R, 1000)
lambdas = [1.0]

#Function for charge distribution
def charge_distribution(r, xi, gamma, R_M):
    ret = 1.0 - 1.0 / xi + gamma * np.tan(gamma * np.log(r / R_M)) / xi
    return ret

#Equation which hase to be solved numerically
def equation1(gamma, r_0, xi, R):
    ret = gamma * np.log(r_0 / R) - np.arctan((1 - xi) / gamma) + np.arctan(1.0 / gamma)
    return ret

def equation2(gamma, R):
    ret = R / np.exp(np.arctan(1.0 / gamma) / gamma)
    return ret

#Loop over different values for lambda
for i in lambdas:
    xi = i * l_B
    #Determine parameter numerically:
    gamma = fsolve(equation1, 1.0, args=(r_0, xi, R), maxfev=10000)
    R_M = equation2(gamma, R)
    
    #Plot charge distribution
    plt.plot(r, charge_distribution(r, xi, gamma, R_M), label=r'Poisson-Boltzmann')

plt.xlim(1.0, R)
plt.ylim(0.0, 1.0)

plt.plot(bin_edges[:-1], real_hist, label=r'Simulation')
plt.xscale('log')
plt.xlabel(r'$r/l_\mathrm{B}$')
plt.ylabel(r'$P(r)$')
plt.legend()
plt.tight_layout()
plt.show()
