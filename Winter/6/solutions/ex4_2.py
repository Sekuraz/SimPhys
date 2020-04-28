import numpy as np
import matplotlib.pyplot as plt
import cising
from scipy.optimize import curve_fit
from tqdm import tqdm

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

def linear_fit(x, a, b):
    return a*x + b

Beta_C = np.log(1.0+np.sqrt(2.0))/2.0
#Beta_C = 1./2.2654037621400978
N_MC = 100000000
N_initial = 100000000

M_mean = []
Ls = [8, 16, 32, 64, 128]

for L in tqdm(Ls):
    Ms = []
    I = cising.IsingModel(Beta_C, L)
    I.try_many_random_flips(N_initial)
    for i in tqdm(range(N_MC)):
        I.try_random_flip()
        Ms.append(I.magnetization())
        
    M_mean.append(np.average(np.abs(Ms)))

popt, pcov = curve_fit(linear_fit, np.log(np.array(Ls)), np.log(np.array(M_mean)))
print(popt)

plt.plot(np.log(Ls), np.log(M_mean), 'o', label=r'simulation results')
plt.plot(np.log(Ls), linear_fit(np.log(Ls), *popt), label=r'linear fit')
plt.ylabel(r'$\mathrm{log}\left(m\left(T=T_\mathrm{c}\right)\right)$')
plt.xlabel(r'$\mathrm{log}\left(L\right)$')
plt.tight_layout()
plt.legend()
plt.show()
