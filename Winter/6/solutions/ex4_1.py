import numpy as np
import matplotlib.pyplot as plt
import cising
import itertools
from scipy import optimize
from scipy import interpolate
from tqdm import tqdm

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

Ts = np.linspace(2.0, 2.4, 21)
Ls = [4, 16, 32]
N_MC = 100000000
N_initial = 100000000
U = dict.fromkeys(Ls)
T_intersection = []
 
def binder_parameter(mu):
    ret = 1 - np.average(mu**4)/(3 * np.average(mu**2)**2)
    return ret

for L in Ls:
    U_temp = []
    for T in tqdm(Ts):
        I = cising.IsingModel(1.0/T, L)
        Ms = []
        I.try_many_random_flips(N_initial)
        for i in range(N_MC):
            I.try_random_flip()
            Ms.append(I.magnetization())
        U_temp.append(binder_parameter(np.array(Ms)))
        
    U[L] = U_temp
    
#U[4] = np.linspace(2.0,2.4, 21)
#U[16] = np.linspace(1.0, 3.4, 21)
#U[32] = np.linspace(0.0, 4.4, 21)
    
for L in Ls:
    plt.plot(Ts, U[L], 'o', label='$L$={}'.format(L))
    
    f = interpolate.interp1d(Ts, U[L], kind='quadratic')
    t = np.linspace(2.0, 2.4, 1000)
    plt.plot(t, f(t), label='interpolation for $L$={}'.format(L))
    
for i,j in itertools.combinations(Ls, 2):
    f1 = interpolate.interp1d(Ts, np.array(U[i])-np.array(U[j]), kind='quadratic')    
    T_intersection.append(optimize.bisect(f1,2.0, 2.4))

plt.ylabel(r'Binder parameter $U$')
plt.xlabel(r'dimensionless temperature $\frac{k_\mathrm{B}T}{J}$')
plt.legend()
plt.tight_layout()
plt.show()

print(np.average(T_intersection))
