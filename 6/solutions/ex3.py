import numpy as np
import matplotlib.pyplot as plt
import cising
from tqdm import tqdm
from numba import njit

Ts = np.linspace(1.0, 5.0, 51)
#Ts = np.linspace(1.0, 50, 1)
N_MC = 100000
N_initial = 100000000

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

#Analytical result for the magnetization
def analytical_magnetization(T):
    if (T < 2/np.log(1+np.sqrt(2))):
        ret = (1 - np.sinh(2 / T) ** (-4)) ** 0.125
    else:
        ret = 0.0
    return ret

#Vectorize the function for the analytical magnetization
analytical_magnetization_vectorized = np.vectorize(analytical_magnetization)

#Function for error analysis
@njit
def compute_act_error(x):
    x = np.asarray(x)
    N = len(x)
    xmean = x.mean()
    xvar = x.var()
    acfs = []
    tau_ints = []
    k = 0
    tau_int = 0.5
    while k < 6*tau_int:
        acf = ((x[:N-k]*x[k:]).mean() - xmean*xmean) / xvar
        tau_int += acf
        N_eff = N/(2*tau_int)
        err_tau = tau_int*np.sqrt(12./N_eff)

        acfs.append(acf)
        tau_ints.append(tau_int)
        k += 1

    err_x = np.sqrt(xvar/N*2.0*tau_int)
    return xmean, err_x

#Exact summation
def compute_energy(sigma):
    L = sigma.shape[0]
    #shifted_ixs = range(1,L) + [0]
    shifted_ixs = np.roll(np.arange(0, L), -1)
    E = -(sigma*sigma[shifted_ixs,:]).sum()
    E -= (sigma*sigma[:,shifted_ixs]).sum()
    return E

def compute_magnetization(sigma):
    return sigma.sum()

def exact_sum(L, Ts):
    # we compute the mean energy and magnetization for all
    # temperatures at once!
    ws = np.zeros_like(Ts)
    Es = np.zeros_like(Ts)
    ms = np.zeros_like(Ts)
    # beta is a NumPy array with len(Ts) elements
    beta = 1./Ts

    V = float(L*L)

    sigma = np.ones((L, L), dtype=int)
    # the bit pattern of the integer "state" is used to generate all
    # possible states of sigma
    for state in range(2**(L*L)):
        # read out the bitpattern
        for i in range(L):
            for j in range(L):
                k = i*L + j
                if state & 2**k > 0:
                    sigma[i,j] = 1
                else:
                    sigma[i,j] = -1

        if state%10000==0: print(state)

        # compute energy and magnetization of this state
        E = compute_energy(sigma)
        mu = compute_magnetization(sigma)

        # this is a vector operation, as beta is a vector
        w = np.exp(-beta*E)
        ws += w
        Es += E/V*w
        ms += np.abs(mu)/V*w

    Emeans = Es/ws
    mmeans = ms/ws
    return Emeans, mmeans

##########################
##### Main program
##########################

#Perform the exact summation
Emeans, mmeans = exact_sum(4, Ts)

for L in [64, 16]:
    E_mean = []
    E_error = []
    M_mean = []
    M_error = []
    for T in tqdm(Ts):
        I = cising.IsingModel(1.0/T, L)
        Es = []
        Ms = []
        I.try_many_random_flips(N_initial)
        for i in range(N_MC):
            I.try_many_random_flips(1000)
            Es.append(I.energy())
            Ms.append(I.magnetization())
            
        mean_energy, err_energy = compute_act_error(np.array(Es)/(L*L))
        mean_m, err_m = compute_act_error(np.abs(Ms))
            
        E_mean.append(mean_energy)
        E_error.append(err_energy)
        M_mean.append(np.average(np.abs(Ms)))
        M_error.append(err_m)
    
    plt.errorbar(Ts, E_mean, yerr=E_error, fmt='o-', label='MC for $L$={}'.format(L))
    #plt.plot(Ts, E_mean)
    plt.plot(Ts, Emeans, 'o-', label='exact summation for $L=4$')
    plt.ylabel(r'mean energy $\frac{E}{J}$')
    plt.xlabel(r'dimensionless temperature $\frac{k_\mathrm{B}T}{J}$')
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    plt.errorbar(Ts, M_mean, yerr = M_error, fmt='o-', label='MC for $L$={}'.format(L))
    #plt.plot(Ts, M_mean)
    plt.plot(Ts, mmeans, 'o-', label='exact summation for $L=4$')
    plt.plot(np.linspace(1.0, 5.0, 1000), analytical_magnetization_vectorized(np.linspace(1.0, 5.0, 1000)), label='Analytical result')
    plt.ylabel(r'mean magnetization $\left|\mu\right|$')
    plt.xlabel(r'dimensionless temperature $\frac{k_\mathrm{B}T}{J}$')
    plt.legend()
    plt.tight_layout()
    plt.show()
