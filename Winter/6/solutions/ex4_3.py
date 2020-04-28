import numpy as np
import matplotlib.pyplot as plt
import cising
from scipy.optimize import curve_fit
from tqdm import tqdm
import pickle

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)


T_C = 2.0/np.log(1.0+np.sqrt(2.0))
N_MC = 1000000
N_initial = 100000000

Ls = [8, 16, 32, 64, 128]
Ts = np.linspace(1.0, 5.0, 51)
Magnetization = dict.fromkeys(Ls)

for L in tqdm(Ls):
    M_mean = []
    for T in tqdm(Ts):
        Ms = []
        I = cising.IsingModel(1/T, L)
        I.try_many_random_flips(N_initial)
        for i in range(N_MC):
            I.try_many_random_flips(100)
            Ms.append(I.magnetization())
        
        M_mean.append(np.average(np.abs(Ms)))
    
    Magnetization[L] = M_mean

state = {}
state['Ls'] = Ls
state['Ts'] = Ts
state['Magnetization'] = Magnetization
    
pickle.dump(state, open("./magnetization.p", "wb"))
