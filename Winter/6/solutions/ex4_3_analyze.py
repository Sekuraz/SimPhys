import numpy as np
import matplotlib.pyplot as plt
import pickle

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)


T_C = 2.0/np.log(1.0+np.sqrt(2.0))
a = 0.2

state = pickle.load(open("./magnetization.p", "rb"))
Ls = state['Ls']
Ts = state['Ts']
Magnetization = state['Magnetization']


for L in Ls:
    plt.plot((1 - Ts/T_C)*L, np.array(Magnetization[L])*(L**a), '.', label='$L={}$'.format(L))

plt.xlim((-20.0,20.0))
plt.xlabel(r'$tL$')
plt.ylabel(r'$mL^a$')
plt.tight_layout()
plt.legend()
plt.show()
