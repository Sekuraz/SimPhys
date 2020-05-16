import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/RDF_TIP3P.csv'
data= np.loadtxt(dir_file,unpack=False).T

plt.plot(data[0,:], data[1,:])

#plt.xlim((-25.0,30.0))
#plt.ylim((0.0,1.05))
plt.xlabel(r'distance $r$ in nm')
plt.ylabel(r'RDF $g(r)$')
#plt.legend()
plt.tight_layout()
plt.show()
