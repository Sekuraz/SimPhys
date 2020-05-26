import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/pos300K.txt'
data= np.genfromtxt(dir_file,unpack=False)
data = data[:,1:]

hbonds = [np.linalg.norm(data[i * 6 +1,:]-data[i * 6 + 5,:]) for i in range(int(len(data[:,0])/6 -1))]
           
print(np.mean(hbonds))
