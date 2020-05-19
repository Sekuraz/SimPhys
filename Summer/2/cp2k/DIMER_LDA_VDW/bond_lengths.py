import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/positions.csv'
data= np.loadtxt(dir_file,unpack=False)

ohbonds = [np.linalg.norm(data[0,:]-data[2,:]), np.linalg.norm(data[1,:]-data[2,:]), np.linalg.norm(data[3,:]-data[5,:]), np.linalg.norm(data[4,:]-data[5,:])]

hbond = np.linalg.norm(data[1,:]-data[5,:])
           
print(ohbonds)
print(np.mean(ohbonds))
print(hbond)
