import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv("C:/Users/skief/Documents/UiT/semester8/space physics/progrTasks/1atmosphere/MSIS.dat",sep=r"\s+", skiprows=16)


#calculate scaleheight from data --> H= k_B*T/(m*g)
k_B = 1.38 * 10e-23     #[m^2*kg/(K*s^2)] boltzmanns constant
g = 9.81 #[m/s] gravitational acceleration
#O weighs 16amu or g/mol
#N2 weighs 28amu or g/mol
#O2 weighs 32amu or g/mol

def calc_H(data_file, height):
     #nr density all molecules at certain height
    n = data_file.iloc[height,1]+ data_file.iloc[height,2]+ data_file.iloc[height,3]
    print(n)
    #average molecule weight
    m = (data_file.iloc[height,1]*16 + data_file.iloc[height,2]*28 + data_file.iloc[height,3]*32)/n *1.66054 * 10e-27 #[kg]
    #scale height
    H = (k_B * data_file.iloc[height,5]) / (m * g)
    return H

H_0 = calc_H(data, 0)
H_10 = calc_H(data, 10)
print(f"scale height at 0km:{H_0}m, at 10km: {H_10}m")

#scale height for single species
#constant g

