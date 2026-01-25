import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv("C:/Users/skief/Documents/UiT/semester8/space physics/progrTasks/1atmosphere/MSIS.dat",sep=r"\s+", skiprows=16)


#calculate scaleheight from data --> H= k_B*T/(m*g)
k_B = 1.38 * 10e-23     #[m^2*kg/(K*s^2)] boltzmanns constant
g0, g = 9.82, 0 #[m/s] gravitational acceleration
R_E = 6378 *10e3 #[m] earth radius
#O weighs 16amu or g/mol
#N2 weighs 28amu or g/mol
#O2 weighs 32amu or g/mol

def calc_H1(data_file, height):
     #nr density all molecules at certain height
    n = data_file.iloc[height,1]+ data_file.iloc[height,2]+ data_file.iloc[height,3]
    #average molecule weight
    m = (data_file.iloc[height,1]*16 + data_file.iloc[height,2]*28 + data_file.iloc[height,3]*32)/n *1.66054 * 10e-27 #[kg]
    #scale height
    H = (k_B * data_file.iloc[height,5]) / (m * g0)
    return H

H_0 = calc_H1(data, 0)
H_10 = calc_H1(data, 10)
print(f"scale height at 0km:{H_0}m, at 10km: {H_10}m")

#scale height for single species

def calc_H2(data_file, height, g, species):
    if species == "O":
        m = 16 * 1.66054 * 10e-27
    elif species == "N2":
        m = 28 * 1.66054 * 10e-27
    elif species == "O2":
        m = 32 * 1.66054 * 10e-27
    
    if g =="g_z":
        g = g0 * R_E**2 /(R_E + height*10e-3)**2
    else:
        g = g0
    H = (k_B * data_file.iloc[height, 5]) / (m * g)
    return H

#constant g:
H_120_O = calc_H2(data, 120, g0, "O")
H_120_N2 = calc_H2(data, 120, g0, "N2")
H_120_O2 = calc_H2(data, 120, g0, "O2")
H_600_O = calc_H2(data, 600, g0, "O")
H_600_N2 = calc_H2(data, 600, g0, "N2")
H_600_O2 = calc_H2(data, 600, g0, "O2")

#non-constant g:
H_120_gz_O = calc_H2(data, 120, g, "O")
H_120_gz_N2 = calc_H2(data, 120, g, "N2")
H_120_gz_O2 = calc_H2(data, 120, g, "O2")
H_600_gz_O = calc_H2(data, 600, g, "O")
H_600_gz_N2 = calc_H2(data, 600, g, "N2")
H_600_gz_O2 = calc_H2(data, 600, g, "O2")



print(f"scale height 120km cst g for O {H_120_O}m for N2 {H_120_N2}m for O2 {H_120_O2}m")
print(f"scale height 600km cst g for O {H_600_O}m for N2 {H_600_N2}m for O2 {H_600_O2}m")

print(f"scale height 120km non-cst g for O {H_120_gz_O}m for N2 {H_120_gz_N2}m for O2 {H_120_gz_O2}m")
print(f"scale height 600km non-cst g for O {H_600_gz_O}m for N2 {H_600_gz_N2}m for O2 {H_600_gz_O2}m")


