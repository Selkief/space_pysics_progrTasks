import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv("C:/Users/skief/Documents/UiT/semester8/space physics/progrTasks/1atmosphere/MSIS.dat",sep=r"\s+", skiprows=16)

#calculate scaleheight from data --> H= k_B*T/(m*g)
k_B = 1.38 * 10e-23     #[m^2*kg/(K*s^2)] boltzmanns constant
g0 = 9.82 #[m/s] gravitational acceleration
R_E = 6378 *10e3 #[m] earth radius
#O weighs 16amu or g/mol
#N2 weighs 28amu or g/mol
#O2 weighs 32amu or g/mol

#define functions
#scale height for all species
def calc_Hall(data_file, height):
     #nr density all molecules at certain height
    n = data_file.iloc[height,1]+ data_file.iloc[height,2]+ data_file.iloc[height,3]
    #average molecule weight
    m = (data_file.iloc[height,1]*16 + data_file.iloc[height,2]*28 + data_file.iloc[height,3]*32)/n *1.66054 * 10e-27 #[kg]
    #scale height
    H = (k_B * data_file.iloc[height,5]) / (m * g0)
    return H

#scale height for single species, cst g or varying g
def calc_Hsingle(data_file, height, g, species):
    if species == "O":
        m = 16 * 1.66054 * 10e-27
    elif species == "N2":
        m = 28 * 1.66054 * 10e-27
    elif species == "O2":
        m = 32 * 1.66054 * 10e-27
    
    if g =="g_z":
        g = g0 * R_E**2 /(R_E + height*10e3)**2
    else:
        g = g0
    H = (k_B * data_file.iloc[height, 5]) / (m * g)
    return H

#experimental scale height H_x = -n_x/(dn_x/dz)
def exp_H(data_file, height, species):
    if species == "O":
        column = 1
    elif species == "N2":
        column = 2
    elif species == "O2":
        column = 3
    H = []
    for i in range(1,height):
        Hz = (-data_file.iloc[i,column])/(data_file.iloc[i, column]-data_file.iloc[i-1, column]/10e3)
        H.append(float(Hz))
    return np.nansum(H)

#calculate results
#scaleheight for all species
H_0 = calc_Hall(data, 0)
H_10 = calc_Hall(data, 10)

#scaleheight for single species
#constant g:
H_120_O = calc_Hsingle(data, 120, g0, "O")
H_120_N2 = calc_Hsingle(data, 120, g0, "N2")
H_120_O2 = calc_Hsingle(data, 120, g0, "O2")
H_600_O = calc_Hsingle(data, 600, g0, "O")
H_600_N2 = calc_Hsingle(data, 600, g0, "N2")
H_600_O2 = calc_Hsingle(data, 600, g0, "O2")
H_cstg =[H_120_O, H_120_N2, H_120_O2, H_600_O, H_600_N2, H_600_O2]

#non-constant g:
H_120_gz_O = calc_Hsingle(data, 120, "g_z", "O")
H_120_gz_N2 = calc_Hsingle(data, 120, "g_z", "N2")
H_120_gz_O2 = calc_Hsingle(data, 120, "g_z", "O2")
H_600_gz_O = calc_Hsingle(data, 600, "g_z", "O")
H_600_gz_N2 = calc_Hsingle(data, 600, "g_z", "N2")
H_600_gz_O2 = calc_Hsingle(data, 600, "g_z", "O2")
H_gz = [H_120_gz_O, H_120_gz_N2, H_120_gz_O2, H_600_gz_O, H_600_gz_N2, H_600_gz_O2]

#difference between scaleheight with constant g and varying g (dependend on z)
diff = []
for i in range(6):
   diff.append(float(H_gz[i] - H_cstg[i]))

#exp scaleheight: #doesnt work yet
expH_120_O = exp_H(data, 120, "O")
exp_H_600_O = exp_H(data, 600, "O")

#print solutions for part 1
print(f"scale height at 0km:{H_0}m, at 10km: {H_10}m")
print(f"scale height 120km cst g for O {H_120_O:.2f}m for N2 {H_120_N2:.2f}m for O2 {H_120_O2:.2f}m")
print(f"scale height 600km cst g for O {H_600_O:.2f}m for N2 {H_600_N2:.2f}m for O2 {H_600_O2:.2f}m")
print(f"scale height 120km non-cst g for O {H_120_gz_O:.2f}m for N2 {H_120_gz_N2:.2f}m for O2 {H_120_gz_O2:.2f}m")
print(f"scale height 600km non-cst g for O {H_600_gz_O:.2f}m for N2 {H_600_gz_N2:.2f}m for O2 {H_600_gz_O2:.2f}m")
print(f"difference between constant and non-constant g (in m):{diff}")

print(f"experimental scaleheight for O at 120km {expH_120_O} and at 600km {exp_H_600_O}")

#part 2: altitude variations of scaleheight
H_all = []
H_O = []
H_N2 = []
H_O2 = []

for k in range(600):
    H_all.append(float(calc_Hall(data, k))*10e-3)
    H_O.append(float(calc_Hsingle(data, k, g0, "O"))*10e-3)
    H_N2.append(float(calc_Hsingle(data, k, g0, "N2"))*10e-3)
    H_O2.append((float(calc_Hsingle(data, k, g0, "O2"))*10e-3))
#define some axis for height
z = np.arange(600)

#plot results
fig, axs = plt.subplots(1,2)
axs[0].plot(H_all,z)
axs[0].grid(True)
axs[0].set_title("scaleheight for all species")
axs[0].set_xlabel("scaleheight in km")
axs[0].set_ylabel("height in km")
axs[1].plot(H_O, z)
axs[1].grid(True)
axs[1].set_title("scaleheight for O")
axs[1].set_xlabel("scaleheight in km")
axs[1].set_ylabel("height in km")
plt.tight_layout()
plt.show()

fig, axs = plt.subplots(1,2)
axs[0].plot(H_N2,z)
axs[0].grid(True)
axs[0].set_title("scaleheight for N2")
axs[0].set_xlabel("scaleheight in km")
axs[0].set_ylabel("height in km")
axs[1].plot(H_O2, z)
axs[1].grid(True)
axs[1].set_title("scaleheight for O2")
axs[1].set_xlabel("scaleheight in km")
axs[1].set_ylabel("height in km")
plt.tight_layout()
plt.show()