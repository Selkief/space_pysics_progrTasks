#important reactions: with e, O+, NO+, O2+, N2+ and NO
#(10 reactions, which are on sl 234 combined notes )
#must write equation of change in density for each particle form reaction table

#coupled ODEs for each altitude
#dn/dt = P - L (sl 218)
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#load data with all the densities
msis = pd.read_csv("C:/Users/skief/Documents/UiT/semester8/space physics/progrTasks/1atmosphere/MSIS.dat",sep=r"\s+", skiprows=16)
iri = pd.read_csv("IRI.dat", sep=r"\s+", skiprows=45)

n_e_data = iri.iloc[:,1].to_numpy()*1e6 #density in m^-3
n_O_data = msis.iloc[100:,1].to_numpy()*1e6 #density in m^-3
n_Oplus_data = iri.iloc[:,4].to_numpy()/100 * n_e_data #density m^-3
n_O2_data = msis.iloc[100:,3].to_numpy()*1e6 #density m^-3
n_O2plus_data = iri.iloc[:,7].to_numpy()/100 * n_e_data #density m^-3
n_N2_data = msis.iloc[100:,2].to_numpy()*1e6 #density m^-3
###dummy values bc i didnt find the values, must replace!!!!
n_N2plus_data = np.zeros_like(n_N2_data)
n_NOplus_data = iri.iloc[:,8].to_numpy()/100 * n_e_data #density m^-3
###dummy values must replace!!!!
n_NO_data = np.zeros_like(n_NOplus_data)

##this equation is still wrong as all the n's should be q's
q_N2plus = n_e_data * 0.92 * n_N2_data / (0.92 * n_N2_data + n_O2_data + 0.56 * n_O_data)

#temperatures
Te = iri.iloc[:,3].to_numpy()
Ti = iri.iloc[:,2].to_numpy()
Tn = msis.iloc[100:,5].to_numpy()
Tr = (Tn + Ti)/2

kB = 1.380649e-23 #[J/K]
#reaction coefficients
alpha1 = 4.2e-13 * np.power((Te/300), -0.85)
alpha2 = 1.9e-13 * np.power((Te/300), -0.5)
alpha3 = 1.8e-13 * np.power((Te/300), -0.39)
alphar = 3.7e-18 * np.power((Te/250), -0.7)
k1 = 1.3e-18 * np.power((Tr/300),-0.5) + 6.8e-16 * np.exp(-1.4e-19*(kB*Tr))
k2 = 2e-17 * np.power((Tr/300), -0.4)
k3 = 4.4e-16
k4 = 5e-22
k5 = 1.4e-16 * np.power((Tr/300), -0.44)
k6 = 5e-17 * np.power((Tr/300), -0.8)

#must change inital conditions --> seperate function?

#define ODEs for a certain height 
def initial_cond(height):
    #initial densities for a certain height
    n_e = n_e_data[ height - 100 ]
    n_Oplus = n_Oplus_data[ height - 100]
    n_O2plus = n_O2plus_data[ height - 100]
    n_N2plus = n_N2plus_data[ height - 100]
    n_NO = n_NO_data[ height - 100]
    n_NOplus = n_NOplus_data[ height - 100]

    return n_e, n_Oplus, n_O2plus, n_N2plus, n_NO, n_NOplus

def reactions(t, z, q, height):

    n_e, n_Oplus, n_O2plus, n_N2plus, n_NOplus, n_NO = z

    idx = height - 100
    n_O = n_O_data[ idx ]
    n_O2 = n_O2_data[ idx ]
    n_N2 = n_N2_data[ idx]

    alpha_1 = alpha1[idx]
    alpha_2 = alpha2[idx]
    alpha_3 = alpha3[idx]
    alpha_r = alphar[idx]
    k_1 = k1[idx]
    k_2 = k2[idx]
    k_3 = k3
    k_4 = k4
    k_5 = k5[idx]
    k_6 = k6[idx]
    
    #set up ODEs !!!MUSt adjust ion_rates in equations?!!
    dn_e = q - n_e * (alpha_1 * n_NOplus + alpha_2 * n_O2plus + alpha_3 * n_N2plus + alpha_r * n_Oplus)
    dn_Oplus = q + k_5 * n_O * n_N2plus - n_Oplus * (alpha_r * n_e + k_1 * n_N2 + k_2 * n_O2)
    dn_O2plus = q + k_2 * n_Oplus * n_O2 + k_6 * n_N2plus * n_O2 - n_O2plus * (alpha_2 * n_e + k_3 * n_NO + k_4 * n_N2)
    dn_N2plus = q - n_N2plus * (alpha_3 * n_e + k_5 * n_O + k_6 * n_O2)
    dn_NOplus = q + k_1 * n_Oplus * n_N2 + k_3 * n_O2plus * n_NO + k_4 * n_O2plus * n_N2 - alpha_1 * n_NOplus *n_e
    dn_NO = q + k_4 * n_O2plus * n_N2 - k_3 * n_O2plus * n_NO

    return [dn_e, dn_Oplus, dn_O2plus, dn_N2plus, dn_NOplus, dn_NO]

IC = initial_cond(110)
sol = solve_ivp(reactions, [0,360], IC, method = "RK45", args=(1e9, 250))


labels = ["e", "O+", "O2+", "N2+", "NO+", "NO"]
fig, axs = plt.subplots(6,1, figsize=(18,18))
for idx, ele in enumerate(sol.y):
    axs[idx].plot(sol.t, ele.T, label = labels[idx])
    axs[idx].set_xlabel("t")
    #axs[idx].plot_legend()
plt.suptitle("densities at height 110km")
plt.tight_layout()
plt.show()


