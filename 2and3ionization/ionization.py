##calculate electron proudction rates for different solar zenith angles (as in task 2)
#  as function of altitude and energy
##calculate photo ionization profiles as function of altitude
##compare total ionization profiles with chapman-profile

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from optical_depth import Tau, Irradiance, watts2photons, irradiance_ph, Xi, n_data, n_all, height, wl

#load photo-ionisation cross sections
phot_ion = pd.read_csv("2and3ionization/phot_ion.dat",sep=r"\s+", skiprows=6)
wl_short = phot_ion.iloc[:,0].to_numpy()

#interpolate photo ionisation cross angles to irradiance data
ion_cs_N2 = np.interp(wl, wl_short, phot_ion.iloc[:,1])
ion_cs_O = np.interp(wl, wl_short, phot_ion.iloc[:,2])
ion_cs_O2 = np.interp(wl, wl_short, phot_ion.iloc[:,3])

#order densities to fit the density file (O-N2-O2) and the irradiance from before
ion_cs = [ion_cs_O, ion_cs_N2, ion_cs_O2]
ion_cs_matrix = np.vstack(ion_cs)

#need to filter out the wavelength that dont ionize
wl_th = np.array([ 91.1, 79.6, 102.6])*1e-9 #ionization threshold wavelength for the different species (m)
filtered_ion_cs = []
for idx, ele in enumerate(wl_th):
    mask = np.where(wl > ele, 0.0, 1.0)
    filtered = mask * ion_cs[idx]
    filtered_ion_cs.append(filtered)


#combined profiles for all species (sums up all seperate ionization rates)
def photo_ion_rate(densities, ion_cs, EUVflux):
    q=0
    for idx,species_cs in enumerate(ion_cs):
        integral = np.sum(EUVflux * species_cs, axis=1)
        q +=  densities[100:,idx]*1e6 * integral 
    return q

total_photoion = []
for X in range(len(Xi)):
    total_photoion.append(photo_ion_rate(n_all, filtered_ion_cs, irradiance_ph[X])) #use SI units (m!)

#make plots
if __name__ == "__main__":
    fig, axs = plt.subplots(1,1)
    ycoord = height[100:]
    for key,ele in enumerate(total_photoion):
        axs.plot(ele, ycoord, label = f"$\chi$ = {Xi[key]}")
    axs.set_xlabel("ionization rate [$m^{-3}$ $s^{-1}$]")
    axs.set_ylabel("height [km]")
    axs.set_xscale("log")
    plt.legend()
    plt.tight_layout()
    plt.show()

plt.plot(wl*1e9, ion_cs_N2, label = "N2")
plt.plot(wl*1e9, ion_cs_O, label="O")
plt.plot(wl*1e9, ion_cs_O2, label = "O2")
plt.xlabel("wavelength [nm]")
plt.title("ionisation cross sections")
plt.legend()
#plt.show()


