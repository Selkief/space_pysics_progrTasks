##task 1
##Make functions that calculate the optical depth as a function of altitude and wavelength
# for vertical incidence, and for variable zenith-angle of the incident light

#can either use two (or three) different functions for the 3 "levels" of incidence angle
#or just use the most complicated one which should be reduced to the simpler ones automatically
#when using lower angles of incidence

#optical depth for vertical incidence:
#tau0 = sum( crosssections* integral( n(z)dz ) )

#make plots with altitude variations of thermospheric densities
#make plots with wavelength variations of cross sections

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

n_data = pd.read_csv("C:/Users/skief/Documents/UiT/semester8/space physics/progrTasks/1atmosphere/MSIS.dat",sep=r"\s+", skiprows=16)
I_data = pd.read_csv("2and3ionization/fism_daily_hr19990216.dat")
I_data.rename(columns={"wavelength (nm)": "wavelength", "irradiance (W/m^2/nm)" : "irradiance"}, inplace=True)

#number densities
n_all = n_data.iloc[:,1] + n_data.iloc[:,2] + n_data.iloc[:,3]

plt.plot(n_data.iloc[90:,1], n_data.iloc[90:,0], label="O")
plt.plot(n_data.iloc[90:,2], n_data.iloc[90:,0], label="N2")
plt.plot(n_data.iloc[90:,3], n_data.iloc[90:,0], label="O2")
plt.plot(n_all[90:], n_data.iloc[90:,0], label = "all species")
plt.xscale("log")
plt.title("number densities in the thermosphere")
plt.xlabel("n [cm^-3]")
plt.ylabel("height [km]")
plt.legend()
plt.show()

#wavelength variations
