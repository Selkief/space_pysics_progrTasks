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
import matplotlib.colors as mcolors
import pandas as pd
import seaborn as sns

n_data = pd.read_csv("C:/Users/skief/Documents/UiT/semester8/space physics/progrTasks/1atmosphere/MSIS.dat",sep=r"\s+", skiprows=16)
I_data = pd.read_csv("2and3ionization/fism_daily_hr19990216.dat")
I_data.rename(columns={"wavelength (nm)": "wavelength", "irradiance (W/m^2/nm)" : "irradiance"}, inplace=True)
cs_data = pd.read_csv("2and3ionization/phot_abs.dat",sep=r"\s+", skiprows=8)

#interpolate cross sections to irradiance-wavelengths
cs_N2 = np.interp( I_data["wavelength"]*1e-9, cs_data.iloc[:,0], cs_data.iloc[:,1] )
cs_O = np.interp( I_data["wavelength"]*1e-9, cs_data.iloc[:,0], cs_data.iloc[:,2] )
cs_O2 = np.interp( I_data["wavelength"]*1e-9, cs_data.iloc[:,0], cs_data.iloc[:,3] )
cs = [cs_O, cs_N2, cs_O2] #ordered to fit the order in the n_data file [m^2]
print(cs_N2)



#relevant angles for Xi
Xi = [0, 15, 30, 45, 60, 75, 85] #[degrees]

def Tau(altitude, angle):
    tau = 0
    for i in range(3):
        integral = 0
        for k in range(altitude, len(n_data.iloc[:,0])):
            integral += n_data.iloc[k,i+1]*1e6*1000
        tau += np.dot(cs[i] , integral)
    return tau

#calculating optical depth for all frequencies and all thermospheric heights
#this is working but also extremely slow for all wavelengths and heights
#Tau_X0 = []
#for m in range(90, 92):
#    T = []
#    for idx, val in enumerate(I_data["wavelength"]):
#        T.append(Tau(idx,m,0))
#    print(T)
#    Tau_X0.append(T)
optical_depth = []
for m in range(90,600):
    optical_depth.append(Tau(m, 0) )

sns.heatmap(optical_depth)
plt.title("optical depth [m]")
plt.show()

sns.heatmap(optical_depth, norm=mcolors.LogNorm())
plt.title("optical depth [m]")
plt.show()

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


#wavelength variations of cross sections
plt.plot(cs_data.iloc[:,0], cs_data.iloc[:,1], label = "N2")
plt.plot(cs_data.iloc[:,0], cs_data.iloc[:,2], label = "O")
plt.plot(cs_data.iloc[:,0], cs_data.iloc[:,3], label = "O2")
plt.title("absorption cross sections for different species")
plt.xlabel("wavelength [m]")
plt.ylabel("absorption cross section [m^2]")
plt.legend()
plt.show()

#wavelength variations of cross sections
plt.plot(cs_N2, label = "N2")
plt.plot(cs_O, label = "O")
plt.plot(cs_O2, label = "O2")
plt.title("absorption cross sections for different species, interpolated")
plt.xlabel("wavelength [m]")
plt.ylabel("absorption cross section [m^2]")
plt.legend()
plt.show()