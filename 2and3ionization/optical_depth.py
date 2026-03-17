##task 1
##Make functions that calculate the optical depth as a function of altitude and wavelength
# for vertical incidence, and for variable zenith-angle of the incident light

#make plots with altitude variations of thermospheric densities
#make plots with wavelength variations of cross sections

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd

#read in all necessary data
#heightprofiles of number densities [cm^-3]
n_data = pd.read_csv("C:/Users/skief/Documents/UiT/semester8/space physics/progrTasks/1atmosphere/MSIS.dat",sep=r"\s+", skiprows=16)
#Irradiance [W/m^2/nm] for different wavelengths [nm]
I_data = pd.read_csv("2and3ionization/fism_daily_hr19990216.dat")
I_data.rename(columns={"wavelength (nm)": "wavelength", "irradiance (W/m^2/nm)" : "irradiance"}, inplace=True)
#crosssections [m^2] for different species for different wavelengths [m]
cs_data = pd.read_csv("2and3ionization/phot_abs.dat",sep=r"\s+", skiprows=8)

#interpolate cross sections to irradiance-wavelengths
cs_N2 = np.interp( I_data["wavelength"]*1e-9, cs_data.iloc[:,0], cs_data.iloc[:,1] )
cs_O = np.interp( I_data["wavelength"]*1e-9, cs_data.iloc[:,0], cs_data.iloc[:,2] )
cs_O2 = np.interp( I_data["wavelength"]*1e-9, cs_data.iloc[:,0], cs_data.iloc[:,3] )
cs = [cs_O, cs_N2, cs_O2] #ordered to fit the order in the n_data file [m^2]

#number densities all species combined for heightprofiles [cm^-3]
n_all = n_data.iloc[:,1] + n_data.iloc[:,2] + n_data.iloc[:,3]

#relevant angles for Xi
Xi = [0, 15, 30, 45, 60] #[degrees]

################define functions###############
#optical depth for different angles of Xi, returns a list over all wavelength present in input data
#this should work for angles up to 60 degrees
def Tau(n_profiles, altitude, angle):
    tau = 0
    for i, cs_species in enumerate(cs):
        integral = 0
        for z in range( altitude, len(n_data.iloc[:,0]) ):
            integral += n_profiles.iloc[z,i+1] * 1e6 * 1000
        tau += cs_species * integral 
    return tau /np.cos( angle * (np.pi/180) )

#Irradiance !!what does it mean "/nm" in file?! need to convert units?!
def Irradiance(opt_depth, Irradiance):
    return np.array(Irradiance) * np.exp(-np.array(opt_depth))

##convert Watts/m^2 into photons/s/m^2 (Plancks law E= hc/wavelength)
def watts2photons(E_total, wavelengths):
    h = 6.62607015e-34 #plancks constant
    c = 2.99792458e8 #speed of light [m/s]
    E_i = h*c/np.array(wavelengths)
    nr_photons = E_total/E_i
    return nr_photons

###############calculate all the stuff#################
#calculate optical depth for heights in the thermosphere (90-600km)
#and for different angles of incidence

optical_depth = [] #in [m]
irradiance = []
irradiance_ph = []
for idx,X in enumerate(Xi):
    T = []
    for m in range(90,600):
        T.append( Tau(n_data, m, X) )
    I = Irradiance(T,I_data["irradiance"])
    I_ph = watts2photons(I, I_data["wavelength"]*1e-9)
    optical_depth.append(T)
    irradiance.append(I)
    irradiance_ph.append(I_ph)



#make plots, only if file is executed directly
if __name__ == "__main__":
    
    xcoord = I_data["wavelength"]
    ycoord = n_data.iloc[91:,0]
    for key,ele in enumerate(optical_depth):
        fig, axs = plt.subplots(2,1)
        plt.suptitle(f"solar zenith angle {Xi[key]} degrees")
        plot1 = axs[0].pcolormesh(xcoord, ycoord, ele, norm=mcolors.LogNorm(vmin = 1e-2, vmax= 1e3), cmap="plasma")
        plot2 = axs[1].pcolormesh(xcoord, ycoord, irradiance[key], norm=mcolors.LogNorm(vmin = 1e-6, vmax = 1e-2), cmap="plasma")
        cbar1 = fig.colorbar(plot1)
        cbar2 = fig.colorbar(plot2)
        cbar1.set_label("optical depth")
        cbar2.set_label("Irradiance [W $m^{-2} s^{-1} nm^{-1}$]")
        axs[0].set_xlabel("wavelength [nm]")
        axs[0].set_ylabel("height [km]")
        axs[1].set_xlabel("wavelength [nm]")
        axs[1].set_ylabel("height [km]")
        plt.tight_layout()
        plt.show()

    #number densities
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
    plt.plot(np.array(cs_data.iloc[:,0])*1e9, cs_data.iloc[:,1], label = "N2")
    plt.plot(np.array(cs_data.iloc[:,0])*1e9, cs_data.iloc[:,2], label = "O")
    plt.plot(np.array(cs_data.iloc[:,0])*1e9, cs_data.iloc[:,3], label = "O2")

    #wavelength variations of cross sections
    plt.plot(I_data["wavelength"],cs_N2, label = "N2 interpolated")
    plt.plot(I_data["wavelength"],cs_O, label = "O interpolated")
    plt.plot(I_data["wavelength"], cs_O2, label = "O2 interpolated")
    plt.title("absorption cross sections for different species")
    plt.xlabel("wavelength [nm]")
    plt.ylabel("absorption cross section [m^2]")
    plt.legend()
    plt.show()