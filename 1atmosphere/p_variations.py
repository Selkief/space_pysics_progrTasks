#calculate atmospheric pressure variation with equation for atmospheric pressure
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