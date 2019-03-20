# Ref: 1973 Electron modecule collision ionization in hydrogen and deuterium
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import constants as const

filename_old = "original/Ionization_D_to_D+1.dat"
filename_new = "data/Ionization_D_to_D+1.dat"
data = np.loadtxt(filename_old)
np.savetxt(filename_new, data, fmt='%1.5e')
plt.plot(data[:,0], data[:,1], label = "Ionization_D_to_D+1")

plt.legend()
plt.xlim((0.0, 500.0))
plt.savefig("fig/Ionization_D.png")


E1 = 20.0 * const.e
me = 9.109382616e-31
V1 = math.sqrt(2.0 * E1 / me)

E2 = 20.0 * const.e
me = 9.109382616e-31
V2 = math.sqrt(2.0 * E2 / me)

# ne: electron density, nn: neutral density
ne = 1.0e20
nn = 1.0e20

dt = 1.0e-12

f_linear = interpolate.interp1d(data[:,0], data[:, 1])
cross_section = f_linear(E1 / const.e)
print("cross_section: ", cross_section)

P = 1.0 - math.exp( -V1 * cross_section * nn * dt )
print("Ionization collision probability: ", P)



plt.show()