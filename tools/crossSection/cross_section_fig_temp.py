import numpy as np
import math
import matplotlib.pyplot as plt


a0 = 5.2917706e-11


filename = "/home/huwanpeng/codes/lizCode/tools/crossSection/C/data/Excitation_C_2s22p23P-2s22p3d3D.dat"
data = np.loadtxt(filename)
data[:,1] = data[:,1] / (3.14 * a0 * a0)
plt.plot(data[:,0], data[:,1], label = "cross section")

#np.savetxt(filename_new, data, fmt='%1.5e')

plt.legend()
plt.xlim((0.0, 500.0))
plt.show()
#plt.savefig("fig/Ionization_H.png")
