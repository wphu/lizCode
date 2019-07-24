#ref: On thermionic emission from plasma-facing components in tokamak-relevant conditions
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import constants as const
import copy

Aeff = 60.0
Wf = 4.55 * const.e

Ts = np.linspace(300.0, 3695.0, 300)
#Ts = np.linspace(300.0, 2740.0, 300)
j = copy.copy(Ts)
j = 1.0e4 * Aeff * Ts * Ts * np.exp(- Wf / (const.Boltzmann * Ts)) / 1.0e3

Ts_characteristic = np.linspace(1.0, 10.0, 6)
Ts_characteristic[0] = 300.0
Ts_characteristic[1] = 500.0
Ts_characteristic[2] = 2600.0
Ts_characteristic[3] = 2900.0
Ts_characteristic[4] = 3400.0
Ts_characteristic[5] = 3695.0
j_characteristic = copy.copy(Ts_characteristic)
j_characteristic = 1.0e4 * Aeff * Ts_characteristic * Ts_characteristic * np.exp(- Wf / (const.Boltzmann * Ts_characteristic)) / 1.0e3
j_characteristic = np.around(j_characteristic)
print("Temperature:        ", Ts_characteristic)
print("thermionic current: ", j_characteristic)
print("work function temperature: ", Wf / const.Boltzmann)


fig = plt.figure()

ax = fig.add_subplot(1,1,1)
ax.plot(Ts, j)

ax.set_xlim(300.0, 3695.0)
#ax.set_ylim()

ax.set_xlabel(r'$T_s\ (K)$')
ax.set_ylabel(r'$Thermionic\ current\ of\ tungsten\ (kA/m^2)$')


plt.show()
fig.savefig("fig/thermionic_current_tungsten.png")