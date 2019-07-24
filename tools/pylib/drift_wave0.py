#ref: On thermionic emission from plasma-facing components in tokamak-relevant conditions
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import constants as const
import copy


font={	'family' : 'serif',
	'weight' : 'bold',
	'size' : 12,
	}

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['lines.linewidth'] = 2.0
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.color'] = "black"

Aeff = 60.0
Wf = 4.55 * const.e

x = np.linspace(0.0, 2.0 * const.pi, 2000)
n = copy.copy(x)
v = copy.copy(x)
nv = copy.copy(x)

n = 1.0 * np.sin(x)
v = 1.0 * np.cos(0.25 * const.pi + x)


sum = 0.0

for i in np.arange(0, x.shape[0]):
    sum = sum + n[i] * v[i]
    nv[i] = n[i] * v[i]

print(sum)



fig = plt.figure()

linestyles = ['-', '--', '-.', ':']
markers = ['x', '^', 's', '*']
labels = ['n', 'v', 'nv']


ax = fig.add_subplot(1,1,1)
ax.plot(x, n, label = labels[0], linestyle = linestyles[0])
ax.plot(x, v, label = labels[1], linestyle = linestyles[1])
ax.plot(x, nv, label = labels[2], linestyle = linestyles[2])

major_ticks = np.arange(0, 2.1 * const.pi, 0.5 * const.pi)
ax.set_xticks(major_ticks)
ax.set_xticklabels(['0', '0.5pi', 'pi', '1.5pi', '2pi'])

major_ticks = np.arange(-1.0, 1.1, 0.5)
ax.set_yticks(major_ticks)


ax.legend(loc = 1, framealpha=1)

ax.set_xlabel(r'$x$')
ax.grid(True)


plt.show()
#fig.savefig("fig.png")
