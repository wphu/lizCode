# Ref: Kinetic Modelling of the Plasma Recombination, Contrib. Plasma Phys. 56, No. 6-8, 698 – 704 (2016) / DOI 10.1002/ctpp.201611004
# The cross section is for Hydrogen
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy import constants as const

# cross section for each energy state, formular (11)
def cross_section_each(E1, E2, n):
    Ry = 13.6
    En = Ry / (n*n)
    me = 9.109382616e-31
    # a0: Bohr radius
    a0 = 5.2917721067e-11
    E1_prime = E1 + E2 + En
    if E1 >= E2:
        cs = 0.5 * cross_section_sub(E1_prime, E2, n)
    else:
        cs = 0.5 * cross_section_sub(E1_prime, E1, n)
    coefficient = 4.0 * math.pow(Ry, 1.5) * math.sqrt(2.0 * me) * n * n * math.pi * math.pi * math.pow(a0, 3) * E1_prime / (E1 *E2)
    cs *= coefficient
    return cs

# formular (13)
def cross_section_sub(E1, E2, n):
    A  = 0.3440
    a  = [-0.014353, 0.75206, -0.29548, 0.056884]
    Ry = 13.6
    En = Ry / (n*n)
    e1 = 1.0 + E1 / En
    e2 = 1.0 + E2 / En
    e12 = e1 - e2

    cs = 1.0 / (e1 * e1) + 1.0 / (e12 * e12) - 1.0 / (e1 * e2) - 1.0 / (e1 * e12)
    sum0 = 0.0
    for i in np.arange(0, 4):
        sum0 += (a[i] / math.pow(e2, i))
    cs += (sum0 * math.log(e1) / (n * math.pow(e2, 3)))
    cs *= (A * math.pow(n, 4) / (En * e1))
    return cs


# total cross section for three body recombination
def cross_section_TBR(E1, E2, nmax):
    cs = 0.0
    for i in np.arange(1, nmax+1):
        cs += cross_section_each(E1, E2, i)
    return cs

dE1 = 1.0e-3
nE1 = 1000

dE2 = 1.0e-2
nE2 = 100

E2 = 0.01
me = 9.109382616e-31
V2 = math.sqrt(2.0 * E2 * const.e / me)

nmax = 20
ne = 1.0e20
x, y = np.mgrid[slice(0.0,dE1*nE1,dE1), slice(0.0,dE2*nE2,dE2)]
print(x.max(), y.max())

cross_section = np.zeros(nE1)
for i in np.arange(0, nE1):
        cross_section[i] = V2 * ne * cross_section_each(dE1*(i+1), E2, nmax)
        #cross_section[i] = cross_section_TBR(dE1*(i+1), E2, nmax)

print(cross_section[1], cross_section[10], cross_section[500])
fig = plt.figure()

ax = fig.add_subplot(1,1,1)
ax.plot(cross_section)

'''
ax = Axes3D(fig)
ax.plot_surface(x, y, cross_section)
ax.set_xlim(0.0, nE1*dE1)
ax.set_ylim(0.0, nE2*dE2)
'''

plt.show()