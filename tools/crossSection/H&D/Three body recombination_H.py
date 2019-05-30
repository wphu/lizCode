# Ref: Kinetic Modelling of the Plasma Recombination, Contrib. Plasma Phys. 56, No. 6-8, 698 – 704 (2016) / DOI 10.1002/ctpp.201611004
# The cross section is for Hydrogen
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy import constants as const
import copy


# formular (13)
def cross_section_sub(E1, E2, n):
    A  = 0.3440
    a  = [-0.014353, 0.75206, -0.29548, 0.056884]
    Ry = 13.6 * const.e
    En = Ry / (n*n)
    e1 = 1.0 + E1 / En
    e2 = 1.0 + E2 / En
    e12 = e1 - e2

    cs = 1.0 / (e1 * e1) + 1.0 / (e12 * e12) - 1.0 / (e1 * e12)
    # this formula may be wrong in the ref
    #cs = 1.0 / (e1 * e1) + 1.0 / (e12 * e12) - 1.0 / (e1 * e12) - 1.0 / (e1 * e2)
    sum0 = 0.0
    for i in np.arange(0, 4):
        sum0 += (a[i] / math.pow(e2, i))
    cs += (sum0 * math.log(e1) / (n * math.pow(e2, 3)))
    cs *= (A * math.pow(n, 4) / (En * e1))
    return cs


# cross section for each energy state, formular (11)
def cross_section_each(E1, E2, n):
    Ry = 13.6 * const.e
    En = Ry / (n*n)
    me = 9.109382616e-31
    # a0: Bohr radius
    a0 = 5.2917721067e-11
    E1_prime = E1 + E2 + En
    cs = 0.0
    if E1 >= E2:
        cs = 0.5 * cross_section_sub(E1_prime, E2, n)
    else:
        cs = 0.5 * cross_section_sub(E1_prime, E1, n)
    cs *= 2.0

    coefficient = 4.0 * math.pow(Ry, 1.5) * math.sqrt(2.0 * me) * n * n * math.pi * math.pi * math.pow(a0, 3) * E1_prime / (E1 *E2)
    cs *= coefficient
    return cs



# formular (14), total cross section for three body recombination
def cross_section_TBR(E1, E2, nmax):
    cs = 0.0
    for i in np.arange(1, nmax+1):
        cs += cross_section_each(E1, E2, i)
    return cs

dE1 = 2.0e-2
nE1 = 50

dE2 = 2.0e-2
nE2 = 50

nmax = 20

x, y = np.mgrid[slice(0.0,dE1*nE1,dE1), slice(0.0,dE2*nE2,dE2)]
cross_section = copy.deepcopy(x)

for i in np.arange(0, nE1):
    for j in np.arange(0, nE2):
        cross_section[i,j] = cross_section_TBR(dE1*(i+1) * const.e, dE2*(j+1) * const.e, nmax)



fig = plt.figure()
# 2d figure

ax = fig.add_subplot(1,1,1)
ax.plot(x[:,40], cross_section[40,:])
#plt.yscale('log')
#plt.xscale('log')


# 3d figure
'''
ax = Axes3D(fig)
ax.plot_surface(x, y, cross_section)
ax.set_xlim(0.0, nE1*dE1)
ax.set_ylim(0.0, nE2*dE2)
'''

'''
E1 = 0.1 * const.e
me = 9.109382616e-31
V1 = math.sqrt(2.0 * E1 / me)

E2 = 0.1 * const.e
me = 9.109382616e-31
V2 = math.sqrt(2.0 * E2 / me)

ne = 1.0e20
ni = 1.0e20

dt = 1.0e-12

nmax = 20

print("cross_section_TBR: ",cross_section_TBR(E1, E2, nmax))
print( - math.sqrt(V1 * V2) * cross_section_TBR(E1, E2, nmax) * math.sqrt(ne * ni) * dt )

# The formular of probability of TB recombination collision below formular (14), the dimension in exp() is not 1.
#P = 1.0 - math.exp( -V1 * V2 * cross_section_TBR(E1, E2, nmax) * ne * ni * dt )

# So the formular is changed as below:
P = 1.0 - math.exp( - math.sqrt(V1 * V2) * cross_section_TBR(E1, E2, nmax) * math.sqrt(ne * ni) * dt )
print("recombination collision probability: ", P)
'''


plt.show()
fig.savefig("fig/Three_body_recombination_H.png")