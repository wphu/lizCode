import math
import random
from scipy import constants as const
import numpy as np

class Particle:
    def __init__(self, name, mass, v):
        self.name = name
        self.mass = mass
        self.v    = v

    def print_velocity(self):
        v_magnitude = math.sqrt(self.v[0] * self.v[0] + self.v[1] * self.v[1] + self.v[2] * self.v[2])
        print(self.name, " velocity: ", v_magnitude, self.v[0], self.v[1], self.v[2])
    def energy(self):
        self.v_squre = self.v[0] * self.v[0] + self.v[1] * self.v[1] + self.v[2] * self.v[2]
        return 0.5 * self.mass * self.v_squre

def collide(p1, p2):
    VRCP = [0.0, 0.0, 0.0]
    VCCM = [0.0, 0.0, 0.0]
    VRC  = [0.0, 0.0, 0.0]
    RML  = 0.5
    RMM  = 0.5

    VRR = 0.0
    for i in np.arange(0,3):
        VRC[i] = p1.v[i] - p2.v[i]
        VRR += (VRC[i] * VRC[i])
    VR = math.sqrt(VRR)

    for i in np.arange(0,3):
        VCCM[i] = RML * p1.v[i] + RMM * p2.v[i]


    RF = random.random()
    B = 2.0 * RF - 1.0
    A = math.sqrt(1.0 - B * B)
    VRCP[0] = B * VR

    RF = random.random()
    C = 2.0 * math.pi * RF
    VRCP[1] = A * math.cos(C) * VR
    VRCP[2] = A * math.sin(C) * VR

    for i in np.arange(0, 3):
        p1.v[i] = VCCM[i] + VRCP[i] * RMM
        p2.v[i] = VCCM[i] - VRCP[i] * RML

random.seed()

mass = 2.0 * 1.67262158e-27
T    = 50
v    = math.sqrt(T * const.e / mass)

p1 = Particle("p1", mass, [0.2 * v, 0.0, 0.0])
p2 = Particle("p2", mass, [-v, 0.0, 0.0])

p1.print_velocity()
p2.print_velocity()
print("total energy: ", p1.energy() + p2.energy())
print("total momentum: ", p1.v[0] + p2.v[0], p1.v[1] + p2.v[1], p1.v[2] + p2.v[2])

collide(p1, p2)

p1.print_velocity()
p2.print_velocity()

print("total energy: ", p1.energy() + p2.energy())
print("total momentum: ", p1.v[0] + p2.v[0], p1.v[1] + p2.v[1], p1.v[2] + p2.v[2])