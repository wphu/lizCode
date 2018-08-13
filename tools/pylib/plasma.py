import math
from scipy import constants as const

Te = 30.0
Ti = 30.0

ne = 2.0e18
ni = 2.0e18

Bmag = 2.0
Bangle = (180.0 - 5.0) * math.pi / 180.0

# mass of H: 1.67262158e-27, D: 2.0 * 1.67262158e-27, T: 3.0 * 1.67262158e-27
me = 9.109382616e-31
mi = 2.0 * 1.67262158e-27
qi = 1

debye_length_e      = math.sqrt(const.epsilon_0 * Te / (ne * const.e))
debye_length_ion    = math.sqrt(const.epsilon_0 * Ti / (ni * qi * const.e))
debye_length        = math.sqrt(const.epsilon_0 * Te *  Ti / ((Te + Ti) * ne * qi * const.e))
period_plasma       = math.sqrt(const.epsilon_0 * me / (ne * const.e * const.e)) * 2.0 * const.pi

ion_sound_speed     = math.sqrt(Te * const.e / mi)
thermal_speed_e     = math.sqrt(Te * const.e / me)
thermal_speed_ion   = math.sqrt(Ti * const.e / mi)

rotation_period_e   = 2.0 * const.pi * me / (const.e * Bmag)
rotation_period_ion = 2.0 * const.pi * mi / (qi * const.e * Bmag)
rotation_radius_e   = me * thermal_speed_e / (const.e * Bmag)
rotation_radius_ion = mi * thermal_speed_ion / (qi * const.e * Bmag)

particle_flux       = ni * thermal_speed_ion * math.sin(Bangle)

print("==========================================")
print("debye_length:        ", debye_length)
print("debye_length_e:      ", debye_length_e)
print("debye_length_ion:    ", debye_length_ion)
print("period_plasma:       ", period_plasma)
print("==========================================")
print(" ")

print("==========================================")
print("rotation_radius_e:    ", rotation_radius_e)
print("rotation_radius_ion:  ", rotation_radius_ion)
print("rotation_period_e:    ", rotation_period_e)
print("rotation_period_ion:  ", rotation_period_ion)
print("==========================================")
print(" ")

print("==========================================")
print("particle_flux:    ", particle_flux)
print("==========================================")
print(" ")


# ===================== iter_langmuir_probe geometry paramters ==================
'''
x1 = 20.0e-3
x2 = 2.0e-3
x3 = 10.0e-3
x4 = 2.0e-3
y1 = 2.0e-3
y2 = 1.0e-3
'''

x1 = 5.0e-3
x2 = 3.0e-3
x3 = 1.5e-3
x4 = 2.0e-3
y1 = 2.0e-3
y2 = 1.5e-3

Bangle = 5.0 * math.pi / 180.0

# the height of intersection between the magnetic and left and right boudnary of the probe
h0 = y1 - x3 * math.tan(Bangle)
h1 = y1 - (x3 + x4) * math.tan(Bangle)

print("==========================================")
print("probe height:    ", y2)
print("probe left  :    ", h0)
print("probe right :    ", h1)
print("==========================================")
print(" ")