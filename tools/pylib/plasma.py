import math
from scipy import constants as const
import numpy as np

Te = 50.0
Ti = 50.0

ne = 1.0e18
ni = 1.0e18

B_mag = 2.0
B_angle_alpha = 10.0 * math.pi / 180.0
B_angle_beta =  90.0 * math.pi / 180.0
B_x =  B_mag * math.sin(B_angle_alpha) * math.cos(B_angle_beta)
B_y = -B_mag * math.cos(B_angle_alpha)
B_z = -B_mag * math.sin(B_angle_alpha) * math.sin(B_angle_beta)

dx = 5.0e-5


# mass of H: 1.67262158e-27, D: 2.0 * 1.67262158e-27, T: 3.0 * 1.67262158e-27
me = 9.109382616e-31
mi = 2.0 * 1.67262158e-27
mr = me * mi / ( me + mi )
qi = 1

debye_length_e      = math.sqrt(const.epsilon_0 * Te / (ne * const.e))
debye_length_ion    = math.sqrt(const.epsilon_0 * Ti / (ni * qi * const.e))
debye_length        = math.sqrt(const.epsilon_0 * Te *  Ti / ((Te + Ti) * ne * qi * const.e))
period_plasma       = math.sqrt(const.epsilon_0 * me / (ne * const.e * const.e)) * 2.0 * const.pi

ion_sound_speed     = math.sqrt(Te * const.e / mi)
thermal_speed_e     = math.sqrt(Te * const.e / me)
thermal_speed_ion   = math.sqrt(Ti * const.e / mi)
ve                  = thermal_speed_e
vi                  = thermal_speed_ion

ion_sound_speed_x =  ion_sound_speed * math.sin(B_angle_alpha) * math.cos(B_angle_beta)
ion_sound_speed_y = -ion_sound_speed * math.cos(B_angle_alpha)
ion_sound_speed_z = -ion_sound_speed * math.sin(B_angle_alpha) * math.sin(B_angle_beta)



Omega_e             = const.e * B_mag / me
Omega_i             = const.e * B_mag / mi
rotation_period_e   = 2.0 * const.pi * me / (const.e * B_mag)
rotation_period_ion = 2.0 * const.pi * mi / (qi * const.e * B_mag)
rotation_radius_e   = me * thermal_speed_e / (const.e * B_mag)
rotation_radius_ion = mi * thermal_speed_ion / (qi * const.e * B_mag)

particle_flux          = ni * ion_sound_speed * math.sin(B_angle_alpha)
ion_saturation_current = ni * ion_sound_speed * const.e
macro_particle_weight  = ne * dx * dx * dx / 100.0

#================== diagnostic ne Te Is ========================
#saturated ion current density A/m^2

Js = np.array([6.364, 5.271, 20.81, 15.47])
Js = Js  * 1.0e4
#divertor temperature and density: tet, net
tet = np.array([1.746, 3.238, 10.93, 5.117]) * 5.0
miH = mi  / 2.0
net = Js * 1.0 / (const.e * np.sqrt(const.e * tet / miH))

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
print("ion_sound_speed:        ", ion_sound_speed)
print("thermal_speed_ion:      ", thermal_speed_ion)
print("particle_flux:          ", particle_flux)
print("ion_saturation_current: ", ion_saturation_current)
print("macro_particle_weight:  ", macro_particle_weight)
print("==========================================")
print(" ")

print("================== parameters for 3d simulation ========================")
print("plasma temperature in K:  ", Te * const.e / const.k)
print("magnetic field strenth:   ", B_mag, "(", B_x, B_y, B_z, ")")
print("ion_sound_speed:          ", ion_sound_speed, "(", ion_sound_speed_x, ion_sound_speed_y, ion_sound_speed_z, ")")
print("========================================================================")
print(" ")


print("================== diagnostic ne Te Is ========================")
print("divertor density:  ", net)
print("========================================================================")
print(" ")


# ===================== classical plasma diffusion coeffiecience in magnetic ==================
#coulomb logarithm from  ref: https://farside.ph.utexas.edu/teaching/plasma/Plasma/node39.html
def fun_coulomb_log(ne, Te):
    if Te < 10.0:
        return 23.0 - math.log( math.pow(ne * 1.0e-6, 0.5) * math.pow(Te, -1.5) )
    else:
        return 24.0 - math.log( math.pow(ne * 1.0e-6, 0.5) * math.pow(Te, -1.0) )

coulomb_log0 = math.log( math.sqrt( const.epsilon_0 * Te * const.e / ( ne * const.e * const.e ) ) 
              / ( qi * const.e * qi * const.e / ( 4.0 * const.pi * const.epsilon_0 * Te * const.e) ) )

coulomb_log0_fun = fun_coulomb_log(ne, Te)

# one particle to one thermal species, niu_k_ee: energy loss frequency for electron-electron collision
# niu_p_ee: momentum loss frequency for electron-electron collision
coef0 = math.pow(const.e, 4) / ( math.pow(4.0 * const.pi * const.epsilon_0, 2) )
niu_k_ee = ne * coef0 * ( 8.0 * const.pi / ( me * me * ve * ve * ve ) ) * coulomb_log0
niu_k_ei = ni * math.pow(qi, 2) * coef0 * ( 8.0 * const.pi / ( me * mi * ve * ve * ve ) ) * coulomb_log0
niu_k_ii = ni * math.pow(qi, 2) * coef0 * ( 8.0 * const.pi / ( mi * mi * vi * vi * vi ) ) * coulomb_log0
niu_k_ie = ne * coef0 * ( 8.0 * const.pi / ( mi * me * vi * vi * vi ) ) * coulomb_log0

niu_p_ee = niu_k_ee
niu_p_ei = niu_k_ei * (mi / (2.0 * me))
niu_p_ii = niu_k_ii
niu_p_ie = niu_k_ie * 0.5

# thermal distribution collisions: one thermal species to another thermal species
# below are all momentum loss frequency
niu_ei = math.sqrt(2.0 / const.pi) * (1.0 / 3.0) * ni * math.pow(qi, 2) * coef0 * ( 4.0 * const.pi / ( math.pow(me, 0.5) * math.pow(Te * const.e, 1.5) ) ) * coulomb_log0
niu_ee = niu_ei / math.pow(2.0, 0.5)
niu_ie = niu_ei * ne * me / (ni * mi)
niu_ii = niu_ei * math.pow(me / mi, 0.5)

# particle diffusion coefficient: De_par for electron in parallel direction, De_per for electron in perpendicular direction
De_per = ve * ve / niu_ei
Di_per = vi * vi / niu_ii
De_par = ve * ve * niu_ei / (Omega_e * Omega_e)
Di_par = vi * vi * niu_ie / (Omega_i * Omega_i)

print("==========================================")
print("coulomb logarithm:                  ", coulomb_log0)
print("coulomb logarithm from function:    ", coulomb_log0_fun)
print("De_per, Di_per:                     ", De_per, Di_per)
print("De_par, Di_par:                     ", De_par, Di_par)
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

B_angle_alpha = 5.0 * math.pi / 180.0

# the height of intersection between the magnetic and left and right boudnary of the probe
h0 = y1 - x3 * math.tan(B_angle_alpha)
h1 = y1 - (x3 + x4) * math.tan(B_angle_alpha)

D_per = 1.0
y_diffusion = math.pow(D_per * x1 / vi, 0.5)
lamada_sol  = math.pow(D_per * 50.0 / vi, 0.5)


print("==========================================")
print("probe height:    ", y2)
print("probe left  :    ", h0)
print("probe right :    ", h1)
print("vi          :    ", vi)
print("y_diffusion :    ", y_diffusion)
print("lamada_sol  :    ", lamada_sol)
print("==========================================")
print(" ")
