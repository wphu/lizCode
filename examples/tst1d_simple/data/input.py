import math

method = 'explicit'

l0 = 0.5e-5
Lsim = [500.*l0]

t0 = 0.5e-12
ns = int(1.0e-9 / t0)
Tsim = 5 * ns
number_output = 5

n_procs = 24

B = 0.0
Bangle = 90.0

plasma_temperature = 20

plasma_density = 1.0e18

particle_number_per_cell = 100
particle_number_per_cell_for_weight = 100



#========================================================
dump_step = int( Tsim / number_output )
timesteps_restore = dump_step
ntime_step_avg = dump_step

ion_step = 1


dim = '1d3v'

number_of_procs = [n_procs]


bc_em_type_x = ['Dirichlet', 'Dirichlet']
#bc_em_type_x = ['Neumann', 'Dirichlet']

bc_em_value_x = [0.0, 0.0]


angle = Bangle * math.pi / 180.0
Bx = B * math.sin(angle)
By = B * math.cos(angle)
Bz = 0.0
externB = [Bx, By, Bz]

ion_sound_velocity = 0.0   #math.sqrt( (plasma_temperature * 1.6021766208e-19) / (2.0 * 1.67262158e-27) )
vx = ion_sound_velocity * math.sin(angle)
vy = ion_sound_velocity * math.cos(angle)
vz = 0.0




random_seed = 0

interpolation_order = 1
projection_order = 1

cell_length = [l0]
sim_length  = Lsim

timestep = t0
n_time = Tsim





Species(
	species_type = 'e',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell',
	ionization_model = 'none',
	#Pusher_type = 'GC0',
	n_part_per_cell = particle_number_per_cell,
	n_part_per_cell_for_weight = particle_number_per_cell_for_weight,
	c_part_max = 1.0,
	mass = 9.109382616e-31,
	charge = -1.6021766208e-19,
	nb_density = plasma_density,
	temperature = [plasma_temperature],
	time_frozen = 0.,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',
)


Species(
	species_type = 'D1',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell',
	ionization_model = 'none',
	timestep_zoom = ion_step,
	n_part_per_cell = particle_number_per_cell,
	n_part_per_cell_for_weight = particle_number_per_cell_for_weight,
	c_part_max = 1.0,
	mass = 2.0 * 1.67262158e-27,
	charge = 1.6021766208e-19,
	nb_density = plasma_density,
	temperature = [plasma_temperature],
	time_frozen = 0.0,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',

	diameter = 2.751E-10,
	ref_temperature = 273.,
	visc_temp_index = 0.75,
	vss_scat_inv = 1.
)
