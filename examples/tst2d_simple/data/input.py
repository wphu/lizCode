import math

method = 'explicit'

l0 = 1.0e-5
nx = 100
ny = 200
Lsim = [nx*l0,ny*l0]

t0 = 1.0e-12
ns = int(1.0e-9 / t0)
Tsim = 5 # * ns
number_output = 5

number_of_procs = [4, 6]

B = 2.0
Bangle = 5.0

plasma_temperature = 20

plasma_density = 1.0e18

particle_number_per_cell = 0
particle_number_per_cell_for_weight = 100


#===============================================================
dump_step = int( Tsim / number_output )
timesteps_restore = dump_step
ntime_step_avg = dump_step

ion_step = 1

is_calVDF = 0

dim = '2d3v'

interpolation_order = 1

cell_length = [l0,l0]
sim_length  = Lsim

timestep = t0
n_time = Tsim

bc_em_type_x = ['periodic']
bc_em_type_y = ['silver-muller']
bc_em_value_x = [0.0, 0.0]


Bx = -B * math.cos(angle)
By = -B * math.sin(angle)
Bz = 0.0
externB = [Bx, By, Bz]

ion_sound_velocity = math.sqrt( (plasma_temperature * 1.6021766208e-19) / (2.0 * 1.67262158e-27) )
vx = -ion_sound_velocity * math.cos(angle)
vy = -ion_sound_velocity * math.sin(angle)
vz = 0.0

random_seed = 0


Grid(
	gridType = "from_file",
)



Species(
	species_type = 'e',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell',
	ionization_model = 'none',
	n_part_per_cell = particle_number_per_cell,
	n_part_per_cell_for_weight = particle_number_per_cell_for_weight,
	c_part_max = 1.0,
	mass = 9.109382616e-31,
	charge = -1.6021766208e-19,
	nb_density = plasma_density,
	temperature = [plasma_temperature],
	mean_velocity = [vx, vy, vz],
	time_frozen = 0.,
	bc_part_type_west  = 'periodic',
	bc_part_type_east  = 'periodic',
	bc_part_type_south = 'supp',
	bc_part_type_north = 'supp'
)


Species(
	species_type = 'D1',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell',
	ionization_model = 'none',
	n_part_per_cell = particle_number_per_cell,
	n_part_per_cell_for_weight = particle_number_per_cell_for_weight,
	c_part_max = 1.0,
	mass = 2.0 * 1.67262158e-27,
	charge = 1.6021766208e-19,
	nb_density = plasma_density,
	temperature = [plasma_temperature],
	mean_velocity = [vx, vy, vz],
	time_frozen = 0.0,
	bc_part_type_west  = 'periodic',
	bc_part_type_east  = 'periodic',
	bc_part_type_south = 'supp',
	bc_part_type_north = 'supp'
)
