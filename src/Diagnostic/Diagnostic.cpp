#include "Diagnostic.h"

Diagnostic::Diagnostic(PicParams &params, SmileiMPI* smpi, vector<Species*>& vecSpecies, vector<Collisions*> &vecCollisions, vector<PSI*>& vecPSI) :
sim_length(params.sim_length),
timestep(params.timestep),
step_dump(params.dump_step),
step_ave(params.ntime_step_avg),
n_species(params.species_param.size()),
n_collision(vecCollisions.size()),
n_psi(vecPSI.size()),
n_dim_field(params.nDim_field),
n_space(params.n_space),
n_space_global(params.n_space_global)

{
	pi_ov_2 = 0.5 * params.const_pi;
	const_e = params.const_e;
	oversize = params.oversize;

	dim.resize( n_dim_field );
    dim_global.resize( n_dim_field );
    for (size_t i=0 ; i<n_dim_field ; i++) 
	{
        dim[i] = n_space[i] + 1 + 2 * oversize[i];
        dim_global[i] = n_space_global[i] + 1;
    }

	double B_magnitude = pow(params.externB[0], 2) + pow(params.externB[1], 2) + pow(params.externB[2], 2);
	B_magnitude = sqrt(B_magnitude);
	sinPhi = params.externB[0] / B_magnitude;
	cosPhi = sqrt(1.0 - sinPhi * sinPhi);
}
