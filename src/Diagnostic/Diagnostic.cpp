#include "Diagnostic.h"

Diagnostic::Diagnostic(PicParams &params) :
n_species(params.species_param.size()),
sim_length(params.sim_length),
dump_step(params.dump_step),
avg_step(params.ntime_step_avg),
timestep(params.timestep),
nDim_field(params.nDim_field),
n_space(params.n_space),
n_space_global(params.n_space_global)
{
	PI_ov_2 = 0.5 * params.const_pi;
	const_e = params.const_e;

	dim_.resize( nDim_field );
    dim_global.resize( nDim_field );

	oversize = params.oversize;

    for (size_t i=0 ; i<nDim_field ; i++) {
        dim_[i] = n_space[i] + 1 + 2 * oversize[i];
        dim_global[i] = n_space_global[i] + 1;
    }

	double B_magnitude = pow(params.externB[0], 2) + pow(params.externB[1], 2) + pow(params.externB[2], 2);
	B_magnitude = sqrt(B_magnitude);
	sinPhi = params.externB[0] / B_magnitude;
	cosPhi = sqrt(1.0 - sinPhi * sinPhi);
}
