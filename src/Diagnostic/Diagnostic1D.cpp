#include "Diagnostic1D.h"
#include "Species.h"
#include "SmileiMPI_Cart1D.h"
#include "ElectroMagn.h"
#include "Collisions.h"
#include <iomanip>
#include <fstream>

using namespace std;

Diagnostic1D::Diagnostic1D(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields, vector<Species*>& vecSpecies, vector<Collisions*> &vecCollisions, vector<PSI*>& vecPSI) :
Diagnostic(params, smpi, vecSpecies, vecCollisions, vecPSI)
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	index_domain_begin = smpi1D->getCellStartingGlobalIndex(0);

	dx_inv_  = 1.0 / params.cell_length[0];

	n_energy = 100;
	energy_max = 200.0;

	ptclNum1D = new Field1D(EMfields->dimPrim, "ptclNum");

	particle_number.resize(n_species);
	total_particle_energy.resize(n_species);

	particle_flux_left.resize(n_species);
	particle_flux_right.resize(n_species);
	heat_flux_left.resize(n_species);
	heat_flux_right.resize(n_species);

	angle_distribution_left.resize(n_species);
	angle_distribution_right.resize(n_species);
	energy_distribution_left.resize(n_species);
	energy_distribution_right.resize(n_species);

	for(int ispec = 0; ispec < n_species; ispec++)
	{
		angle_distribution_left[ispec].resize(90);
		angle_distribution_right[ispec].resize(90);
	}

	for(int ispec = 0; ispec < n_species; ispec++)
	{
		energy_distribution_left[ispec].resize(n_energy);
		energy_distribution_right[ispec].resize(n_energy);
	}

	psi_rate_left.resize(n_psi);
	psi_rate_right.resize(n_psi);

	radiative_energy_collision.resize(n_collision);
	radiative_energy_collision_global.resize(n_collision);
	for(int i_collision = 0; i_collision < n_collision; i_collision++)
	{
		radiative_energy_collision[i_collision].resize(n_space_global[0]);
		radiative_energy_collision_global[i_collision].resize(n_space_global[0]);
	}
}

Diagnostic1D::~Diagnostic1D()
{
	delete ptclNum1D;
}



void Diagnostic1D::run( SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int itime )
{
	Species *s1;
	Particles *p1;
	double v_square, v_magnitude, energy;
	double mass_ov_2;
	double wlt;			//weight * cell_length / time for calculating flux
	int i_angle, i_energy;
	double flux_temp;
	vector<double> angle_temp;

	angle_temp.resize(90);

	//reset diagnostic parameters to zero
	if( ((itime - 1) % step_dump) == 0 ) 
	{
		for(int ispec = 0; ispec < n_species; ispec++)
		{
			particle_flux_left[ispec] = 0.0;
			particle_flux_right[ispec] = 0.0;
			heat_flux_left[ispec] = 0.0;
			heat_flux_right[ispec] = 0.0;
			for(int i_angle = 0; i_angle < 90; i_angle++)
			{
				angle_distribution_left[ispec][i_angle] = 0.0;
				angle_distribution_right[ispec][i_angle] = 0.0;
			}
			for(int i_energy = 0; i_energy < n_energy; i_energy++)
			{
				energy_distribution_left[ispec][i_energy] = 0.0;
				energy_distribution_right[ispec][i_energy] = 0.0;
			}
		}

		for(int i_psi = 0; i_psi < n_psi; i_psi++)
		{
			psi_rate_left[i_psi] = 0.0;
			psi_rate_right[i_psi] = 0.0;
		}

		for(int i_collision = 0; i_collision < n_collision; i_collision++)
		{
			for(int i_bin = 0; i_bin < radiative_energy_collision[i_collision].size(); i_bin++)
			{
				radiative_energy_collision[i_collision][i_bin] = 0.0;
			}	
		}
	}

	//calculate particle flux, heat flux, angle distribution, energy distribution
	if( (itime % step_dump) > (step_dump - step_ave) || (itime % step_dump) == 0 )
	{
		for(int ispec = 0; ispec < n_species; ispec++)
		{
			s1 = vecSpecies[ispec];
			p1 = &(s1->psi_particles);
			mass_ov_2 = 0.5 * s1->species_param.mass;

			for(int i_particle = 0; i_particle < p1->size(); i_particle++)
			{
				//cout<<"particle number:  "<<p1->position(0,i_particle)<<endl;
				v_square = p1->momentum(0,i_particle) * p1->momentum(0,i_particle) + p1->momentum(1,i_particle) * p1->momentum(1,i_particle) + p1->momentum(2,i_particle) * p1->momentum(2,i_particle);
				v_magnitude = sqrt(v_square);
				energy = mass_ov_2 * v_square;
				i_angle = 90.0 * acos( abs(p1->momentum(0,i_particle)) / v_magnitude ) / pi_ov_2;
				i_energy = n_energy * energy / energy_max;
				if(i_energy >= n_energy)
				{
					i_energy = n_energy + 1;
				}
				if( p1->position(0,i_particle) < 0.0 ) 
				{
					particle_flux_left[ispec]++;
					heat_flux_left[ispec] += energy;
					if (i_angle >= 0 && i_angle < 90)
					{
						angle_distribution_left[ispec][i_angle]++;
					}
					else
					{
						WARNING("i_angle left out of range: i_angle = "<<i_angle);
					}
					energy_distribution_left[ispec][i_energy]++;					

				}
				else if( p1->position(0,i_particle) > sim_length[0] ) 
				{
					particle_flux_right[ispec]++;
					heat_flux_right[ispec] += energy;
					if (i_angle >= 0 && i_angle < 90)
					{
						angle_distribution_right[ispec][i_angle]++;
					}
					else
					{
						WARNING("i_angle right out of range: i_angle = "<<i_angle);
					}
					energy_distribution_right[ispec][i_energy]++;
				}
			}
		}
	}


	//MPI gather diagnostic parameters to master
	if( (itime % step_dump) == 0 ) 
	{
		vector<double> flux_temp;
		vector<double> angle_distribution_temp;
		vector<double> energy_distribution_temp;
		vector<double> psi_rate_temp;

		flux_temp.resize(n_species);
		angle_distribution_temp.resize(90);
		energy_distribution_temp.resize(n_energy);

		for(int ispec = 0; ispec < n_species; ispec++)
		{
			s1 = vecSpecies[ispec];
			wlt = s1->species_param.weight / (dx_inv_ * timestep * step_ave);
			particleFlux[ispec][0] *= wlt;
			particleFlux[ispec][1] *= wlt;
			heatFlux[ispec][0] 	   *= wlt;
			heatFlux[ispec][1]     *= wlt;

			for(int i_direction = 0; i_direction < 2; i_direction++)
			{
				MPI_Allreduce( &particleFlux[ispec][i_direction], &flux_temp, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
				particleFlux[ispec][i_direction] = flux_temp;

				MPI_Allreduce( &heatFlux[ispec][i_direction], &flux_temp, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
				heatFlux[ispec][i_direction] = flux_temp;

				MPI_Allreduce( &angleDist[ispec][i_direction][0], &angle_temp[0], 90, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
				angleDist[ispec][i_direction] = angle_temp;
			}
		}

		s1 = vecSpecies[0];
		wlt = s1->species_param.weight / (dx_inv_ * timestep * step_ave);
		for(int i_collision = 0; i_collision < n_collision; i_collision++)
		{
			smpi->reduce_sum_double(&radiative_energy_collision[i_collision][0], &radiative_energy_collision_global[i_collision][0], radiative_energy_collision[i_collision].size());
			for(int i_bin = 0; i_bin < radiative_energy_collision_global[i_collision].size(); i_bin++)
			{
				radiative_energy_collision_global[i_collision][i_bin] *= wlt;
			}
		}
	}

	//calculate velocity and temperature of each species
	calVT(smpi, vecSpecies, EMfields, itime);

	//calTotalEnergy(smpi, vecSpecies, EMfields, itime);



}


void Diagnostic1D::calVT(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime)
{
	Species *s1;
	Particles *p1;
	int i_temp;
	double xjn,xjmxi;
	double m_ov_3e;
	double m_ov_2;
	double vx, vy, vz;
	double step_ave_inv_;
	double v_square;
	vector<int> number_temp;
	vector<double> energy_temp;

	number_temp.resize(n_species);
	energy_temp.resize(n_species);
	step_ave_inv_ = 1.0 / step_ave;
	for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
	{
		s1 = vecSpecies[iSpec];
		p1 = &(s1->particles);
		Field1D* Vx1D_s = static_cast<Field1D*>(EMfields->Vx_s[iSpec]);
		Field1D* Vy1D_s = static_cast<Field1D*>(EMfields->Vy_s[iSpec]);
		Field1D* Vz1D_s = static_cast<Field1D*>(EMfields->Vz_s[iSpec]);
		Field1D* Vp1D_s = static_cast<Field1D*>(EMfields->Vp_s[iSpec]);

		Field1D* Vx1D_s_avg = static_cast<Field1D*>(EMfields->Vx_s_avg[iSpec]);
		Field1D* Vy1D_s_avg = static_cast<Field1D*>(EMfields->Vy_s_avg[iSpec]);
		Field1D* Vz1D_s_avg = static_cast<Field1D*>(EMfields->Vz_s_avg[iSpec]);
		Field1D* Vp1D_s_avg = static_cast<Field1D*>(EMfields->Vp_s_avg[iSpec]);

		Field1D* T1D_s = static_cast<Field1D*>(EMfields->T_s[iSpec]);
		Field1D* T1D_s_avg = static_cast<Field1D*>(EMfields->T_s_avg[iSpec]);

		if( ((itime - 1) % step_dump) == 0 ) {
			Vx1D_s_avg->put_to(0.0);
			Vy1D_s_avg->put_to(0.0);
			Vz1D_s_avg->put_to(0.0);
			Vp1D_s_avg->put_to(0.0);
			T1D_s_avg ->put_to(0.0);
		}

		//calculate macroscopic velocity (average velocity) and particle number at grid points
		if( (itime % step_dump) > (step_dump - step_ave) || (itime % step_dump) == 0 )
		{
			m_ov_3e = s1->species_param.mass / ( const_e * 3.0 );
			m_ov_2 = s1->species_param.mass / 2.0;
			ptclNum1D->put_to(0.0);
			Vx1D_s->put_to(0.0);
			Vy1D_s->put_to(0.0);
			Vz1D_s->put_to(0.0);
			Vp1D_s->put_to(0.0);
			T1D_s->put_to(0.0);

			//reset particleNumber and kineticEnergy
			particleNumber[iSpec] = 0;
			kineticEnergy[iSpec] = 0.0;
			//get particleNumber
			particleNumber[iSpec] = p1->size();

			for(int i_particle = 0; i_particle < p1->size(); i_particle++)
			{
				//Locate particle on the grid
				xjn    = p1->position(0, i_particle) * dx_inv_;  	//normalized distance to the first node
				i_temp      = floor(xjn);                   		//index of the central node
				xjmxi  = xjn - (double)i_temp;              		//normalized distance to the nearest grid point

				i_temp -= index_domain_begin;
				(*ptclNum1D)(i_temp) 	+= 1.0;
				(*Vx1D_s)(i_temp) 		+= p1->momentum(0, i_particle);
				(*Vy1D_s)(i_temp) 		+= p1->momentum(1, i_particle);
				(*Vz1D_s)(i_temp) 		+= p1->momentum(2, i_particle);
				(*Vp1D_s)(i_temp) 		+= (p1->momentum(0, i_particle) * sinPhi + p1->momentum(1, i_particle) * cosPhi);

				//calculate total kineticEnergy
				v_square = p1->momentum(0, i_particle) * p1->momentum(0, i_particle) + p1->momentum(1, i_particle) * p1->momentum(1, i_particle) + p1->momentum(2, i_particle) * p1->momentum(2, i_particle);
				kineticEnergy[iSpec] += ( m_ov_2 * v_square );
			}
			for(int i = 0; i < ptclNum1D->dims_[0]; i++)
			{
				if( (*ptclNum1D)(i) != 0.0 )
				{
					(*Vx1D_s)(i) /= (*ptclNum1D)(i);
					(*Vy1D_s)(i) /= (*ptclNum1D)(i);
					(*Vz1D_s)(i) /= (*ptclNum1D)(i);
					(*Vp1D_s)(i) /= (*ptclNum1D)(i);
				}
			}

			//calculate temperature
			for(int i_particle = 0; i_particle < p1->size(); i_particle++)
			{
				//Locate particle on the grid
				xjn    = p1->position(0, i_particle) * dx_inv_;  	//normalized distance to the first node
				i_temp      = floor(xjn);                   		//index of the central node
				xjmxi  = xjn - (double)i_temp;              		//normalized distance to the nearest grid point

				i_temp -= index_domain_begin;
				//vx = p1->momentum(0, i_particle);
				//vy = p1->momentum(1, i_particle);
				//vz = p1->momentum(2, i_particle);
				vx = p1->momentum(0, i_particle) - (*Vx1D_s)(i_temp);
				vy = p1->momentum(1, i_particle) - (*Vy1D_s)(i_temp);
				vz = p1->momentum(2, i_particle) - (*Vz1D_s)(i_temp);
				(*T1D_s)(i_temp) 		+= ( vx * vx + vy * vy + vz * vz );
			}
			for(int i = 0; i < ptclNum1D->dims_[0]; i++)
			{
				if( (*ptclNum1D)(i) != 0.0 )
				{
					(*T1D_s)(i) = (*T1D_s)(i) * m_ov_3e / (*ptclNum1D)(i);
				}

			}

			//sum velocity and temperature
			for(int i = 0; i < ptclNum1D->dims_[0]; i++)
			{
				(*Vx1D_s_avg)(i) += (*Vx1D_s)(i);
				(*Vy1D_s_avg)(i) += (*Vy1D_s)(i);
				(*Vz1D_s_avg)(i) += (*Vz1D_s)(i);
				(*Vp1D_s_avg)(i) += (*Vp1D_s)(i);
				(*T1D_s_avg)(i)  += (*T1D_s)(i);
			}
		}


		//Calculate the average parameters and MPI gather
		if( (itime % step_dump) == 0 )
		{
			for(int i = 0; i < ptclNum1D->dims_[0]; i++)
			{
				(*Vx1D_s_avg)(i) *= step_ave_inv_;
				(*Vy1D_s_avg)(i) *= step_ave_inv_;
				(*Vz1D_s_avg)(i) *= step_ave_inv_;
				(*Vp1D_s_avg)(i) *= step_ave_inv_;
				(*T1D_s_avg)(i)  *= step_ave_inv_;
			}

			//another way: firstly gather V, T, ptclNum, then calculate V_global, T_global
			SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
			smpi1D->gatherField2( static_cast<Field1D*>(EMfields->Vx_s_global_avg[iSpec]), Vx1D_s_avg );
			smpi1D->gatherField2( static_cast<Field1D*>(EMfields->Vy_s_global_avg[iSpec]), Vy1D_s_avg );
			smpi1D->gatherField2( static_cast<Field1D*>(EMfields->Vz_s_global_avg[iSpec]), Vz1D_s_avg );
			smpi1D->gatherField2( static_cast<Field1D*>(EMfields->Vp_s_global_avg[iSpec]), Vp1D_s_avg );
			smpi1D->gatherField2( static_cast<Field1D*>(EMfields->T_s_global_avg [iSpec]), T1D_s_avg );

		}

	}
	//gather particleNumber and kineticEnergy to master process
	if( (itime % step_dump) == 0 )
	{
		MPI_Allreduce( &particleNumber[0], &number_temp[0], n_species , MPI_INT,MPI_SUM, MPI_COMM_WORLD);
		particleNumber = number_temp;
		MPI_Allreduce( &kineticEnergy[0], &energy_temp[0], n_species , MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
		kineticEnergy = energy_temp;
	}


}


void Diagnostic1D::calTotalEnergy(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime)
{
		Species *s1;
		Particles *p1;
		double m_ov_2, v_square;
		int i;
		double totalElectricFieldEnergy_temp;
		vector<double> totalParticleEnergy_temp;
		totalParticleEnergy_temp.resize(n_species);

		for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
		{
				s1 = vecSpecies[iSpec];
				p1 = &(s1->particles);

				m_ov_2 = s1->species_param.mass / 2.0;

				//reset particleNumber and kineticEnergy
				totalParticleEnergy[iSpec] = 0.0;

				for(int i_particle = 0; i_particle < p1->size(); i_particle++)
				{
					//calculate total kineticEnergy
					v_square = p1->momentum(0, i_particle) * p1->momentum(0, i_particle) + p1->momentum(1, i_particle) * p1->momentum(1, i_particle) + p1->momentum(2, i_particle) * p1->momentum(2, i_particle);
					totalParticleEnergy[iSpec] += ( m_ov_2 * v_square );
				}

		}

		totalElectricFieldEnergy = 0.0;
		Field1D* Ex1D = static_cast<Field1D*>(EMfields->Ex_);
		for(int i = oversize[0]; i < Ex1D->dims_[0] - oversize[0]; i++)
		{
				totalElectricFieldEnergy += (*Ex1D)(i) * (*Ex1D)(i);
		}

		smpi->reduce_sum_double(&totalParticleEnergy[0], &totalParticleEnergy_temp[0], n_species);
		smpi->reduce_sum_double(&totalElectricFieldEnergy, &totalElectricFieldEnergy_temp, 1);
		totalParticleEnergy = totalParticleEnergy_temp;
		totalElectricFieldEnergy = totalElectricFieldEnergy_temp;

		if(smpi->isMaster())
		{
				ofstream outfile;
				outfile.open("totalEnergy.txt", ios::app);
				outfile<<setw(15)<<itime<<setw(15)<<totalElectricFieldEnergy<<setw(15)<<totalParticleEnergy[0]<<setw(15)<<totalParticleEnergy[1]
							 <<setw(15)<<totalParticleEnergy[0]+totalParticleEnergy[1]<<endl;
				outfile.close();
		}

}
