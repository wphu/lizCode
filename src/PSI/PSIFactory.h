#ifndef PSIFACTORY_H
#define PSIFACTORY_H

#include "PSI1D_SEE.h"
#include "PSI1D_Sputtering.h"
#include "PSI1D_RefDep.h"
#include "PSI1D_Recycling.h"

#include "PSI2D_SEE.h"
#include "PSI2D_Sputtering.h"
#include "PSI2D_RefDep.h"


#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>
using namespace std;

class PSIFactory {

public:

	// Reads the input file and creates the PSI objects accordingly
	static vector<PSI*> create(PicParams& params, InputData &ifile, vector<Species*>& vecSpecies, SmileiMPI* smpi)
	{
	    vector<PSI*> vecPSI;

		string PSI_type;
		string emitKind;
	    vector<string> sg1, sg2, sg3;
	    vector<unsigned int> sgroup1, sgroup2, sgroup3;
	    string psiPos;
		double emitTemp;
		double weight_const;
		double emitOffset;
		double a_FN;
		double b_FN;
		double work_function;
		double emitJ;
		double emitFlux;
		unsigned int nPartEmit;
		double recycling_factor;
		bool is_self_consistent;
		string relSpecies;


	    bool intra, debye_length_required = false;
	    int debug_every;
	    ostringstream mystream;


	    // Loop over each binary PSI group and parse info
	    unsigned int numPSI=ifile.nComponents("PSI");
	    for (unsigned int i_psi = 0; i_psi < numPSI; i_psi++) 
		{
			ifile.extract("PSI_type",PSI_type,"PSI",i_psi);
			if(params.geometry == "1d3v" && PSI_type == "sputtering")
			{
		        ifile.extract("species1",sg1,"PSI",i_psi);
				ifile.extract("species2",sg2,"PSI",i_psi);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
				sgroup2 = params.FindSpecies(sg2);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << i_psi);
				if (sgroup2.size()==0) ERROR("No valid `species2` requested in PSI #" << i_psi);

				// Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",i_psi);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",i_psi);

				// is_self_consistent
				is_self_consistent = false;
				ifile.extract("is_self_consistent",is_self_consistent,"PSI",i_psi);

		        vecPSI.push_back( new PSI1D_Sputtering(params, smpi, vecSpecies, i_psi, sgroup1[0],sgroup2[0],is_self_consistent, psiPos, emitTemp) );

			}

			else if(params.geometry == "1d3v" && PSI_type == "RefDep")
			{
		        ifile.extract("species1",sg1,"PSI",i_psi);
				ifile.extract("species2",sg2,"PSI",i_psi);
				ifile.extract("species3",sg3,"PSI",i_psi);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
				sgroup2 = params.FindSpecies(sg2);
				sgroup3 = params.FindSpecies(sg3);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << i_psi);
				if (sgroup2.size()==0) ERROR("No valid `species2` requested in PSI #" << i_psi);
				if (sgroup3.size()==0) WARNING("No valid `species3` maybe requested in PSI #" << i_psi);

				// Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",i_psi);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",i_psi);

				// is_self_consistent
				is_self_consistent = false;
				ifile.extract("is_self_consistent",is_self_consistent,"PSI",i_psi);

		        vecPSI.push_back( new PSI1D_RefDep(params, smpi, vecSpecies, i_psi, sgroup1[0], sgroup2[0], sgroup3[0], is_self_consistent, psiPos, emitTemp) );
			}

			else if(params.geometry == "1d3v" && PSI_type == "Recycling")
			{
		        ifile.extract("species1",sg1,"PSI",i_psi);
				ifile.extract("species2",sg2,"PSI",i_psi);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
				sgroup2 = params.FindSpecies(sg2);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << i_psi);
				if (sgroup2.size()==0) ERROR("No valid `species1` requested in PSI #" << i_psi);

				// Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",i_psi);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",i_psi);

				recycling_factor = 0.0; // default
				ifile.extract("recycling_factor",recycling_factor,"PSI",i_psi);

				// is_self_consistent
				is_self_consistent = false;
				ifile.extract("is_self_consistent",is_self_consistent,"PSI",i_psi);

		        vecPSI.push_back( new PSI1D_Recycling(params, smpi, i_psi, sgroup1[0], sgroup2[0], is_self_consistent, psiPos, emitTemp, recycling_factor) );
			}

			else if(params.geometry == "1d3v" && PSI_type == "SEE")
			{
		        ifile.extract("species1",sg1,"PSI",i_psi);
				ifile.extract("species2",sg2,"PSI",i_psi);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
				sgroup2 = params.FindSpecies(sg2);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << i_psi);
				if (sgroup2.size()==0) ERROR("No valid `species2` requested in PSI #" << i_psi);

				// Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",i_psi);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",i_psi);

				// secondary electron emission yield
		        double SEEYield = 0.0; // default
		        ifile.extract("SEEYield",SEEYield,"PSI",i_psi);

				// is_self_consistent
				is_self_consistent = false;
				ifile.extract("is_self_consistent",is_self_consistent,"PSI",i_psi);
				
		        vecPSI.push_back( new PSI1D_SEE(params, smpi, i_psi, sgroup1[0],sgroup2[0], is_self_consistent, psiPos, emitTemp, SEEYield) );

			}

			else {
				ERROR("no PSI_type match: "<<PSI_type);
			}
	    }

	    return vecPSI;
	};

};

#endif
