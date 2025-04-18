/**
 * @file OmegaCV.cpp
 *
 * @brief Implementation of OmegaCV.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "Mixture.h"
#include "TransferModel.h"

// Test headers
// #include "RateLaws.h"
// #include "RateManager.h"
// #include "Reaction.h"
// #include "Kinetics.h"
#include "Constants.h"

// using namespace Mutation::Thermodynamics;
// using namespace Mutation;
// using namespace Mutation::Thermodynamics;
using namespace Mutation::Kinetics;


namespace Mutation {
    namespace Transfer {

/**
* Represents a coupling between chemistry and vibrational energy modes.
*/
class OmegaCV : public TransferModel
{
/**
* The necessary functions for the computation of the  source terms for
* Vibration-Chemistry-Vibration coupling are implemented.
*
* The dissociation of molecules is a procedure which removes energy from the vibrational mode
* of them and re-introduces it as chemical energy of the produced atoms and molecules. In order to accurately determine
* the exact energy lost from the vibrational mode in a dissociation reaction one must use state-to-state models in order to
* account for the fact that not all the vibrational states have the exact same probability to dissociate. Generally, the higher
* vibrational levels have a greater probability to dissociate compared to the lower ones because the latter ones should ladder climb
* before dissociating. Nonetheless, the state-to-state models are not easy to use in a CFD code due to the increase of the computational
* cost. Therefore, one should rely on accurate Boltzmann equilibrium models in order to take into account for this mode of energy transfer.
*
* Taking into account the above, the models for the Vibration-Chemistry-Vibration energy transfer can be generaly placed in to two different categories,
* the non-preferential and the preferential models. In the first category, all the energy states have the same probability to dissociate compared to the
* preferential ones which assume higher probability of dissociation for the higher energy levels. Currently, only the simplest non-preferential model has been
* implemented, but the code gives the opportunity to be extended for more complicated Vibration-Chemistry-Vibration energy transfer models.
*
*/

public:
	OmegaCV(Mutation::Mixture& mix)
		: TransferModel(mix)
	{
		m_ns = mixture().nSpecies();
		mp_wrk1 = new double [m_ns];
		mp_wrk2 = new double [m_ns];
	};

	~OmegaCV()
	{
		delete [] mp_wrk1;
		delete [] mp_wrk2;
	};
/**
 * Computes the source terms of the Vibration-Chemistry energy transfer in \f$ [J/(m^3\cdot s)] \f$
 *
 * \f[ \Omega^{CV}_{mi} = \sum c_1 e^V_{mi} \dot{\omega}_i\f]
 *
 */
	double source()
	{
		// static int i_transfer_model = 0;
		static int i_transfer_model = (dynamic_cast<const MMT*>(mixture().reactions()[1].rateLaw()) != NULL);
		switch (i_transfer_model){
		   case 0:
			  return compute_source_Candler();
		   case 1:
			  return compute_source_MMT();
		  break;
		   default:
			  std::cerr << "The selected Chemistry-Vibration-Chemistry model is not implemented yet";
			  return 0.0;
		}
	}

private:
	int m_ns;
	double* mp_wrk1;
	double* mp_wrk2;
	double* mp_wrk3;
	double* mp_wrk4;
	double Ttr;
	double Tv;

	double const compute_source_Candler();
	double const compute_source_MMT();
};

 /**
 * Non-preferential Model according to Candler with
 *
 * with \f$ c_1 \f$ equal to 1 for non-preferential models. For more information:
 *
 * G. V. Candler, MacCormack, Computation of weakly ionized hypersonic flows in thermochemical nonequilibrium, Journal of
 * Thermophysics and Heat Transfer, 1991, 5(11):266
 *
 */

double const OmegaCV::compute_source_Candler()
{


	 // Getting Vibrational Energy
	 mixture().speciesHOverRT(NULL, NULL, NULL, mp_wrk1, NULL, NULL);

	 // Getting Production Rate
	 mixture().netProductionRates(mp_wrk2);

	 // Inner Product
	 double c1 = 1.0E0;
	 double sum = 0.E0;

	 for(int i = 0 ; i < m_ns; ++i)
		 sum += mp_wrk1[i]*mp_wrk2[i]/mixture().speciesMw(i);

	 return(c1*sum*mixture().T()*RU);
 }

 /**
 * Preferential Model according to MMT 
 */

 double const OmegaCV::compute_source_MMT()
 {
	Ttr = mixture().T(); // Translational temperature
	Tv = mixture().Tv(); // Vibrational temperature
	const int nr = mixture().nReactions();

	// Getting Production Rate
	mixture().speciesHOverRT(NULL, NULL, NULL, mp_wrk1, NULL, NULL);
	
	// Get reaction enthalpies
	std::fill(mp_wrk2, mp_wrk2+nr, 0.0);
	mixture().getReactionDelta(mp_wrk1,mp_wrk2);

	// Get molar rates of progress
	mixture().netRatesOfProgress(mp_wrk3);

	/// for r in reactions:
	// if MMT:
	/// do below
	// if Arrhenius
	//	Call netProgressRates(mp_wrk2)
	// Add reactants
	// Subtract products
	// * R * T / M
	
	//Attempt to get data for each reaction
	for(int i=0; i<nr; ++i) {
		if (dynamic_cast<const MMT*>(mixture().reactions()[1].rateLaw()) != NULL){
			const Reaction& r = mixture().reactions()[i];
			const MMT& rate = dynamic_cast<const MMT&>(*(r.rateLaw()));
			const double TD = rate.T();
			const double thetaV = rate.thetaV();
			const double aU = rate.a();
			const double Us = rate.U_s();
			const double U = 1 / (aU / Ttr + 1 / Us);
			const double TF = 1 / (1 / Tv - 1 / Ttr - 1 / U);
			const double e_vib = -1.0 * KB * (thetaV / std::exp(thetaV/TF)) - (TD / std::exp(TD/TF));
			mp_wrk4[i] = e_vib; // J units
		}
    // m_reactants.decrSpecies(mp_rop, p_wdot);
    // m_rev_prods.incrSpecies(mp_rop, p_wdot);
    // m_irr_prods.incrSpecies(mp_rop, p_wdot);		}
	};

	double Qv = 0.0;
	// const double eps = 1.0E-15;
	// const int ns = mixture().nSpecies();
	// for(int i=0; i<ns; ++i){
		// if Arrhenius: Vibrational Energy * Production Rate (Candler)
		for(int j=0; j<nr; ++j){
			if (dynamic_cast<const MMT*>(mixture().reactions()[j].rateLaw()) != NULL){
				Qv += mp_wrk1[j]; // evib_d * R_j
			}
			else { // Arrhenius reactions
				Qv += mp_wrk2[j]*mp_wrk3[j] *RU*mixture().T();
				break;
			}
		}
	// }



	/*
	Q = 0
	ind = 0
	for s in species:
		for r in reaction:
			if mixture().Y()[ind] == 0: // or < eps
				continue
			Q += e_vib(r,s) * R(s)
		ind += 1
  	*/
 
	  // Getting Vibrational Energy
	//   mixture().speciesHOverRT(NULL, NULL, NULL, mp_wrk1, NULL, NULL);

	  
 
	  // Inner Product
	//   double c1 = 1.0E0;
	//   double sum = 0.E0;
 
	//   for(int i = 0 ; i < m_ns; ++i)
	// 	  sum += mp_wrk1[i]*mp_wrk2[i]; // /mixture().speciesMw(i);
 
	  return(Qv);
  }
 

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaCV, TransferModel> omegaCV("OmegaCV");
      
    } // namespace Transfer
} // namespace Mutation
