/**
 * @file OmegaVT.cpp
 *
 * @brief Implementation of OmegaVT.
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

#include "MillikanWhite.h"
#include "HarmonicOscillator.h"
#include "Mixture.h"
#include "TransferModel.h"
#include <cmath>

#include <Eigen/Dense>

using namespace Mutation;
using namespace Mutation::Kinetics;

namespace Mutation {
    namespace Transfer {

/**
 * Combines vibrator energy model with relaxation time model to form a 
 * vibration-translation energy relaxation model for a single vibrator.
 */
class Vibrator
{
public:
    Vibrator(const Thermodynamics::HarmonicOscillator&, const MillikanWhiteModel&);

    /// Index of species in mixture.
    size_t speciesIndex() const { 
        return m_relaxation_model.speciesIndex(); 
    }

    /// Molecular weight of the vibrator in kg/mol
    double molecularWeight() const { 
        return m_relaxation_model.molecularWeight(); 
    }
    
    /// Energy in K of the vibrator at the given temperature.
    double energy(double T) const { 
        return m_energy_model.energy(T); 
    }

    /// Relaxation time of the vibrator.
    /// @todo Find better solution to limit dependency to just thermodynamic
    /// state information, not whole mixture.
    double relaxationTime(const class Thermodynamics::Thermodynamics& thermo) const { 
        return m_relaxation_model.relaxationTime(thermo);
    }

private:
    Thermodynamics::HarmonicOscillator m_energy_model;
    MillikanWhiteModel m_relaxation_model;
};


Vibrator::Vibrator(
    const Thermodynamics::HarmonicOscillator& ho, const MillikanWhiteModel& mw
) : m_energy_model(ho), m_relaxation_model(mw) { }


/**
 * @brief Vibration-translation energy relaxation model.
 * 
 * Implements the Millikan-White model with the Park correction based on
 * Eqs. (55-58) in NASA-TP-2867, Gnoffo et al. (1989).
 */
class OmegaVT : public TransferModel
{
public:
    OmegaVT(Mixture&);
    // OmegaVT(Mixture&) : TransferModel(mix)
    // {
    //     m_ns = mixture().nSpecies();
    //     hv_T = new double [m_ns];
    //     hv_Tv = new double [m_ns];
    //     hel_T = new double [m_ns];
    //     hel_Tv = new double [m_ns];
    // };
    // virtual ~OmegaVT() = default;
    ~OmegaVT()
	{
		delete [] hv_T;
		delete [] hv_Tv;
		delete [] hel_T;
		delete [] hel_Tv;
	};
    // virtual double source() override;
    double source()
	{
		static int i_transfer_model = 0;
		const int nr = mixture().nReactions();
		for(int i=0; i<nr; ++i) {
			i_transfer_model = (dynamic_cast<const MMT*>(mixture().reactions()[i].rateLaw()) != NULL);
			if (i_transfer_model){
				break;
			}
		}
		// static int i_transfer_model = (dynamic_cast<const MMT*>(mixture().reactions()[1].rateLaw()) != NULL);
		switch (i_transfer_model){
		   case 0:
			  return compute_source_non_preferential();
		   case 1:
              return compute_source_preferential();
		  break;
		   default:
			  std::cerr << "The selected vibration-translation energy relaxation model is not implemented yet";
			  return 0.0;
		}
	}
private:
    std::vector<Vibrator> m_vibrators;
    int m_ns;
	double* hv_T;
	double* hv_Tv;
	double* hel_T;
	double* hel_Tv;

    double compute_source_non_preferential();
	double compute_source_preferential();
};


// Register the transfer model
Utilities::Config::ObjectProvider<OmegaVT, TransferModel> omegaVT("OmegaVT");


OmegaVT::OmegaVT(Mixture& mix) : 
    TransferModel(mix)
{
    

    Thermodynamics::HarmonicOscillatorDB hodb;
    MillikanWhiteModelDB mwdb(mix);

    for (const auto& species: mix.species())
    {
        if (species.type() == Thermodynamics::MOLECULE)
        {
            // std::cout << "Species: " << species.name() << std::endl;
            auto ho = hodb.create(species.name());
            // std::cout << "Temperatures: " << ho.characteristicTemperatures()[0] << std::endl;
            auto mw = mwdb.create(species.name(), ho.characteristicTemperatures()[0]);
            m_vibrators.emplace_back(ho, mw);
        }   
        else
        {
            // std::cout << "Species: " << species.name() << std::endl;
            auto ho = hodb.create(species.name());
            // std::cout << "Temperatures: " << 0 << std::endl;
            auto mw = mwdb.create(species.name(), 0);
            m_vibrators.emplace_back(ho, mw);
        }   
    }
    
}


// double OmegaVT::source()
// {
//     auto Y = mixture().Y();
//     auto T = mixture().T();
//     auto Tv = mixture().Tv();
//     m_ns = mixture().nSpecies();
//     hv_T = new double[m_ns](); 
//     hv_Tv = new double[m_ns](); 
//     hel_T = new double[m_ns](); 
//     hel_Tv = new double[m_ns](); 

//     // double src = 0.0;

//     // for (auto& vibrator: m_vibrators)
//     // {
//     //     const auto iv = vibrator.speciesIndex();
//     //     const auto tau = vibrator.relaxationTime(mixture());
//     //     const auto mw = vibrator.molecularWeight();
//     //     src += Y[iv]*(vibrator.energy(T) - vibrator.energy(Tv))/(mw*tau);
//     // }
//     // std::cout << src << std::endl;
//     // return src * mixture().density() * RU;


//     auto X = mixture().X();
//     const int ns = mixture().nSpecies();
//     mixture().speciesHOverRT(T, T, T, T, T, NULL, NULL, NULL, hv_T, hel_T, NULL);
//     mixture().speciesHOverRT(T, Tv, T, Tv, Tv, NULL, NULL, NULL, hv_Tv, hel_Tv, NULL);

//     double src = 0.0;
//     double tau_sum = 0.0;
//     double X_sum = 0.0;
//     double ET, ETv;

//     for (auto& vibrator: m_vibrators)
//     {
//         const auto iv = vibrator.speciesIndex();
//         const auto mw = vibrator.molecularWeight();

//         // ET = T * RU * (hv_T[iv] + hel_T[iv]); // J/mol
//         // ETv = T * RU * (hv_Tv[iv] + hel_Tv[iv]); // J/mol
        
//         // std::cout << "T: " << ET << " Tv: " << ETv << std::endl;
//         // src += Y[iv] * (ET - ETv) / (mw);
//         // if (vibrator.energy(Tv) > 1) // molecule
//         if (mixture().species()[iv].type() == Thermodynamics::MOLECULE)
//         {
//             const auto tau = vibrator.relaxationTime(mixture());
//             src += Y[iv]*(vibrator.energy(T) - vibrator.energy(Tv))/(mw) * RU;
//             // std::cout << "tau= " << tau << " X= " << X[iv] << std::endl;
//             tau_sum += X[iv] / tau;
//             X_sum += X[iv];
//             // std::cout << "tau_sum: " << tau_sum << " X_sum: " << X_sum << std::endl;
//         }
//         else {
//             // std::cout << "T: " << hv_T[iv] << " Tv: " << hv_Tv[iv] << std::endl;
//             ET = T * RU * (hv_T[iv] + hel_T[iv]); // J/mol
//             ETv = T * RU * (hv_Tv[iv] + hel_Tv[iv]); // J/mol
//             src += Y[iv] * (ET - ETv) / (mw);
//         }
//     }
//     double tau_denom = X_sum / tau_sum;
//     src = src / tau_denom;
//     // std::cout << src << std::endl;

//     return src * mixture().density();// * RU;
// }

double OmegaVT::compute_source_non_preferential()
{
    auto Y = mixture().Y();
    auto T = mixture().T();
    auto Tv = mixture().Tv();

    double src = 0.0;

    for (auto& vibrator: m_vibrators)
    {
        const auto iv = vibrator.speciesIndex();
        const auto tau = vibrator.relaxationTime(mixture());
        const auto mw = vibrator.molecularWeight();
        src += Y[iv]*(vibrator.energy(T) - vibrator.energy(Tv))/(mw*tau);
    }

    return src * mixture().density() * RU;
}




double OmegaVT::compute_source_preferential()
{
    auto Y = mixture().Y();
    auto T = mixture().T();
    auto Tv = mixture().Tv();
    m_ns = mixture().nSpecies();
    hv_T = new double[m_ns](); 
    hv_Tv = new double[m_ns](); 
    hel_T = new double[m_ns](); 
    hel_Tv = new double[m_ns](); 

    // double src = 0.0;

    // for (auto& vibrator: m_vibrators)
    // {
    //     const auto iv = vibrator.speciesIndex();
    //     const auto tau = vibrator.relaxationTime(mixture());
    //     const auto mw = vibrator.molecularWeight();
    //     src += Y[iv]*(vibrator.energy(T) - vibrator.energy(Tv))/(mw*tau);
    // }
    // std::cout << src << std::endl;
    // return src * mixture().density() * RU;


    auto X = mixture().X();
    const int ns = mixture().nSpecies();
    mixture().speciesHOverRT(T, T, T, T, T, NULL, NULL, NULL, hv_T, hel_T, NULL);
    mixture().speciesHOverRT(T, Tv, T, Tv, Tv, NULL, NULL, NULL, hv_Tv, hel_Tv, NULL);

    double src = 0.0;
    double tau_sum = 0.0;
    double X_sum = 0.0;
    double ET, ETv;

    for (auto& vibrator: m_vibrators)
    {
        const auto iv = vibrator.speciesIndex();
        const auto mw = vibrator.molecularWeight();

        // ET = T * RU * (hv_T[iv] + hel_T[iv]); // J/mol
        // ETv = T * RU * (hv_Tv[iv] + hel_Tv[iv]); // J/mol
        
        // std::cout << "T: " << ET << " Tv: " << ETv << std::endl;
        // src += Y[iv] * (ET - ETv) / (mw);
        // if (vibrator.energy(Tv) > 1) // molecule
        if (mixture().species()[iv].type() == Thermodynamics::MOLECULE)
        {
            const auto tau = vibrator.relaxationTime(mixture());
            src += Y[iv]*(vibrator.energy(T) - vibrator.energy(Tv))/(mw) * RU;
            // std::cout << "tau= " << tau << " X= " << X[iv] << std::endl;
            tau_sum += X[iv] / tau;
            X_sum += X[iv];
            // std::cout << "tau_sum: " << tau_sum << " X_sum: " << X_sum << std::endl;
        }
        else {
            // std::cout << "T: " << hv_T[iv] << " Tv: " << hv_Tv[iv] << std::endl;
            ET = T * RU * (hv_T[iv] + hel_T[iv]); // J/mol
            ETv = T * RU * (hv_Tv[iv] + hel_Tv[iv]); // J/mol
            src += Y[iv] * (ET - ETv) / (mw);
        }
    }
    double tau_denom = X_sum / tau_sum;
    src = src / tau_denom;
    // std::cout << src << std::endl;

    return src * mixture().density();// * RU;
}





/**
 * Represents a coupling between vibrational and translational energy modes.
 */
// class OldOmegaVT : public TransferModel
// {
// public:

//     OldOmegaVT(Mixture& mix)
//         : TransferModel(mix), m_mw(mix)
//     {
//         m_const_Park_correction = std::sqrt(PI*KB/(8.E0*NA));
//         m_ns              = m_mixture.nSpecies();
//         m_transfer_offset = m_mixture.hasElectrons() ? 1 : 0;

//         mp_Mw = new double [m_ns];
//         for(int i = 0; i < m_ns; ++i)
//             mp_Mw[i] = m_mixture.speciesMw(i);
//         mp_hv = new double [m_ns];
//         mp_hveq = new double [m_ns];
//     }

//     virtual ~OldOmegaVT()
//     {
//         delete [] mp_Mw;
//         delete [] mp_hv;
//         delete [] mp_hveq;
//     }

//     /**
//      * Computes the source terms of the Vibration-Translational energy transfer in \f$ [J/(m^3\cdot s)] \f$
//      * using a Landau-Teller formula taking into account Park's correction (default; can be disabled by making zero the appropriate flag, see below):
//      * \f[ \Omega^{VT}_m = \rho_m \frac{e^V_m\left(T\right)-e^V_m\left(T_{vm}\right)} {\tau^{VT}_m} \f]
//      *
//      * The average relaxation time \f$ \tau^{VT}_m \f$ is given by the expression:
//      *
//      * \f[ \tau^{VT}_m = \frac{ \sum_{j \in \mathcal{H}} \rho_j / M_j}{ \sum_{j \in \mathcal{H}} \rho_j / (M_j \tau^{VT}_{mj}) } \f]
//      *
//      * More information about the above model can be found in @cite Panesi 2008
//      *
//      * @todo compare the average relaxation time with the one suggested by Park93
//      *
//      */

//     double source()
//     {
//         const double * p_Y = m_mixture.Y();
//         double rho = m_mixture.density();
//         double T = m_mixture.T();
//         double Tv = m_mixture.Tv();

//         m_mixture.speciesHOverRT(T, T, T, T, T, NULL, NULL, NULL, mp_hveq, NULL, NULL);
//         m_mixture.speciesHOverRT(T, Tv, T, Tv, Tv, NULL, NULL, NULL, mp_hv, NULL, NULL);

//         int inv = 0;
//         double src = 0.0;
//         for (int iv = 0; iv-inv < m_mw.nVibrators(); ++iv){
//             if(m_mixture.species(iv).type() != Mutation::Thermodynamics::MOLECULE){
//                 inv++;
//             } else {
//                 src += p_Y[iv]*rho*RU*T/mp_Mw[iv]*(mp_hveq[iv] - mp_hv[iv])/compute_tau_VT_m(iv-inv);
//             }
//         }
//         return src;
//     }

// private:

//     MillikanWhite m_mw;

//    /**
//     * @brief This function computes the Millikan and White relaxation time
//     * for each diatomic collision
//     *
//     * @param vibrator index
//     * @param partner index
//     *
//     * @return Millikan and White relaxation time \tau_{m,j}^{MW}
//     */

//    inline double const compute_tau_VT_mj(int const, int const);

//     /**
//      * @brief Computes the frequency average over heavy particles.
//      * @param vibrator index
//      *
//      * @return
//      */
//     double compute_tau_VT_m(int const);

//     /**
//      * @brief This function computes the Park correction
//      * for vibrational-translational energy transfer
//      * for each collision pair.
//      *
//      * @param vibrator index
//      * @param partner index
//      *
//      * @return Park correction \tau_{m,j}^P
//      */

//     inline double const compute_Park_correction_VT(int const, int const);

//     /**
//      * Necessary variables
//      */
//     int m_ns;
//     int m_transfer_offset;

//     double* mp_Mw;
//     double* mp_hv;
//     double* mp_hveq;

//     double m_const_Park_correction;
// };
      
// // Implementation of the Vibrational-Translational Energy Transfer.
      
// inline double const OmegaVT::compute_Park_correction_VT(int const i_vibrator, int const i_partner)
// {
//     // Limiting cross section for Park's Correction
//     double P = m_mixture.P();
//     double T = m_mixture.T();

//     double sigma;
//     if (T > 20000.0) {
//       sigma = m_mw[i_vibrator].omega() * 6.25 ; // 6.25 = (50000/20000)^2
//     } else {
//       sigma = m_mw[i_vibrator].omega() *(2.5E9/(T*T));
//     }
    
//     return(m_const_Park_correction * sqrt(m_mw[i_vibrator][i_partner].mu()*T)/(sigma*P));
// }

 
// inline double const OmegaVT::compute_tau_VT_mj(int const i_vibrator, int const i_partner)
// {
// //    Enable in the future for multiple vibrational temperatures
//       double P = m_mixture.P();
//       double T = m_mixture.T();

//       return( exp( m_mw[i_vibrator][i_partner].a() * (pow(T,-1.0/3.0) - m_mw[i_vibrator][i_partner].b()) -18.421) * ONEATM / P );
// }
      
// double OmegaVT::compute_tau_VT_m(int const i_vibrator)
// {
//     const double * p_Y = m_mixture.Y();
    
//     double sum1 = 0.0;
//     double sum2 = 0.0;
    
//     // Partner offset
//     for (int i_partner = m_transfer_offset; i_partner < m_ns; ++i_partner){
//         double tau_j = compute_tau_VT_mj(i_vibrator, i_partner-m_transfer_offset) + compute_Park_correction_VT(i_vibrator, i_partner-m_transfer_offset);
//         sum1 += p_Y[i_partner]/(mp_Mw[i_partner]);
//         sum2 += p_Y[i_partner]/(mp_Mw[i_partner]*tau_j);
//     }
  
//     return(sum1/sum2);
// }


// // Register the transfer model
// Utilities::Config::ObjectProvider<
//     OmegaVT, TransferModel> omegaVT("OmegaVT");

    } // namespace Transfer
} // namespace Mutation 
