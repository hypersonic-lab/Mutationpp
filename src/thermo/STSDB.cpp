/**
 * @file STSDB.cpp
 *
 * @brief Provides a N3 bins thermodynamic database.
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

#include "Constants.h"
#include "ThermoDB.h"
#include "Species.h"
#include "ParticleRRHO.h"
#include "AutoRegistration.h"
#include "Functors.h"
#include "LookupTable.h"
#include "Utilities.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>

using namespace std;
using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

// Simple for loop over the species just to make life a little easier
#define LOOP(__op__)\
for (int i = 0; i < m_ns; ++i) {\
    __op__ ;\
}

#define LOOP_ATOMS(__op__)\
for (int i = 0, j = 0; i < m_na; ++i) {\
    j = mp_indices[i];\
    __op__ ;\
} // j = mp_indices[i]

// Loops over molecules. Inside the loop, the index i is zero based and indexes
// internal molecular data.  The index j is the index corresponding to the 
// original species data.
#define LOOP_MOLECULES(__op__)\
for (int i = 0, j = 0; i < m_nm; ++i) {\
    j = mp_indices[m_na+i];\
    __op__ ;\
}
    // j = mp_indices[m_na+i];

// Loops over heavy particles (non electron species).  Inside the loop, index i
// is zero based and indexes internal molecular data. The index j is the index
// corresponding to the original species data.
#define LOOP_HEAVY(__op__)\
for (int i = 0, j = 0; i < m_na + m_nm; ++i) {\
    j = (m_has_electron ? i+1 : i);\
    __op__ ;\
}

typedef struct {
    double ln_omega_t;  // ln(omega^(2 / L) * theta)
    double linearity;   // L / 2
} RotData;

// typedef struct {
//     double g;           // degeneracy
//     double theta;       // characteristic temperature
// } ElecLevel;

// typedef struct {
//     unsigned int offset;
//     unsigned int nheavy;
//     unsigned int nlevels;
//     int* p_nelec;
//     ElecLevel* p_levels;
// } ElectronicData;


/**
 * A thermodynamic database that uses the Rigid-Rotator Harmonic-Oscillator
 * model for computing species thermodynamic properties.  See the individual
 * thermodynamic functions for specific descriptions of the model.
 */
    
/**** TODO: Figure out what states exist
         Establish array in that form
         Read in energy constants from file
*/
class STSDB : public ThermoDB
{
public:

    STSDB(int arg) : ThermoDB(298.15, 101325.0)
//          ,m_has_electron(false),
//           m_use_tables(false),
//           m_last_bfacs_T(0.0)
    { }
    /**
     * Destructor.
     */
    ~STSDB()
    {
//         delete [] mp_lnqtmw;
        delete [] mp_hform;
        delete [] mp_indices;
        delete [] mp_rot_data;
        delete [] mp_nvib;
        delete [] mp_vib_temps;

//         delete [] m_elec_data.p_nelec;
//         delete [] m_elec_data.p_levels;
        delete [] mp_part_sst;
//         delete [] mp_el_bfacs;

//         if (m_use_tables) {
//             delete mp_el_bfac_table;
//         }
    }

    /**
     * Computes the unitless species specific heat at constant pressure
     * \f$ C_{P,i} / R_U\f$ in thermal nonequilibrium.
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param cp   - on return, the array of species non-dimensional \f$C_P\f$'s
     * @param cpt  - if not NULL, the array of species translational \f$C_P\f$'s
     * @param cpr  - if not NULL, the array of species rotational \f$C_P\f$'s
     * @param cpv  - if not NULL, the array of species vibrational \f$C_P\f$'s
     * @param cpel - if not NULL, the array of species electronic \f$C_P\f$'s
     *
     * @todo Compute \f$C_P\f$ directly instead of using finite-differencs.
     */
    void cp(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const cp, double* const cpt, double* const cpr, 
        double* const cpv, double* const cpel)
    {
        for (int i = 0; i < m_ns; i++)
        {
            cp[i] = 0.0;
        }

        // Special case if we only want total Cp
        if (cp != NULL && cpt == NULL && cpr == NULL && cpv == NULL && 
            cpel == NULL)
        {
            cpT(cp, Eq());
            cpR(cp, PlusEq());
            cpV(Tv, cp, PlusEq());
            cpE(Tel, cp, PlusEq());
            return;
        }
        
        // Otherwise we have to compute each component directly.
        // Translation
        if (cpt == NULL) {
            if (cp != NULL)
                //cpT(cp, PlusEq());
                cpT(cp, Eq());
        } else {
            cpT(cpt, Eq());
            if (cp != NULL) 
                LOOP(cp[i] = cpt[i]);
        }
        
        // Rotation
        if (cpr == NULL) {
            if (cp != NULL)
                cpR(cp, PlusEq());
        } else {
            cpR(cpr, Eq());
            if (cp != NULL) 
                LOOP_MOLECULES(cp[j] += cpr[j]);
        }
        
        // Vibration
        if (cpv == NULL) {
            if (cp != NULL)
                cpV(Tv, cp, PlusEq());
        } else {
            cpV(Tv, cpv, Eq());
            if (cp != NULL) 
                LOOP_MOLECULES(cp[j] += cpv[j]);
        }
        
        // Electronic
        if (cpel == NULL) {
            if (cp != NULL)
                cpE(Tel, cp, PlusEq());
        } else {
            cpE(Tel, cpel, Eq());
            if (cp != NULL)
                LOOP_HEAVY(cp[j] += cpel[j]);
        }
    }

    /**
     * Computes the unitless species enthalpy \f$ h_i/R_U T_h\f$ of each
     * species in thermal nonequilibrium, which is non-dimensionalized by the 
     * heavy particle translational temperature. 
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param h   - on return, the array of species non-dimensional enthalpies
     * @param ht  - if not NULL, the array of species translational enthalpies
     * @param hr  - if not NULL, the array of species rotational enthalpies
     * @param hv  - if not NULL, the array of species vibrational enthalpies
     * @param hel - if not NULL, the array of species electronic enthalpies
     * @param hf  - if not NULL, the array of the species formation enthalpies
     */
    void enthalpy(
        double Th, double Te, double Tr, double Tv, double Tel, double* const h,
        double* const ht, double* const hr, double* const hv, double* const hel,
        double* const hf)
    {
        // Setting to zero
        for (int i = 0; i < m_ns; i++)
        {
            h[i] = 0.0;
        }

        // Special case where we only want the total enthalpy
        if (ht == NULL && hr == NULL && hv == NULL && hel == NULL && 
            hf == NULL && h != NULL) 
        {
            hT(Th, Te, h, Eq());
            hR(Tr, h, PlusEq());
            hV(Tv, h, PlusEq());
            hE(Tel, h, PlusEq());
            hF(h, PlusEq());
            LOOP(h[i] /= Th);
            return;
        }
        
        // Otherwise selectively choose what we want
        // Translational enthalpy
        if (ht == NULL) {
            //hT(Th, Te, h, Eq());
            if (h != NULL)
                hT(Th, Te, h, EqDiv(Th));
        } else {
            hT(Th, Te, ht, EqDiv(Th));
            if (h != NULL)
                LOOP(h[i] = ht[i]);
        }
        
        // Rotatonal enthalpy
        if (hr == NULL) {
            if (h != NULL)
                hR(Tr, h, PlusEqDiv(Th));
        } else {
            LOOP(hr[i] = 0.0);
            hR(Tr, hr, EqDiv(Th));
            if (h != NULL)
                LOOP_MOLECULES(h[j] += hr[j]);
        }
        
        // Vibrational enthalpy
        if (hv == NULL) {
            if (h != NULL)
                hV(Tv, h, PlusEqDiv(Th));
        } else {
            LOOP(hv[i] = 0.0);
            hV(Tv, hv, EqDiv(Th));
            if (h != NULL)
                LOOP_MOLECULES(h[j] += hv[j]);
        }
            
        // Electronic enthalpy
        if (hel == NULL) {
            if (h != NULL)
                hE(Tel, h, PlusEqDiv(Th));
        } else {
            LOOP(hel[i] = 0.0);
            hE(Tel, hel, EqDiv(Th));
            if (h != NULL)
                LOOP(h[i] += hel[i]);
        }
        
        // Formation enthalpy
        if (hf == NULL) {
            if (h != NULL)
                hF(h, PlusEqDiv(Th));
        } else {
            hF(hf, EqDiv(Th));
            if (h != NULL)
                LOOP(h[i] += hf[i]);
        }
    }

    /**
     * Computes the unitless species entropy \f$s_i/R_u\f$ allowing for thermal 
     * nonequilibrium.
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param s   - on return, the array of species non-dimensional entropies
     * @param st  - if not NULL, the array of species translational entropies
     * @param sr  - if not NULL, the array of species rotational entropies
     * @param sv  - if not NULL, the array of species vibrational entropies
     * @param sel - if not NULL, the array of species electronic entropies
     */
    void entropy(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const s, double* const st, double* const sr, double* const sv,
        double* const sel)
    {
        // Special case where we only want the total entropy
        if (st == NULL && sr == NULL && sv == NULL && sel == NULL) {
            sT(Th, Te, P, s, Eq());
            sR(Tr, s, PlusEq());
            sV(Tv, s, PlusEq());
            sE(Tel, s, PlusEq());
            
            // Include spin contribution for free electron entropy
            if (m_has_electron)
                s[0] += std::log(2.0);
            
            return;
        }
        
        // Otherwise collect individual components
        // Translational entropy
        if (st == NULL)
            sT(Th, Te, P, s, Eq());
        else {
            sT(Th, Te, P, st, Eq());
            LOOP(s[i] = st[i]);
        }
        
        // Rotational entropy
        if (sr == NULL)
            sR(Tr, s, PlusEq());
        else {
            LOOP(sr[i] = 0.0);
            sR(Tr, sr, Eq());
            LOOP_MOLECULES(s[j] += sr[j]);
        }
        
        // Vibrational entropy
        if (sv == NULL)
            sV(Tv, s, PlusEq());
        else {
            LOOP(sv[i] = 0.0);
            sV(Tv, sv, Eq());
            LOOP_MOLECULES(s[j] += sv[j]);
        }
        
        // Electronic entropy
        if (sel == NULL)
            sE(Tel, s, PlusEq());
        else {
            LOOP(sel[i] = 0.0);
            sE(Tel, sel, Eq());
            LOOP(s[i] += sel[i]);
        }
        
        // Include spin contribution for free electron entropy
        if (m_has_electron)
            s[0] += std::log(2.0);
    }
    
    /**
     * Computes the unitless Gibbs free energy of each species i, 
     * \f$G_i / R_u T_h\f$ where \f$G_i\f$ is non-dimensionalized by the heavy
     * particle translational temperature.
     *
     * @todo Compute the individual components of the Gibbs function directly
     * instead of H - TS.
     */
    void gibbs(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const g, double* const gt, double* const gr, double* const gv,
        double* const gel)
    {        
        // First compute the non-dimensional enthalpy
        enthalpy(Th, Te, Tr, Tv, Tel, g, NULL, NULL, NULL, NULL, NULL);

        // Subtract the entropies
        sT(Th, Te, P, g, MinusEq());
        sR(Tr, g, MinusEq());
        sV(Tv, g, MinusEq());
        sE(Tel, g, MinusEq());

        // Account for spin of free electrons
        if (m_has_electron)
            g[0] -= std::log(2.0);
    }
    
protected:
    /**
     * Loads all of the species from the RRHO database.
     */
    virtual void loadAvailableSpecies(std::list<Species>& species)
    {
        IO::XmlDocument species_doc(databaseFileName("species.xml", "thermo"));
        IO::XmlElement::const_iterator species_iter = species_doc.root().begin();

        for ( ; species_iter != species_doc.root().end(); ++species_iter) {
            // Add the species to the list
            species.push_back(*species_iter);


            // We can also add all of the excited states as implicitly defined
            // species
            IO::XmlElement::const_iterator rrho_iter =
                species_iter->findTagWithAttribute(
                    "thermodynamics", "type", "STS");

            if (rrho_iter == species_iter->end()){
                continue;}

            Species& ground_state = species.back();
            ParticleRRHO rrho(*rrho_iter);
            // for (size_t i = 0; i < rrho.nElectronicLevels(); ++i)
            // {
            //     species.push_back(Species(ground_state, i));
            // }
            for (size_t i = 0; i < rrho.nVibrationalLevels(); ++i)
            {
                species.push_back(Species(ground_state, i));
            }
        }

        // // @todo: 1/18/2023
        // // Find the species in the species.xml database. Not needed.

        // // @debug
        // for (const auto& sp : species) {
        //     std::cout << "Species name = " << sp.name() << std::endl;
        // }
        // double in; std::cin >> in;
    }

    /**
     * Load thermodynamic data from the species list.
     */
    virtual void loadThermodynamicData()
    {
        m_ns = species().size();
		m_has_electron = (species()[0].type() == ELECTRON);

		// Load the RRHO models for each of the needed species
        IO::XmlDocument species_doc(databaseFileName("species.xml", "thermo"));
        
        vector<ParticleRRHO> rrhos;
        map<std::string, const ParticleRRHO*> to_expand;
        
        for (int i = 0; i < m_ns; ++i) {
            if (species()[i].name() == species()[i].groundStateName()) {
                rrhos.push_back(*(species_doc.root().findTagWithAttribute(
                    "species", "name", species()[i].groundStateName())->
                        findTagWithAttribute("thermodynamics", "type", "STS")));
            }
            else {
                const ParticleRRHO* p_rrho = to_expand[species()[i].groundStateName()];
                if (p_rrho == NULL) {
                    p_rrho = new ParticleRRHO(
                        *(species_doc.root().findTagWithAttribute(
                            "species", "name", species()[i].groundStateName())->
                                findTagWithAttribute("thermodynamics", "type", "STS")));
                    to_expand[species()[i].groundStateName()] = p_rrho;
                }
                rrhos.push_back(ParticleRRHO(*p_rrho, species()[i].level()));
            }
        }

		map<std::string, const ParticleRRHO*>::iterator iter =
		to_expand.begin();
        while (iter != to_expand.end()) {
            delete iter->second;
            iter++;
        }

		vector<int> atom_indices;
        vector<int> molecule_indices;

        LOOP(
            switch(species()[i].type()) {
                case ATOM:
                    atom_indices.push_back(i);
                    break;
                case MOLECULE:
                    molecule_indices.push_back(i);
                    break;
                default:
                    break;
            }
        )

        m_na = atom_indices.size();
        m_nm = molecule_indices.size();

		// Order the atoms first followed by the molecules
        mp_indices = new int [m_na + m_nm];
        copy(atom_indices.begin(), atom_indices.end(), mp_indices);
        copy(molecule_indices.begin(), molecule_indices.end(), mp_indices+m_na);

        mp_hform = new double [m_ns];
        LOOP(mp_hform[i] = rrhos[i].formationEnthalpy() / RU)

		// Store the molecule's rotational energy parameters
        mp_rot_data = new RotData [m_nm];
        LOOP_MOLECULES(
            const ParticleRRHO& rrho = rrhos[j];
            int linear = rrho.linearity();
            mp_rot_data[i].linearity  = linear / 2.0;
            mp_rot_data[i].ln_omega_t = 
                std::log(rrho.rotationalTemperature()) + 2.0 / linear *
                std::log(rrho.stericFactor());
        )

        // Store the vibrational temperatures of all the molecules in a compact
        // form
        mp_nvib = new int [m_nm];
        int nvib = 0;
        LOOP_MOLECULES(
            mp_nvib[i] = rrhos[j].nVibrationalLevels();
            nvib += mp_nvib[i];
        )
        
        mp_vib_temps = new double [nvib];
        int ilevel = 0;
        LOOP_MOLECULES(
            const ParticleRRHO& rrho = rrhos[j];
            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++)
                mp_vib_temps[ilevel] = rrho.vibrationalEnergy(0);
        )

        mp_part_sst = new double [m_ns];
        hT(Tss, Tss, mp_part_sst, Eq());
        hR(Tss, mp_part_sst, PlusEq());
        hV_f(Tss, mp_part_sst, PlusEq());
        hE(Tss, mp_part_sst, PlusEq());
    }

private:
    typedef Equals<double> Eq;
    typedef EqualsYDivAlpha<double> EqDiv;
    typedef PlusEqualsYDivAlpha<double> PlusEqDiv;
    typedef PlusEquals<double> PlusEq;
    typedef MinusEquals<double> MinusEq;

    int m_ns;
    int m_na;
    int m_nm;

	bool m_has_electron;

	double* mp_part_sst;
	double* mp_hform;

	int*       mp_indices;
	RotData*   mp_rot_data;

    int*       mp_nvib;
    double*    mp_vib_temps;

	double g0_O2 = 3.0;
	double g1_O2 = 2.0;
	double theta_1_O2 = 11900;
	double g0_O = 5.0;
	double g1_O = 4.0;
	double theta_1_O = 270;
	double cm2J = 1.98630e-23;
    double Tss = standardTemperature();
    double ThetaR = 2.08; //char temp rot O2

	std::array<double, 47> m_energy = { // Varandas energy ladder cm-1
    786.0234, 2343.573, 3881.3038,
	5398.5964, 6894.8782, 8369.5118,
	9822.0457, 11251.754, 12658.072,
	14040.4352, 15398.3596, 16731.2,
	18038.3111, 19319.2091, 20573.1679,
	21799.7037, 22998.2518, 24168.0058,
	25308.4816, 26419.0341, 27498.9373,
	28547.546, 29564.2148, 30548.1374,
	31498.6683, 32414.9205, 33296.168,
	34141.6041, 34950.3419, 35721.3326,
    36453.8504, 37146.6854, 37798.8698,
	38409.2744, 38976.5281, 39499.5822,
	39976.9044, 40407.0429, 40788.6264,
	41119.9613, 41399.6763, 41626.4003,
	41799.6494, 41920.3914, 41993.5464,
	42029.6803, 42042.9885
	};

    /**
     * Computes the translational Cp/Ru for each species.
     */
    template <typename OP>
    void cpT(double* const cp, const OP& op) {
        LOOP(op(cp[i], 2.5));
    }

    /**
     * Computes the rotational Cp/Ru for each species.
     */
    template <typename OP>
    void cpR(double* const cp, const OP& op) {
        op(cp[0], 0.0);
        LOOP_ATOMS(op(cp[j], 0.0));
        LOOP_MOLECULES(op(cp[j], 1.0));
    }

    /**
     * Computes the vibratinoal Cp/Ru for each species.
     */
    template <typename OP>
    void cpV(double Tv, double* const cp, const OP& op) {
        op(cp[0], 0.0);
        LOOP_ATOMS(op(cp[j], 0.0));
        LOOP_MOLECULES(op(cp[j], 0.0));
    }

	    template <typename OP>
    void cpE(double T, double* const p_cp, const OP& op)
    {
		// double g0_O2 = 3.0;
		// double g1_O2 = 2.0;
		// double theta_1_O2 = 11900;
		// double g0_O = 5.0;
		// double g1_O = 4.0;
		// double theta_1_O = 270;

		op(p_cp[0], 0.0);

		for (int i = 0; i < m_ns; i++){
			if (i == 0){
				op(p_cp[i],pow((theta_1_O/T),2.0) * (g1_O/g0_O * exp(-theta_1_O / T)) / (pow(1.0+g1_O/g0_O * exp(-theta_1_O / T),2.0))); // Ground state
			}
			else {
				op(p_cp[i],pow((theta_1_O2/T),2.0) * (g1_O2/g0_O2 * exp(-theta_1_O2 / T)) / (pow(1.0+g1_O2/g0_O2 * exp(-theta_1_O2 / T),2.0)));
			}
		}
    }

    /**
     * Computes the translational enthalpy of each species in K.
     */
    template <typename OP>
    void hT(double T, double Te, double* const h, const OP& op) {
        if (m_has_electron)
            op(h[0], 2.5 * Te);
        LOOP_HEAVY(op(h[j], 2.5 * T))
    }

    /**
     * Computes the rotational enthalpy of each species in K.
     */
    template <typename OP>
    void hR(double T, double* const h, const OP& op) {
        LOOP_MOLECULES(op(h[j], mp_rot_data[i].linearity * T))
    }

    /**
     * Computes the vibrational enthalpy of each species in K.
     */
    template <typename OP>
    void hV(double T, double* const h, const OP& op) {
        if (T < 10.0) {
            LOOP_MOLECULES(op(h[j], 0.0));
        } else {
            LOOP_MOLECULES(op(h[j], mp_vib_temps[i] * cm2J / (KB)))
            // LOOP_MOLECULES(op(h[j], m_energy[i] * cm2J / (KB)))
        }
    }

    template <typename OP>
    void hV_f(double T, double* const h, const OP& op) {
        if (T < 10.0) {
            LOOP_MOLECULES(op(h[j], 0.0));
        } else {
            LOOP_MOLECULES(op(h[j], mp_vib_temps[0] * cm2J / (KB)))
            // LOOP_MOLECULES(op(h[j], m_energy[0] * cm2J / (KB)))
        }
    }

	/**
     * Computes the electronic enthalpy of each species in K and applies the
     * value to the enthalpy array using the given operation.
     */
    template <typename OP>
    void hE(double T, double* const p_h, const OP& op)
    {

		op(p_h[0], 0.0);

		for (int i = 0; i < m_ns; i++){
			if (i == 0){
				op(p_h[i],((theta_1_O) * g1_O/g0_O * exp(-theta_1_O / T)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / T))); // Ground state
			}
			else {
				op(p_h[i],((theta_1_O2) * g1_O2/g0_O2 * exp(-theta_1_O2 / T)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / T)));
			}
		}
    }

	/**
     * Computes the formation enthalpy of each species in K.
     */
    template <typename OP>
    void hF(double* const h, const OP& op) {
        LOOP(op(h[i], mp_hform[i] - mp_part_sst[i]))
    }


    /**
     * Computes the unitless translational entropy of each species.
     */
    template <typename OP>
    void sT(double Th, double Te, double P, double* const s, const OP& op) {
        // double fac = 2.5 * (1.0 + std::log(Th)) - std::log(P);
        // if (m_has_electron)
        //     op(s[0], 2.5 * std::log(Te / Th) + fac + mp_lnqtmw[0]);
        for (int i = (m_has_electron ? 1 : 0); i < m_ns; ++i)
            op(s[i], 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5);
    }

    /**
     * Computes the unitless rotational entropy of each species.
     */
    template <typename OP>
    void sR(double T, double* const s, const OP& op) {
        LOOP_MOLECULES(
            op(s[j], (log(T / (2 * ThetaR)) + 1.0));
        )
    }

    /**
     * Computes the unitless vibrational entropy of each species.
     */
    template <typename OP>
    void sV(double T, double* const s, const OP& op) {
        int ilevel = 0;
        double fac, sum1, sum2;
        LOOP_MOLECULES(
            op(s[j], 0.0);
        )
    }

    /**
     * Computes the unitless electronic entropy of each species.
     */
    template <typename OP>
    void sE(double T, double* const p_s, const OP& op) {
		op(p_s[0], 0.0);
		for (int i = 0; i < m_ns; i++){
			if (i == 0){
				op(p_s[i],(log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/T)) + (g1_O/g0_O*theta_1_O/T*exp(-theta_1_O/T))/(1+(g1_O/g0_O)*exp(-theta_1_O/T)))); // Ground state
			}
			else {
				op(p_s[i],(log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/T)) + (g1_O2/g0_O2*theta_1_O2/T*exp(-theta_1_O2/T))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/T))));
			}
		}
    }
}; // class STSDB
#undef LOOP
#undef LOOP_HEAVY
#undef LOOP_MOLECULES


// Register the STSDB model with the other thermodynamic databases
Utilities::Config::ObjectProvider<STSDB, ThermoDB> stsDB("STS");

    } // namespace Thermodynamics
} // namespace Mutation


