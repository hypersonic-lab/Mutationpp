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

    STSDB(int arg) : ThermoDB(298.15, 101325.0){
//          ,m_has_electron(false),
//           m_use_tables(false),
//           m_last_bfacs_T(0.0) {} //
//             ~STSDB()
//     {
//         delete [] mp_lnqtmw;
        delete [] mp_hform;
        delete [] mp_indices;
        delete [] mp_rot_data;
// //!        delete [] mp_nvib;
// //!        delete [] mp_vib_temps;

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

        // Setting to zero
    for (int i = 0; i < m_ns; i++){
        cp[i] = 0.;
    }

    // if (cp != NULL && cpt == NULL && cpr == NULL && cpv == NULL && 
    //         cpel == NULL)
    //     {
    //         cpT(cp, Eq());
    //         cpR(cp, PlusEq());
    //         cpV(Tv, cp, PlusEq());
    //         // cpE(Tel, cp, PlusEq());
    //         return;
    //     }

     // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
    //  if (cpt != NULL) {
    //     // cpT(cpt, Eq());  // Setting to 2.5
    //         if (cp != NULL){
    //             // LOOP(cp[i] += cpt[i]);
    //             for (int i = 0; i < m_ns; i++){
    //                 cpt[i] += 2.5; // Cv = 3/2 R; Cp = Cv + R
    //                 cp[i] += 2.5;
    //             }
    //         }
    //  } else {
    //     if (cp != NULL)
    //             //cpT(cp, PlusEq());
    //             // cpT(cp, Eq());   // Setting to 2.5
    //      for (int i = 0; i < m_ns; i++){
    //          cp[i] += 2.5;
    //      }
    //  }

	    if (cpt == NULL) {
            if (cp != NULL)
                //cpT(cp, PlusEq());
                cpT(cp, Eq());
        } else {
            cpT(cpt, Eq());
            if (cp != NULL) 
                LOOP(cp[i] = cpt[i]);
        }

     // Rotation. Assuming fulling active rotational mode
    //  if (cpr != NULL) {
    //      for (int i = 0; i < m_ns; i++){
    //         // cpR(cpr, Eq()); // Setting to 0
    //         if (cp != NULL){
    //             // LOOP_MOLECULES(cp[j] += cpr[j]);
    //          if (i == 0) {
    //             cpr[i] = 0.0; // Ground state
    //             cp[i] += 0.0;
    //              continue; } // Ground state
    //          cpr[i] += 1.0; //iCv = R; Cp = Cv + R
    //          cp[i] += 1.0;}
    //      }

    //  } else {
    //     if (cp != NULL){
    //             // cpR(cp, PlusEq()); // Add cpr to cp
    //      for (int i = 0; i < m_ns; i++){
    //          if (i == 0) {
    //             cp[i] += 0.0; // Ground state
    //              continue; } // Ground state
    //          cp[i] += 1.0;
    //      }}
    //  }

	// Rotation
	if (cpr == NULL) {
		if (cp != NULL)
			cpR(cp, PlusEq());
	} else {
		cpR(cpr, Eq());
		if (cp != NULL) 
			LOOP_MOLECULES(cp[j] += cpr[j]);
	}

     // etc...

     // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
    //  if (cpv != NULL) {
    //     // cpV(Tv, cpv, Eq());
    //         if (cp != NULL) {
    //             // LOOP_MOLECULES(cp[j] += cpv[j]);
    //      for (int i = 0; i < m_ns; i++){
    //         if (i == 0)
    //         {
    //             cpv[i] += 0.0; // Setting as zero for now. Need to think
    //          cp[i] += 0.0;
    //          continue;
    //         }
    //          cpv[i] += 0.0; // Setting as zero for now. Need to think
    //          cp[i] += 0.0;
    //      }
    //         }
    //  } else {
    //     if (cp != NULL){
    //             // cpV(Tv, cp, PlusEq());
    //      for (int i = 0; i < m_ns; i++){
    //         if (i == 0)
    //         {
    //          cp[i] += 0.0;
    //          continue;
    //         }
    //          cp[i] += 0.0;
    //      }
    //  }
    //  }

	        // Vibration
        if (cpv == NULL) {
            if (cp != NULL)
                cpV(Tv, cp, PlusEq());
        } else {
            cpV(Tv, cpv, Eq());
            if (cp != NULL) 
                LOOP_MOLECULES(cp[j] += cpv[j]);
        }


    // double g0_O2 = 3.0;
    // double g1_O2 = 2.0;
    // double theta_1_O2 = 11900;
    // double g0_O = 5.0;
    // double g1_O = 4.0;
    // double theta_1_O = 270;

    //  if (cpel != NULL) {
    //      for (int i = 0; i < m_ns; i++){
    //          if (i == 0) {
    //             // cpel[i] = 0.0; // Ground state
    //             cpel[i] = pow((theta_1_O/Tel),2.0) * (g1_O/g0_O * exp(-theta_1_O / Tel)) / (pow(1.0+g1_O/g0_O * exp(-theta_1_O / Tel),2.0)); // Ground state
    //             cp[i] += cpel[i];
    //              continue; } // Ground state
    //         //  cpel[i] += 0.0; // Boyd p. 110
    //          cpel[i] += pow((theta_1_O2/Tel),2.0) * (g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (pow(1.0+g1_O2/g0_O2 * exp(-theta_1_O2 / Tel),2.0)); // Boyd p. 110
    //          cp[i] += cpel[i];
    //      }

    //  } else {
    //     if (cp != NULL){
    //             // cpR(cp, PlusEq()); // Add cpr to cp
    //      for (int i = 0; i < m_ns; i++){
    //          if (i == 0) {
    //             // cp[i] = 0.0; // Ground state
    //             cp[i] += pow((theta_1_O/Tel),2.0) * (g1_O/g0_O * exp(-theta_1_O / Tel)) / (pow(1.0+g1_O/g0_O * exp(-theta_1_O / Tel),2.0)); // Ground state
    //              continue; } // Ground state
    //         //  cp[i] += 0.0;
    //          cp[i] += pow((theta_1_O2/Tel),2.0) * (g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (pow(1.0+g1_O2/g0_O2 * exp(-theta_1_O2 / Tel),2.0));
    //      }}
    //  }

	// Electronic
	if (cpel == NULL) {
		if (cp != NULL)
			cpE(Tel, cp, PlusEq());
	} else {
		cpE(Tel, cpel, Eq());
		if (cp != NULL)
			LOOP_HEAVY(cp[j] += cpel[j]);
	}

     // Electronic. For now setting as zero
     // cpel[0] = 0.0;
     // cpel[1] = 0.0;
     // cpel[2] = 0.0;

    }
// TODO: Here I am, 10/19/23 --> 4:50pm
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
        double Th, double Te, double Tr, double Tv, double Tel,
        double* const h, double* const ht, double* const hr,
        double* const hv, double* const hel, double* const hf)
    {
        // Given Ts calculate h, ht, hr, hv, hel, hf
        // Old equations, before generalize
    //    h[0] = 2.5 + m_vhf[0];
    //    h[1] = m_vh[1]*Th + m_vhf[1];
    //    h[2] = m_vh[2]*Th + m_vhf[2];

        // const int SIZE = 37;
        // double state[SIZE];
        // double blank[SIZE];
        // double energy[SIZE];
        // string inFileName = "../../data/thermo/oxygen_energy.txt";
        // ifstream inFile;
        // inFile.open(inFileName.c_str());
        // if (inFile.is_open())
        // {
        //     for (int i = 0; i < SIZE; i++)
        //     {
        //         inFile >> state[i];
        //         inFile >> blank[i];
        //         inFile >> energy[i];
        //     }

        //     inFile.close(); // CLose input file
        // }

    // double h_r[50]; 

    // int ctr = 0;
    // for (int i = 0; i < m_nm; ++i){
    //     double val1 = 0.0;
    //     double val2 = 0.0;
    //         for (int j = 0; j < numRoStates[i]; ++j)
    //         {
    //             val1 += (2*j+1)*exp(-1 * ro_energy[ctr][3] / (KB * Tr)) * ro_energy[ctr][3]/(KB * Tr * Tr);
    //             val2 += (2*j+1)*exp(-ro_energy[ctr][3]/(KB*Tr));
    //             ctr += 1;
    //         };
    //     val1 *= Tr * Tr;
    //     double h_val = val1 / val2;;
    //     h_r[i] = h_val;
    // }

        
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

        
        
        // Setting to zero
        for (int i = 0; i < m_ns; i++){
            h[i] = 0.;
        }

        // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        // if (ht != NULL) {
        //     // hT(Th, Te, ht, EqDiv(Th));
        //     // if (h != NULL)
        //     //     LOOP(h[i] = ht[i]);
        //     for (int i = 0; i < m_ns; i++){
        //         ht[i] = 2.5 * Th / Th; // Is this non-dimensional too? Taking in work flow too. Otherwise it would be 1.5.
        //         m_ht[i] = 2.5 * Th / Th;
        //         h[i] = ht[i];
        //     }

        // } else {
        //     //hT(Th, Te, h, Eq());
        //     if (h != NULL){
        //         // hT(Th, Te, h, EqDiv(Th));
        //     for (int i = 0; i < m_ns; i++){
        //         m_ht[i] = 2.5 * Th / Th;
        //         h[i] = 2.5 * Th / Th;
        //     }}
        // }

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

        // Rotation. Assuming fully active rotational mode
        // if (hr != NULL) {
        //     // LOOP(hr[i] = 0.0);
        //     // hR(Tr, hr, EqDiv(Th));
        //     if (h != NULL){
        //         // LOOP_MOLECULES(h[j] += hr[j]);
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             m_hr[i] = 0.0;
        //             hr[i] = 0.0; // Ground state
        //             h[i] += 0.0;
        //             continue; }
        //         // hr[i] = h_r[i-1] / Th;
        //         // m_hr[i] = h_r[i-1] / Th;
        //         // h[i] += h_r[i-1] / Th;
        //         hr[i] = 1.0 * Tr / Th;
        //         m_hr[i] = 1.0 * Tr / Th;
        //         h[i] += 1.0 * Tr / Th;
        //     }}

        // } else {
        //     if (h != NULL){
        //         // hR(Tr, h, PlusEqDiv(Th));
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             m_hr[i] = 0.0;
        //             h[i] += 0.0; // Ground state
        //             continue; }
        //         // m_hr[i] = h_r[i-1] / Th;
        //         // h[i] += h_r[i-1] / Th;
        //         m_hr[i] = 1.0 * Tr / Th;
        //         h[i] += 1.0 * Tr / Th;
        //     }}
        // }

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

        // etc...

        // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
//         if (hv != NULL) {
//             // LOOP(hv[i] = 0.0);
//             // hV(Tv, hv, EqDiv(Th));
//             if (h != NULL){
//                 // LOOP_MOLECULES(h[j] += hv[j]);
//             for (int i = 0; i < m_ns; i++){
//                 if (i == 0) {
//                     hv[i] = 0.0;
//                     h[i] += 0.0; // Ground state
//                     m_hv[i] = 0.0;
//                     continue; }
//                 // exp(-\nu \eps_i / T)
//                 // hv[i] = (i-1) * energy[i-1] * 1.42879 / Th; //* exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
//                 // m_hv[i] = (i-1) * energy[i-1] * 1.42879 / Th;
//                 double cm2J = 1.98630e-23;
//                 hv[i] = m_energy[i-1] * cm2J / (KB * Th); //* exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 m_hv[i] = m_energy[i-1] * cm2J / (KB * Th);
//                 // hv[i] = m_energy[i-1] * 1.42879 / (Th) / (exp(-1.0 * m_energy[i-1]* 1.42879  / Th) - 1.0); //* exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 h[i] += m_energy[i-1] * cm2J / (KB * Th);
//             }}

//         } else {
//             if (h != NULL){
//                 // hV(Tv, h, PlusEqDiv(Th));
//             for (int i = 0; i < m_ns; i++){
//                 if (i == 0) {
//                     h[i] += 0.0;
//                     m_hv[i] = 0.0;
// //                    h[i] = 0.0; // Ground state
//                     continue; }
//                 // exp(-\nu \eps_i / T)
//                 double cm2J = 1.98630e-23;
//                 // m_hv[i] += (i-1) * m_energy[i-1] * 1.42879 / Th;   // * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 m_hv[i] += m_energy[i-1] * cm2J / (KB * Th);   // * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 // m_hv[i] = m_energy[i-1] * 1.42879 / (Th) / (exp(1.0 * m_energy[i-1]* 1.42879  / Th) - 1.0);   // * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 // h[i] += (i-1) * m_energy[i-1] * 1.42879 / Th; //m_energy[i-1] * 1.42879 / Th * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 h[i] += m_energy[i-1] * cm2J / (KB * Th); //m_energy[i-1] * 1.42879 / Th * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//             }}
//         }
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


        // Electronic. For now setting as zero
        // hel[0] = 0.0;
        // hel[1] = 0.0;
        // hel[2] = 0.0;
//         double g0_O2 = 3.0;
//         double g1_O2 = 2.0;
//         double theta_1_O2 = 11900;
//         double g0_O = 5.0;
//         double g1_O = 4.0;
//         double theta_1_O = 270;

//         if (hel != NULL) {
//             // LOOP(hv[i] = 0.0);
//             // hV(Tv, hv, EqDiv(Th));
//             if (h != NULL){
//                 // LOOP_MOLECULES(h[j] += hv[j]);
//             for (int i = 0; i < m_ns; i++){
//                 if (i == 0) {
//                     hel[i] = ((theta_1_O) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel)) / Th;
//                     m_hel[i] = ((theta_1_O) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel)) / Th;
//                     h[i] += ((theta_1_O) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel)) / Th; // Ground state
//                     continue; }
//                 m_hel[i] = ((theta_1_O2) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / Th;   // * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 hel[i] = ((theta_1_O2) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / Th; //* exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 h[i] += ((theta_1_O2) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / Th;
//             }}

//         } else {
//             if (h != NULL){
//                 // hV(Tv, h, PlusEqDiv(Th));
//             for (int i = 0; i < m_ns; i++){
//                 if (i == 0) {
//                     h[i] += ((theta_1_O) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel)) / Th;
//                     m_hel[i] = ((theta_1_O) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel)) / Th;
// //                    h[i] = 0.0; // Ground state
//                     continue; }
//                 m_hel[i] = ((theta_1_O2) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / Th;   // * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//                 h[i] += ((theta_1_O2) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / Th; //m_energy[i-1] * 1.42879 / Th * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
//             }}
//         }

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


        // if (hf != NULL) {
        //     // LOOP(hv[i] = 0.0);
        //     // hV(Tv, hv, EqDiv(Th));
        //     // if (h != NULL)
        //         // LOOP_MOLECULES(h[j] += hv[j]);
        //     for (int i = 0; i < m_ns; i++){
        //         hf[i] = 0.0; //* exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
        //         h[i] += hf[i];
        //     }

        // } else {
        //     if (h != NULL)
        //         // hV(Tv, h, PlusEqDiv(Th));
        //     for (int i = 0; i < m_ns; i++){
        //         h[i] += 0.0; //m_energy[i-1] * 1.42879 / Th * exp(-1*m_energy[i-1] * 1.42879 / Th); // See KMH notes
        //     }
        // }

        // See ln 806 and 655 in RRHO. Adjust enthalpy of formation to any temp
        // double Tss = standardTemperature();
        // // Formation enthalpy
        // if (hf == NULL) {
        //     if (h != NULL){
        //         for (int i = 0; i < m_ns; i++){
        //     if (i == 0){
        //         h[i] += (249229.0 / RU - 2.5*Tss - ((theta_1_O/Tss) * g1_O/g0_O * exp(-theta_1_O / Tss)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tss))) / Th;
		// 		continue;
        //     }
		// 	double cm2J = 1.98630e-23;
        //     h[i] += (0 - (2.5 * Tss + 1.0 * Tss + m_energy[0] * cm2J / (KB) + ((theta_1_O2) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tss)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tss)))) / Th; 
        // }}
        // } else {
        //     for (int i = 0; i < m_ns; i++){
        //     if (i == 0){
        //         hf[i] = (249229.0 / RU - 2.5*Tss - ((theta_1_O/Tss) * g1_O/g0_O * exp(-theta_1_O / Tss)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tss))) / Th;
		// 		continue;
        //     }
		// 	double cm2J = 1.98630e-23;
        //     hf[i] = (0 - (2.5 * Tss + 1.0 * Tss + m_energy[0] * cm2J / (KB) + ((theta_1_O2) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tss)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tss)))) / Th; 
        // }
        //     if (h != NULL){
        //         for (int i = 0; i < m_ns; i++){
        //             h[i] += hf[i];
        //         }
        //         // LOOP(h[i] += hf[i]);
        //     }
        // }
        // std::cout << "Method 1: " << hf[0] << " " << hf[1] << endl;

        
        // Formation enthalpy
        if (hf == NULL) {
            if (h != NULL)
                hF(h, PlusEqDiv(Th));
        } else {
            hF(hf, EqDiv(Th));
            if (h != NULL)
                LOOP(h[i] += hf[i]);
        }
        // std::cout << "Method 2: " << hf[0] << " " << hf[1] << endl;



        // Old equations, before generalize
        // I think this below is already calculated above by the sums
        // h[0] = ht[0] + hr[0] + hv[0] + hel[0] + m_vhf[0];
        // h[1] = ht[1] + hr[1] + hv[1] + hel[1] + m_vhf[1];
        // h[2] = ht[2] + hr[2] + hv[2] + hel[2] + m_vhf[2];
        //h[0] = 0.0; // m_vh[0]*Th + m_vhf[0];
        //h[1] = 1.0; // m_vh[1]*Th + m_vhf[1];
        //h[2] = 1.0; // m_vh[2]*Th + m_vhf[2];

      //  hf[0] = 1.;
      //  hf[1] = 2.;
      //  hf[2] = 3.;
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
        double* const sel)//, double* const hr)
    {
        // Given Ts calculate s, st, sr, sv, sel
        // Note: Check if NULL

            // Following similar approach as enthalpy
            // Setting to zero

        //  if (st == NULL && sr == NULL && sv == NULL && sel == NULL) {
        //     sT(Th, Te, P, s, Eq());
        //     sR(Tr, s, PlusEq());
        //     sV(Tv, s, PlusEq());
        //     // sE(Tel, s, PlusEq());

        //     // Include spin contribution for free electron entropy
        //     if (m_has_electron)
        //         s[0] += std::log(2.0);

        //     return;
        // }
    // double s_r[50]; 

    // int ctr = 0;
    // for (int i = 0; i < m_nm; ++i){
    //     double val = 0.0;
    //     double val2 = 0.0;
    //         for (int j = 0; j < numRoStates[i]; ++j)
    //         {
    //             val += (2*j+1)*exp(-1 * ro_energy[ctr][3] / (KB * Tr));
    //             val2 += (2*j+1)*exp(-ro_energy[ctr][3]/(KB*Tr)) * ro_energy[ctr][3]/(KB*Tr*Tr);
    //             ctr += 1;
    //         };
    //     double val1 = log(0.5 * val);
    //     val2 *= Tr;
    //     double val3 = val;
    //     double s_val = val1 + val2 / val3;
    //     s_r[i] = s_val;
    // }

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


        for (int i = 0; i < m_ns; i++){
            s[i] = 0.;
        }
        
        // enthalpy(Th, Te, Tr, Tv, Tel, s, NULL, NULL, NULL, NULL, NULL);
        // gibbs(Th, Te, Tr, Tv, Tel, P, s, NULL, NULL, NULL, NULL);

        // double test[50];

        // enthalpy(Th, Te, Tr, Tv, Tel, test, NULL, NULL, NULL, NULL, NULL);
        // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        // Will need to upload masses of each species

        
        // if (st != NULL) {
        //     // sT(Th, Te, P, st, Eq());
        //     // LOOP(s[i] = st[i]);
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // EQ 3.90 of Boyd. // Ground state
        //             m_st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //             s[i] +=  2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //             continue; }
        //         st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //         s[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //         m_st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //    }

        // } else {
        //     // sT(Th, Te, P, s, Eq());
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             s[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // EQ 3.90 of Boyd. // Ground state
        //             m_st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //             // s[i] += st[i];
        //             continue; }
        //         s[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //         m_st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //     }
        // }

        if (st == NULL)
            sT(Th, Te, P, s, Eq());
        else {
            sT(Th, Te, P, st, Eq());
            LOOP(s[i] = st[i]);
        }

        // Rotation. Assuming fulling active rotational mode
        // double ThetaR = 2.08; //char temp rot O2
        // if (sr != NULL) {
        //     // LOOP(sr[i] = 0.0);
        //     // sR(Tr, sr, Eq());
        //     // LOOP_MOLECULES(s[j] += sr[j]);
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             sr[i] += 0.0; // Ground state
        //             s[i] += 0.0;
        //             m_sr[i] += 0.0;
        //             continue; }
        //         // sr[i] += s_r[i-1]; // From slide 20 of Magin, need to check units
        //         // s[i] += s_r[i-1];
        //         // m_sr[i] += s_r[i-1]; // From Magin above;
		// 		// ln(omegaT) = log(thetaR) + 2/Lin * log*(steric)
		// 		// :::: Linearity * (1 + ln(T) - log(thetaR) + 2/Lin * log*(steric))
		// 		// :::: 1 + ln(T) - log(thetaR) - 2/Lin * log*(steric)
		// 		// :::: 1 + ln(T) - log(thetaR) - log(2)
		// 		// :::: 1 + ln(T/(2thetaR))
        //         sr[i] += (log(Tr / (2 * ThetaR)) + 1.0); // From slide 20 of Magin, need to check units
        //         s[i] += (log(Tr / (2 * ThetaR)) + 1.0);
        //         m_sr[i] += (log(Tr / (2 * ThetaR)) + 1.0); // From Magin above;

        //     }
            
        //     // Old equations, before generalize
        //     // sr[1] = 1.0 + log((0.5 * Th / 2.1) / N ) + 1.0; // Eq. 3.78 of Boyd. Need to define N or substitute
        //     // sr[2] = 1.0 + log((0.5 * Th / 2.1) / N ) + 1.0; // Eq. 3.78 of Boyd. Need to define N or substitute

        // } else {
        //     // sR(Tr, s, PlusEq());
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             s[i] += 0.0;
        //             m_sr[i] += 0.0;
        //             continue; }
        //         // s[i] += s_r[i-1]; // From Magin above;
        //         // m_sr[i] += s_r[i-1]; // From Magin above;
        //         s[i] += (log(Tr / (2 * ThetaR)) + 1.0); // From Magin above;
        //         m_sr[i] += (log(Tr / (2 * ThetaR)) + 1.0); // From Magin above;
        //     }

        // }
        
        // Rotational entropy
        if (sr == NULL)
            sR(Tr, s, PlusEq());
        else {
            LOOP(sr[i] = 0.0);
            sR(Tr, sr, Eq());
            LOOP_MOLECULES(s[j] += sr[j]);
        }
        // etc...

        // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
        // if (sv != NULL) {
        //     // LOOP(sv[i] = 0.0);
        //     // sV(Tv, sv, Eq());
        //     // LOOP_MOLECULES(s[j] += sv[j]);
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             sv[i] += 0.0; // Ground state
        //             m_sv[i] += 0.0;
        //             s[i] += 0.0;
        //          continue; } // Ground state
        //         sv[i] += 0.0;//m_hv[i] / Th - 1/Th * log(1-exp(-1.0*m_energy[i-1] * 1.42879/Th)); // Setting to 0 based on discussion with George -- no degeneracy, don't lose any info since sts
        //         s[i] += 0.0;
        //     }
        //     // Old equations, before generalize
        //     // sv[1] = 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th; // Eq. 3.78 of Boyd. Need to define N or substitute
        //     // sv[2] =  1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        // } else {
        //     // sV(Tv, s, PlusEq());
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             s[i] += 0.0; // Ground state
        //             m_sv[i] += 0.0;
        //          continue; } // Ground state
        //         s[i] += 0.0;//m_hv[i] / Th - 1/Th * log(1-exp(-1.0*m_energy[i-1] * 1.42879/Th));
        //         m_sv[i] += 0.0;//m_hv[i] / Th - 1/Th * log(1-exp(-1.0*m_energy[i-1] * 1.42879/Th));
        //     }
        //     // Old equations, before generalize
        //     // s[1] += 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th;
        //     // s[2] += 1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        // }

        // Vibrational entropy
        if (sv == NULL)
            sV(Tv, s, PlusEq());
        else {
            LOOP(sv[i] = 0.0);
            sV(Tv, sv, Eq());
            LOOP_MOLECULES(s[j] += sv[j]);
        }

        // double g0_O2 = 3.0;
        // double g1_O2 = 2.0;
        // double theta_1_O2 = 11900;
        // double g0_O = 5.0;
        // double g1_O = 4.0;
        // double theta_1_O = 270;


        // if (sel != NULL) {
        //     // LOOP(sv[i] = 0.0);
        //     // sV(Tv, sv, Eq());
        //     // LOOP_MOLECULES(s[j] += sv[j]);
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             sel[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel))); // Ground state
        //             m_sel[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel)));
        //             s[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel)));
        //          continue; } // Ground state
        //         sel[i] += (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel))); // Setting to 0 based on discussion with George -- no degeneracy, don't lose any info since sts
        //         m_sel[i] +=  (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
        //         s[i] += (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
        //     }
        //     // Old equations, before generalize
        //     // sv[1] = 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th; // Eq. 3.78 of Boyd. Need to define N or substitute
        //     // sv[2] =  1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        // } else {
        //     // sV(Tv, s, PlusEq());
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             s[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel))); // Ground state
        //             m_sel[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel)));
        //          continue; } // Ground state
        //         s[i] += (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
        //         m_sel[i] += (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
        //     }
        //     // Old equations, before generalize
        //     // s[1] += 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th;
        //     // s[2] += 1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        // }

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




        // Electronic. For now setting as zero
        // sel[0] = 0.0;
        // sel[1] = 0.0;
        // sel[2] = 0.0;

        //h[0] += m_vhf[0];
        //h[1] += m_vhf[1];
        // h[2] += m_vhf[2];
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
        // enthalpy(Th, Te, Tr, Tv, Tel, g, NULL, NULL, NULL, NULL, NULL);

        // Subtract the entropies
        // sT(Th, Te, P, g, MinusEq());
        // sR(Tr, g, MinusEq());
        // sV(Tv, g, MinusEq());
        // sE(Tel, g, MinusEq());

        // Account for spin of free electrons
        // if (m_has_electron)
            // g[0] -= std::log(2.0);
    // }

        // // Given Ts calculate g, gt, gr, gv, gel
        // // Note: Check if NULL

        //            // Following similar approach as enthalpy
        //     // Setting to zero

        for (int i = 0; i < m_ns; i++){
            g[i] = 0.;
        }

        // for (int i = 0; i < m_ns; i++){
        //     g[i] = 0.;
        // }
        enthalpy(Th, Te, Tr, Tv, Tel, g, NULL, NULL, NULL, NULL, NULL);

        // Subtract the entropies
        sT(Th, Te, P, g, MinusEq());
        sR(Tr, g, MinusEq());
        sV(Tv, g, MinusEq());
        sE(Tel, g, MinusEq());

        // for (int i = 0; i < m_ns; i++){
        //     g[i] = 0.;
        // }
        // entropy(Th, Te, Tr, Tv, Tel, P, g, NULL, NULL, NULL, NULL);
        // for (int i = 0; i < m_ns; i++){
        //     g[i] = 0.;
        // }
        // // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        // // Will need to upload masses of each species
        // if (gt != NULL) {
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             gt[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //             continue;
        //         }
        //         gt[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // G = H - TS
        //         g[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //     }
        // } else {
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0){
        //             g[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
        //             continue;
        //         }
        //         g[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // G = H - TS
        //     }
        // }

        // double ThetaR = 2.08; //char temp rot O2
        // if (gr != NULL) {
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             gr[i] -= 0.0; // Ground state
        //             g[i] -= 0.0;
        //             continue; }
        //         gr[i] -= (log(Tr / (2 * ThetaR)) + 1.0); // G = H - TS
        //         g[i] -= (log(Tr / (2 * ThetaR)) + 1.0);
        //     }
        // } else {
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             g[i] -= 0.0;
        //             continue; }
        //         g[i] -= (log(Tr / (2 * ThetaR)) + 1.0);  // G = H - TS
        //     }
        // }

        // if (gv != NULL) {
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             gv[i] -= 0.0; // Ground state
        //             g[i] -= 0.0;
        //          continue; } // Ground state
        //         gv[i] -= 0.0; // G = H - TS // Tv?
        //         g[i] -= 0.0;
        //     }
        // } else {
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             g[i] -= 0.0;
        //          continue; } // Ground state
        //         g[i] -= 0.0;
        //     }
        // }
        
        // double g0_O2 = 3.0;
        // double g1_O2 = 2.0;
        // double theta_1_O2 = 11900;
        // double g0_O = 5.0;
        // double g1_O = 4.0;
        // double theta_1_O = 270;
        // if (gel != NULL) {
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             gel[i] -= (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel))); // Ground state
        //             g[i] -= (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel)));
        //          continue; } // Ground state
        //         gel[i] -= (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel))); // Setting to 0 based on discussion with George -- no degeneracy, don't lose any info since sts
        //         g[i] -= (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
        //     }
        // } else {
        //     for (int i = 0; i < m_ns; i++){
        //         if (i == 0) {
        //             g[i] -= (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel))); // Ground state
        //          continue; } // Ground state
        //         g[i] -= (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
        //     }
        // }
        // // Electronic. For now setting as zero
    //     // sel[0] = 0.0;
    //     // sel[1] = 0.0;
    //     // sel[2] = 0.0;

    //     //h[0] += m_vhf[0];
    //     //h[1] += m_vhf[1];
    //     // h[2] += m_vhf[2];
    }
    
//Calculate enthalpy from Gibbs

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
            for (size_t i = 0; i < rrho.nElectronicLevels(); ++i)
                species.push_back(Species(ground_state, i));
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
//        m_ns = 3; // Number of energy states

        m_ns = species().size();
		m_has_electron = (species()[0].type() == ELECTRON);

		        // Compute the contribution of the partition functions at the standard
        // state temperature to the species enthalpies

        // Store the species formation enthalpies in K

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


        mp_part_sst = new double [m_ns];
        hT(Tss, Tss, mp_part_sst, Eq());
        hR(Tss, mp_part_sst, PlusEq());
        hV_f(Tss, mp_part_sst, PlusEq());
        hE(Tss, mp_part_sst, PlusEq());










        m_vh.resize(m_ns);
        m_vhf.resize(50);
        m_ht.resize(m_ns);
        m_hr.resize(m_ns);
        m_hv.resize(m_ns);
        m_hel.resize(m_ns);
        m_hf.resize(m_ns);
        m_st.resize(m_ns);
        m_sr.resize(m_ns);
        m_sv.resize(m_ns);
        m_sel.resize(m_ns);




        // Add ht, hr, hv...

        // Nitrogen m_vhf
//        m_vhf[0] = 472440;
//        m_vhf[1] = 19425.13;
//        m_vhf[2] = 71231.59;
//        m_vhf[3] = 143937.47;
//        m_vhf[4] = 233797.83;
//        m_vhf[5] = 338000.97;
//        m_vhf[6] = 454210.59;
//        m_vhf[7] = 582378.25;
//        m_vhf[8] = 721541.19;
//        m_vhf[9] = 873713.56;

        m_vh[0] = 0.; // Atomic oxygen
        m_vh[1] = 7.87380953594E+02;
        m_vh[2] = 1.4E+03;
        // m_vh[3] = 7.87380953594E+02;

        //USE THESE
        // I think these are enthlapy of formation
        // m_vhf[0] = 100.59; // Atomic oxygen
        // m_vhf[1] = 7.87380953594E+02;
        // m_vhf[2] = 1.0+02;
        // m_vh[3] = 7.87380953594E+02;

        // @todo: 1/18/2023
        // Load the necessary thermodynamic data
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

	double g0_O2 = 3.0;
	double g1_O2 = 2.0;
	double theta_1_O2 = 11900;
	double g0_O = 5.0;
	double g1_O = 4.0;
	double theta_1_O = 270;
	double cm2J = 1.98630e-23;
    double Tss = standardTemperature();
    double ThetaR = 2.08; //char temp rot O2



	std::array<double, 47> m_energy = {
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


    // Store here only the necessary data for calculating species thermodynamics
    // const int m_ns = 48; // need to see how to recognize number of states from M++
    // const int m_na = 1; // need to see how to recognize number of states from M++
    // const int m_nm = 47; // need to see how to recognize number of states from M++
    // double m_vh[m_ns];
    // double m_vhf[m_ns];
    // double hv[m_ns];
    // double ht[m_ns];
    // double hr[m_ns];
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
            LOOP_MOLECULES(op(h[j], m_energy[i] * cm2J / (KB)))
        }
    }

    template <typename OP>
    void hV_f(double T, double* const h, const OP& op) {
        if (T < 10.0) {
            LOOP_MOLECULES(op(h[j], 0.0));
        } else {
            LOOP_MOLECULES(op(h[j], m_energy[0] * cm2J / (KB)))
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
        // const double onelnT = 1.0 + std::log(T);
        LOOP_MOLECULES(
            // op(s[j], mp_rot_data[i].linearity * (onelnT -
            //     mp_rot_data[i].ln_omega_t));
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
//!            sum1 = sum2 = 0.0;
//!            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
//!                fac  =  std::exp(mp_vib_temps[ilevel] / T);
//!                sum1 += mp_vib_temps[ilevel] / (fac - 1.0);
//!                sum2 += std::log(1.0 - 1.0 / fac);
//!            }
//!            op(s[j], (sum1 / T - sum2));
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


        // updateElecBoltzmannFactors(T);
        // op(p_s[0], 0.0);

        // double* facs = mp_el_bfacs;
        // for (int i = 0; i < m_elec_data.nheavy; ++i, facs += 3) {
        //     if (facs[0] > 0)
        //         op(p_s[i+m_elec_data.offset],
        //             (facs[1]/(facs[0]*T) + std::log(facs[0])));
        //     else
        //         op(p_s[i+m_elec_data.offset], 0.0);
        // }
    }


    // double sv[m_ns];
    // double st[m_ns];
    // double sr[m_ns];
    std::vector<double> m_vh {};
    std::vector<double> m_vhf {};
    std::vector<double> m_hv {}; //should this be in private?? //m_ for private
    std::vector<double> m_ht {}; //should this be in private??
    std::vector<double> m_hr {}; //should this be in private??
    std::vector<double> m_hel {}; //should this be in private??
    std::vector<double> m_hf {}; //should this be in private??
    std::vector<double> m_sv {}; //should this be in private??
    std::vector<double> m_st {}; //should this be in private??
    std::vector<double> m_sr {}; //should this be in private??
    std::vector<double> m_sel {}; //should this be in private??

}; // class STSDB
#undef LOOP
#undef LOOP_HEAVY
#undef LOOP_MOLECULES


// Register the STSDB model with the other thermodynamic databases
Utilities::Config::ObjectProvider<STSDB, ThermoDB> stsDB("STS");

    } // namespace Thermodynamics
} // namespace Mutation


