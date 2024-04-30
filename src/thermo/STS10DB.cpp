/**
 * @file STS10DB.cpp
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
    j = i;\
    __op__ ;\
} // j = mp_indices[i]

// Loops over molecules. Inside the loop, the index i is zero based and indexes
// internal molecular data.  The index j is the index corresponding to the 
// original species data.
#define LOOP_MOLECULES(__op__)\
for (int i = 0, j = 0; i < m_nm; ++i) {\
    j = m_na+i;\
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

// typedef struct {
//     double ln_omega_t;  // ln(omega^(2 / L) * theta)
//     double linearity;   // L / 2
// } RotData;

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
class STS10DB : public ThermoDB
{
public:

    STS10DB(int arg) : ThermoDB(298.15, 101325.0){}
//          ,m_has_electron(false),
//           m_use_tables(false),
//           m_last_bfacs_T(0.0) {} //
//             ~STS10DB()
//     {
//         delete [] mp_lnqtmw;
//         delete [] mp_hform;
//         delete [] mp_indices;
// //!        delete [] mp_rot_data;
// //!        delete [] mp_nvib;
// //!        delete [] mp_vib_temps;

//         delete [] m_elec_data.p_nelec;
//         delete [] m_elec_data.p_levels;
//         delete [] mp_part_sst;
//         delete [] mp_el_bfacs;

//         if (m_use_tables) {
//             delete mp_el_bfac_table;
//         }
//     }

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
     if (cpt != NULL) {
        // cpT(cpt, Eq());  // Setting to 2.5
            if (cp != NULL){
                // LOOP(cp[i] += cpt[i]);
                for (int i = 0; i < m_ns; i++){
                    cpt[i] += 2.5; // Cv = 3/2 R; Cp = Cv + R
                    cp[i] += cpt[i];
                }
            }
     } else {
        if (cp != NULL)
                //cpT(cp, PlusEq());
                // cpT(cp, Eq());   // Setting to 2.5
         for (int i = 0; i < m_ns; i++){
             cp[i] += 2.5;
         }
     }

     // Rotation. Assuming fulling active rotational mode
     if (cpr != NULL) {
         for (int i = 0; i < m_ns; i++){
            // cpR(cpr, Eq()); // Setting to 0
            if (cp != NULL){
                // LOOP_MOLECULES(cp[j] += cpr[j]);
             if (i == 0) {
                cpr[i] = 0.0; // Ground state
                cp[i] += cpr[i];
                 continue; } // Ground state
             cpr[i] += 2.0; //iCv = R; Cp = Cv + R
             cp[i] += cpr[i];}
         }

     } else {
        if (cp != NULL){
                // cpR(cp, PlusEq()); // Add cpr to cp
         for (int i = 0; i < m_ns; i++){
             if (i == 0) {
                cp[i] += 0.0; // Ground state
                 continue; } // Ground state
             cp[i] += 2.0;
         }}
     }

     // etc...

     // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
     if (cpv != NULL) {
        // cpV(Tv, cpv, Eq());
            if (cp != NULL) {
                // LOOP_MOLECULES(cp[j] += cpv[j]);
         for (int i = 0; i < m_ns; i++){
            if (i == 0)
            {
                cpv[i] += 0.0; // Setting as zero for now. Need to think
             cp[i] += cpv[i];
             continue;
            }
             cpv[i] += 1.0; // Setting as zero for now. Need to think
             cp[i] += cpv[i];
         }
            }
     } else {
        if (cp != NULL){
                // cpV(Tv, cp, PlusEq());
         for (int i = 0; i < m_ns; i++){
            if (i == 0)
            {
                cpv[i] += 0.0; // Setting as zero for now. Need to think
             cp[i] += cpv[i];
             continue;
            }
             cp[i] += 1.0;
         }
     }
     }

    // Rotation. Assuming fulling active rotational mode
    double g0_O2 = 3.0;
    double g1_O2 = 2.0;
    double theta_1_O2 = 11900;
    double g0_O = 5.0;
    double g1_O = 4.0;
    double theta_1_O = 270;

     if (cpel != NULL) {
         for (int i = 0; i < m_ns; i++){
             if (i == 0) {
                // cpel[i] = 0.0; // Ground state
                cpel[i] = 1.0/Th * pow((theta_1_O/Tel),2.0) * (g1_O/g0_O * exp(-theta_1_O / Tel)) / (pow(1.0+g1_O/g0_O * exp(-theta_1_O / Tel),2.0)) + 1.0; // Ground state
                cp[i] += cpel[i];
                 continue; } // Ground state
            //  cpel[i] += 0.0; // Boyd p. 110
             cpel[i] += 1.0/Th * pow((theta_1_O2/Tel),2.0) * (g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (pow(1.0+g1_O2/g0_O2 * exp(-theta_1_O2 / Tel),2.0)) + 1.0; // Boyd p. 110
             cp[i] += cpel[i];
         }

     } else {
        if (cp != NULL){
                // cpR(cp, PlusEq()); // Add cpr to cp
         for (int i = 0; i < m_ns; i++){
             if (i == 0) {
                // cp[i] = 0.0; // Ground state
                cp[i] += 1.0/Th * pow((theta_1_O/Tel),2.0) * (g1_O/g0_O * exp(-theta_1_O / Tel)) / (pow(1.0+g1_O/g0_O * exp(-theta_1_O / Tel),2.0)) + 1.0; // Ground state
                 continue; } // Ground state
            //  cp[i] += 0.0;
             cp[i] += 1.0/Th * pow((theta_1_O2/Tel),2.0) * (g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (pow(1.0+g1_O2/g0_O2 * exp(-theta_1_O2 / Tel),2.0)) + 1.0;
         }}
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
        double energy[50];
        energy[0] = 786.0234;
        energy[1] = 2343.573;
        energy[2] = 3881.3038;
        energy[3] = 5398.5964;
        energy[4] = 6894.8782;
        energy[5] = 8369.5118;
        energy[6] = 9822.0457;
        energy[7] = 11251.754;
        energy[8] = 12658.072;
        energy[9] = 14040.4352;
        energy[10] = 15398.3596;
        energy[11] = 16731.2;
        energy[12] = 18038.3111;
        energy[13] = 19319.2091;
        energy[14] = 20573.1679;
        energy[15] = 21799.7037;
        energy[16] = 22998.2518;
        energy[17] = 24168.0058;
        energy[18] = 25308.4816;
        energy[19] = 26419.0341;
        energy[20] = 27498.9373;
        energy[21] = 28547.546;
        energy[22] = 29564.2148;
        energy[23] = 30548.1374;
        energy[24] = 31498.6683;
        energy[25] = 32414.9205;
        energy[26] = 33296.168;
        energy[27] = 34141.6041;
        energy[28] = 34950.3419;
        energy[29] = 35721.3326;
        energy[30] = 36453.8504;
        energy[31] = 37146.6854;
        energy[32] = 37798.8698;
        energy[33] = 38409.2744;
        energy[34] = 38976.5281;
        energy[35] = 39499.5822;
        energy[36] = 39976.9044;
        energy[37] = 40407.0429;
        energy[38] = 40788.6264;
        energy[39] = 41119.9613;
        energy[40] = 41399.6763;
        energy[41] = 41626.4003;
        energy[42] = 41799.6494;
        energy[43] = 41920.3914;
        energy[44] = 41993.5464;
        energy[45] = 42029.6803;
        energy[46] = 42042.9885;
//        double energy[11];
//        energy[0] = 787.380953594;
//        energy[1] = 2343.76026609;
//        energy[2] = 3876.56829159;
//        energy[3] = 5386.02874609;
//        energy[4] = 6872.33479359;
//        energy[5] = 8335.64904609;
//        energy[6] = 9776.10356359;
//        energy[7] = 11193.7998541;
//        energy[8] = 12588.8088736;
//        energy[9] = 13961.1710261;
        // energy[10] = 15310.8961636;
        // energy[11] = 16637.9635861;
        // energy[12] = 17942.3220416;
        // energy[13] = 19223.8897261;
        // energy[14] = 20482.5542836;
        // energy[15] = 21718.1728061;
        // energy[16] = 22930.5718336;
        // energy[17] = 24119.5473541;
        // energy[18] = 25284.8648036;
        // energy[19] = 26426.2590661;
        // energy[20] = 27543.4344736;
        // energy[21] = 28636.0648061;
        // energy[22] = 29703.7932916;
        // energy[23] = 30746.2326061;
        // energy[24] = 31762.9648736;
        // energy[25] = 32753.5416661;
        // energy[26] = 33717.4840036;
        // energy[27] = 34654.2823541;
        // energy[28] = 35563.3966336;
        // energy[29] = 36444.2562061;
        // energy[30] = 37296.2598836;
        // energy[31] = 38118.7759261;
        // energy[32] = 38911.1420416;
        // energy[33] = 39672.6653861;
        // energy[34] = 40402.6225636;
        // energy[35] = 41100.2596261;
        // energy[36] = 41764.7920736;
        // energy[2] = 3876.56829159;
        
                // Special case where we only want the total enthalpy
        // if (ht == NULL && hr == NULL && hv == NULL && hel == NULL && 
        //     hf == NULL && h != NULL) 
        // {
        //     hT(Th, Te, h, Eq());
        //     hR(Tr, h, PlusEq());
        //     hV(Tv, h, PlusEq());
        //     // hE(Tel, h, PlusEq());
        //     // hF(h, PlusEq());
        //     LOOP(h[i] /= Th);
        //     return;
        // }

        
        
        // Setting to zero
        for (int i = 0; i < m_ns; i++){
            h[i] = 0.;
        }

        // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        if (ht != NULL) {
            // hT(Th, Te, ht, EqDiv(Th));
            // if (h != NULL)
            //     LOOP(h[i] = ht[i]);
            for (int i = 0; i < m_ns; i++){
                ht[i] += 2.5; // Is this non-dimensional too? Taking in work flow too. Otherwise it would be 1.5.
                m_ht[i] += 2.5;
                h[i] += ht[i];
            }

        } else {
            //hT(Th, Te, h, Eq());
            if (h != NULL){
                // hT(Th, Te, h, EqDiv(Th));
            for (int i = 0; i < m_ns; i++){
                m_ht[i] += 2.5;
                h[i] += 2.5;
            }}
        }

        // Rotation. Assuming fully active rotational mode
        if (hr != NULL) {
            // LOOP(hr[i] = 0.0);
            // hR(Tr, hr, EqDiv(Th));
            if (h != NULL){
                // LOOP_MOLECULES(h[j] += hr[j]);
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    m_hr[i] += 0.0;
                    hr[i] += 0.0; // Ground state
                    h[i] += 0.0;
                    continue; }
                hr[i] += 1.0;
                m_hr[i] += 1.0;
                h[i] += 1.0;
            }}

        } else {
            if (h != NULL){
                // hR(Tr, h, PlusEqDiv(Th));
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    m_hr[i] += 0.0;
                    h[i] += 0.0; // Ground state
                    continue; }
                m_hr[i] += 1.0;
                h[i] += 1.0;
            }}
        }

        // etc...

        // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
        if (hv != NULL) {
            // LOOP(hv[i] = 0.0);
            // hV(Tv, hv, EqDiv(Th));
            if (h != NULL){
                // LOOP_MOLECULES(h[j] += hv[j]);
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    hv[i] += 0.0;
                    h[i] += 0.0; // Ground state
                    m_hv[i] += 0.0;
                    continue; }
                hv[i] += (i-1) * energy[i-1] * 1.42879 / (Tv); //* exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
                m_hv[i] += (i-1) * energy[i-1] * 1.42879 / (Tv);
                // hv[i] = energy[i-1] * 1.42879 / (Th) / (exp(-1.0 * energy[i-1]* 1.42879  / Th) - 1.0); //* exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
                h[i] += (i-1) * energy[i-1] * 1.42879 / (Tv);
            }}

        } else {
            if (h != NULL){
                // hV(Tv, h, PlusEqDiv(Th));
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    h[i] += 0.0;
                    m_hv[i] += 0.0;
//                    h[i] = 0.0; // Ground state
                    continue; }
                m_hv[i] += (i-1) * energy[i-1] * 1.42879 / (Tv);   // * exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
                // m_hv[i] = energy[i-1] * 1.42879 / (Th) / (exp(1.0 * energy[i-1]* 1.42879  / Th) - 1.0);   // * exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
                h[i] += (i-1) * energy[i-1] * 1.42879 / (Tv); //energy[i-1] * 1.42879 / Th * exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
            }}
        }

        // Electronic. For now setting as zero
        // hel[0] = 0.0;
        // hel[1] = 0.0;
        // hel[2] = 0.0;
        double g0_O2 = 3.0;
        double g1_O2 = 2.0;
        double theta_1_O2 = 11900;
        double g0_O = 5.0;
        double g1_O = 4.0;
        double theta_1_O = 270;

        if (hel != NULL) {
            // LOOP(hv[i] = 0.0);
            // hV(Tv, hv, EqDiv(Th));
            if (h != NULL){
                // LOOP_MOLECULES(h[j] += hv[j]);
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    hel[i] += ((theta_1_O/Tel) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel));
                    m_hel[i] += ((theta_1_O/Tel) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel));
                    h[i] += ((theta_1_O/Tel) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel)); // Ground state
                    continue; }
                m_hel[i] += ((theta_1_O2/Tel) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel));   // * exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
                hel[i] += ((theta_1_O2/Tel) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)); //* exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
                h[i] += ((theta_1_O2/Tel) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel));
            }}

        } else {
            if (h != NULL){
                // hV(Tv, h, PlusEqDiv(Th));
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    h[i] += ((theta_1_O/Tel) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel));
                    m_hel[i] += ((theta_1_O/Tel) * g1_O/g0_O * exp(-theta_1_O / Tel)) / (1.0 + g1_O/g0_O * exp(-theta_1_O / Tel));
//                    h[i] = 0.0; // Ground state
                    continue; }
                m_hel[i] += ((theta_1_O2/Tel) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel));   // * exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
                h[i] += ((theta_1_O2/Tel) * g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)) / (1.0 + g1_O2/g0_O2 * exp(-theta_1_O2 / Tel)); //energy[i-1] * 1.42879 / Th * exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
            }}
        }

        // if (hf != NULL) {
        //     // LOOP(hv[i] = 0.0);
        //     // hV(Tv, hv, EqDiv(Th));
        //     // if (h != NULL)
        //         // LOOP_MOLECULES(h[j] += hv[j]);
        //     for (int i = 0; i < m_ns; i++){
        //         hf[i] = 0.0; //* exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
        //         h[i] += hf[i];
        //     }

        // } else {
        //     if (h != NULL)
        //         // hV(Tv, h, PlusEqDiv(Th));
        //     for (int i = 0; i < m_ns; i++){
        //         h[i] += 0.0; //energy[i-1] * 1.42879 / Th * exp(-1*energy[i-1] * 1.42879 / Th); // See KMH notes
        //     }
        // }

//         for (int i = 0; i < m_ns; i++){
// //            h[i] += energy[i];
//             h[i] += m_vhf[i];
//         }

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
        double energy[50];
        energy[0] = 786.0234;
        energy[1] = 2343.573;
        energy[2] = 3881.3038;
        energy[3] = 5398.5964;
        energy[4] = 6894.8782;
        energy[5] = 8369.5118;
        energy[6] = 9822.0457;
        energy[7] = 11251.754;
        energy[8] = 12658.072;
        energy[9] = 14040.4352;
        energy[10] = 15398.3596;
        energy[11] = 16731.2;
        energy[12] = 18038.3111;
        energy[13] = 19319.2091;
        energy[14] = 20573.1679;
        energy[15] = 21799.7037;
        energy[16] = 22998.2518;
        energy[17] = 24168.0058;
        energy[18] = 25308.4816;
        energy[19] = 26419.0341;
        energy[20] = 27498.9373;
        energy[21] = 28547.546;
        energy[22] = 29564.2148;
        energy[23] = 30548.1374;
        energy[24] = 31498.6683;
        energy[25] = 32414.9205;
        energy[26] = 33296.168;
        energy[27] = 34141.6041;
        energy[28] = 34950.3419;
        energy[29] = 35721.3326;
        energy[30] = 36453.8504;
        energy[31] = 37146.6854;
        energy[32] = 37798.8698;
        energy[33] = 38409.2744;
        energy[34] = 38976.5281;
        energy[35] = 39499.5822;
        energy[36] = 39976.9044;
        energy[37] = 40407.0429;
        energy[38] = 40788.6264;
        energy[39] = 41119.9613;
        energy[40] = 41399.6763;
        energy[41] = 41626.4003;
        energy[42] = 41799.6494;
        energy[43] = 41920.3914;
        energy[44] = 41993.5464;
        energy[45] = 42029.6803;
        energy[46] = 42042.9885;

        for (int i = 0; i < m_ns; i++){
            s[i] = 0.;
        }
        
        // enthalpy(Th, Te, Tr, Tv, Tel, s, NULL, NULL, NULL, NULL, NULL);
        // gibbs(Th, Te, Tr, Tv, Tel, P, s, NULL, NULL, NULL, NULL);

        // double test[50];

        // enthalpy(Th, Te, Tr, Tv, Tel, test, NULL, NULL, NULL, NULL, NULL);
        // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        // Will need to upload masses of each species
        if (st != NULL) {
            // sT(Th, Te, P, st, Eq());
            // LOOP(s[i] = st[i]);
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // EQ 3.90 of Boyd. // Ground state
                    m_st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                    s[i] +=  2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                    continue; }
                st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                s[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                m_st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
           }

        } else {
            // sT(Th, Te, P, s, Eq());
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    s[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // EQ 3.90 of Boyd. // Ground state
                    m_st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                    // s[i] += st[i];
                    continue; }
                s[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                m_st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
            }
        }

        // Rotation. Assuming fulling active rotational mode
        double ThetaR = 2.08; //char temp rot O2
        if (sr != NULL) {
            // LOOP(sr[i] = 0.0);
            // sR(Tr, sr, Eq());
            // LOOP_MOLECULES(s[j] += sr[j]);
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    sr[i] += 0.0; // Ground state
                    s[i] += 0.0;
                    m_sr[i] += 0.0;
                    continue; }
                sr[i] += (log(Tr / (2 * ThetaR)) + 1.0); // From slide 20 of Magin, need to check units
                s[i] += (log(Tr / (2 * ThetaR)) + 1.0);
                m_sr[i] += (log(Th / (2 * ThetaR)) + 1.0); // From Magin above;

            }
            
            // Old equations, before generalize
            // sr[1] = 1.0 + log((0.5 * Th / 2.1) / N ) + 1.0; // Eq. 3.78 of Boyd. Need to define N or substitute
            // sr[2] = 1.0 + log((0.5 * Th / 2.1) / N ) + 1.0; // Eq. 3.78 of Boyd. Need to define N or substitute

        } else {
            // sR(Tr, s, PlusEq());
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    s[i] += 0.0;
                    m_sr[i] += 0.0;
                    continue; }
                s[i] += (log(Th / (2 * ThetaR)) + 1.0); // From Magin above;
                m_sr[i] += (log(Th / (2 * ThetaR)) + 1.0); // From Magin above;
            }

        }

        // etc...

        // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
        if (sv != NULL) {
            // LOOP(sv[i] = 0.0);
            // sV(Tv, sv, Eq());
            // LOOP_MOLECULES(s[j] += sv[j]);
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    sv[i] += 0.0; // Ground state
                    m_sv[i] += 0.0;
                    s[i] += 0.0;
                 continue; } // Ground state
                sv[i] += 0.0;//m_hv[i] / Th - 1/Th * log(1-exp(-1.0*energy[i-1] * 1.42879/Th)); // Setting to 0 based on discussion with George -- no degeneracy, don't lose any info since sts
                s[i] += 0.0;
            }
            // Old equations, before generalize
            // sv[1] = 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th; // Eq. 3.78 of Boyd. Need to define N or substitute
            // sv[2] =  1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        } else {
            // sV(Tv, s, PlusEq());
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    s[i] += 0.0; // Ground state
                    m_sv[i] += 0.0;
                 continue; } // Ground state
                s[i] += 0.0;//m_hv[i] / Th - 1/Th * log(1-exp(-1.0*energy[i-1] * 1.42879/Th));
                m_sv[i] += 0.0;//m_hv[i] / Th - 1/Th * log(1-exp(-1.0*energy[i-1] * 1.42879/Th));
            }
            // Old equations, before generalize
            // s[1] += 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th;
            // s[2] += 1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        }

        double g0_O2 = 3.0;
        double g1_O2 = 2.0;
        double theta_1_O2 = 11900;
        double g0_O = 5.0;
        double g1_O = 4.0;
        double theta_1_O = 270;

        if (sel != NULL) {
            // LOOP(sv[i] = 0.0);
            // sV(Tv, sv, Eq());
            // LOOP_MOLECULES(s[j] += sv[j]);
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    sel[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Th*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel))); // Ground state
                    m_sel[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel)));
                    s[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel)));
                 continue; } // Ground state
                sel[i] += (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Th*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel))); // Setting to 0 based on discussion with George -- no degeneracy, don't lose any info since sts
                m_sel[i] +=  (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
                s[i] += (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Th*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
            }
            // Old equations, before generalize
            // sv[1] = 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th; // Eq. 3.78 of Boyd. Need to define N or substitute
            // sv[2] =  1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        } else {
            // sV(Tv, s, PlusEq());
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    s[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel))); // Ground state
                    m_sel[i] += (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel)));
                 continue; } // Ground state
                s[i] += (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
                m_sel[i] += (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
            }
            // Old equations, before generalize
            // s[1] += 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th;
            // s[2] += 1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        }


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
        double energy[50];
        energy[0] = 786.0234;
        energy[1] = 2343.573;
        energy[2] = 3881.3038;
        energy[3] = 5398.5964;
        energy[4] = 6894.8782;
        energy[5] = 8369.5118;
        energy[6] = 9822.0457;
        energy[7] = 11251.754;
        energy[8] = 12658.072;
        energy[9] = 14040.4352;
        energy[10] = 15398.3596;
        energy[11] = 16731.2;
        energy[12] = 18038.3111;
        energy[13] = 19319.2091;
        energy[14] = 20573.1679;
        energy[15] = 21799.7037;
        energy[16] = 22998.2518;
        energy[17] = 24168.0058;
        energy[18] = 25308.4816;
        energy[19] = 26419.0341;
        energy[20] = 27498.9373;
        energy[21] = 28547.546;
        energy[22] = 29564.2148;
        energy[23] = 30548.1374;
        energy[24] = 31498.6683;
        energy[25] = 32414.9205;
        energy[26] = 33296.168;
        energy[27] = 34141.6041;
        energy[28] = 34950.3419;
        energy[29] = 35721.3326;
        energy[30] = 36453.8504;
        energy[31] = 37146.6854;
        energy[32] = 37798.8698;
        energy[33] = 38409.2744;
        energy[34] = 38976.5281;
        energy[35] = 39499.5822;
        energy[36] = 39976.9044;
        energy[37] = 40407.0429;
        energy[38] = 40788.6264;
        energy[39] = 41119.9613;
        energy[40] = 41399.6763;
        energy[41] = 41626.4003;
        energy[42] = 41799.6494;
        energy[43] = 41920.3914;
        energy[44] = 41993.5464;
        energy[45] = 42029.6803;
        energy[46] = 42042.9885;

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
        // for (int i = 0; i < m_ns; i++){
        //     g[i] = 0.;
        // }
        // entropy(Th, Te, Tr, Tv, Tel, P, g, NULL, NULL, NULL, NULL);
        // for (int i = 0; i < m_ns; i++){
        //     g[i] = 0.;
        // }
        // // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        // // Will need to upload masses of each species
        if (gt != NULL) {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    gt[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                    continue;
                }
                gt[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // G = H - TS
                g[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
            }
        } else {
            for (int i = 0; i < m_ns; i++){
                if (i == 0){
                    g[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0159994 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                    continue;
                }
                g[i] -= 2.5 * log(Th) - log(P) + log(pow((2*PI*0.0319988 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // G = H - TS
            }
        }

        double ThetaR = 2.08; //char temp rot O2
        if (gr != NULL) {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    gr[i] -= 0.0; // Ground state
                    g[i] -= 0.0;
                    continue; }
                gr[i] -= (log(Tr / (2 * ThetaR)) + 1.0); // G = H - TS
                g[i] -= (log(Tr / (2 * ThetaR)) + 1.0);
            }
        } else {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    g[i] -= 0.0;
                    continue; }
                g[i] -= (log(Tr / (2 * ThetaR)) + 1.0);  // G = H - TS
            }
        }

        if (gv != NULL) {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    gv[i] -= 0.0; // Ground state
                    g[i] -= 0.0;
                 continue; } // Ground state
                gv[i] -= 0.0; // G = H - TS // Tv?
                g[i] -= 0.0;
            }
        } else {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    g[i] -= 0.0;
                 continue; } // Ground state
                g[i] -= 0.0;
            }
        }
        
        double g0_O2 = 3.0;
        double g1_O2 = 2.0;
        double theta_1_O2 = 11900;
        double g0_O = 5.0;
        double g1_O = 4.0;
        double theta_1_O = 270;
        if (gel != NULL) {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    gel[i] -= (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Th*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel))); // Ground state
                    g[i] -= (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Th*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel)));
                 continue; } // Ground state
                gel[i] -= (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Th*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel))); // Setting to 0 based on discussion with George -- no degeneracy, don't lose any info since sts
                g[i] -= (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Th*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
            }
        } else {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    g[i] -= (log(g0_O) + log(1.0+g1_O/g0_O*exp(-theta_1_O/Tel)) + (g1_O/g0_O*theta_1_O/Tel*exp(-theta_1_O/Tel))/(1+(g1_O/g0_O)*exp(-theta_1_O/Tel))); // Ground state
                 continue; } // Ground state
                g[i] -= (log(g0_O2) + log(1.0+g1_O2/g0_O2*exp(-theta_1_O2/Tel)) + (g1_O2/g0_O2*theta_1_O2/Tel*exp(-theta_1_O2/Tel))/(1+(g1_O2/g0_O2)*exp(-theta_1_O2/Tel)));
            }
        }
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
                    "thermodynamics", "type", "RRHO");

            if (rrho_iter == species_iter->end())
                continue;

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
        m_vh.resize(m_ns);
        m_vhf.resize(m_ns);
        m_ht.resize(m_ns);
        m_hr.resize(m_ns);
        m_hv.resize(m_ns);
        m_hel.resize(m_ns);
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
        m_vhf[0] = 0.0;
        m_vhf[1] = 0.0;
        m_vhf[2] = 0.0;
        m_vhf[3] = 0.0;
        m_vhf[4] = 0.0;
        m_vhf[5] = 0.0;
        m_vhf[6] = 0.0;
        m_vhf[7] = 0.0;
        m_vhf[8] = 0.0;
        m_vhf[9] = 0.0;
        m_vhf[10] = 0.0;
        // m_vhf[11] = 0.0;
        // m_vhf[12] = 0.0;
        // m_vhf[13] = 0.0;
        // m_vhf[14] = 0.0;
        // m_vhf[15] = 0.0;
        // m_vhf[16] = 0.0;
        // m_vhf[17] = 0.0;
        // m_vhf[18] = 0.0;
        // m_vhf[19] = 0.0;
        // m_vhf[20] = 0.0;
        // m_vhf[21] = 0.0;
        // m_vhf[22] = 0.0;
        // m_vhf[23] = 0.0;
        // m_vhf[24] = 0.0;
        // m_vhf[25] = 0.0;
        // m_vhf[26] = 0.0;
        // m_vhf[27] = 0.0;
        // m_vhf[28] = 0.0;
        // m_vhf[29] = 0.0;
        // m_vhf[30] = 0.0;
        // m_vhf[31] = 0.0;
        // m_vhf[32] = 0.0;
        // m_vhf[33] = 0.0;
        // m_vhf[34] = 0.0;
        // m_vhf[35] = 0.0;
        // m_vhf[36] = 0.0;
        // m_vhf[37] = 0.0;
        // m_vhf[38] = 0.0;
        // m_vhf[39] = 0.0;
        // m_vhf[40] = 0.0;
        // m_vhf[41] = 0.0;
        // m_vhf[42] = 0.0;
        // m_vhf[43] = 0.0;
        // m_vhf[44] = 0.0;
        // m_vhf[45] = 0.0;
        // m_vhf[46] = 0.0;
        // m_vhf[0] = 9448.87932;
        // m_vhf[1] = 37795.51728;
        // m_vhf[2] = 85039.91388000001;
        // m_vhf[3] = 151182.06912;
        // m_vhf[4] = 236221.983;
        // m_vhf[5] = 340159.65552000003;
        // m_vhf[6] = 462995.08668000007;
        // m_vhf[7] = 604728.27648;
        // m_vhf[8] = 765359.22492;
        // m_vhf[9] = 944887.9319999999;
        // m_vhf[10] = 1143314.39772;
        // m_vhf[11] = 1360638.6220800001;
        // m_vhf[12] = 1596860.6050800001;
        // m_vhf[13] = 1851980.3467200003;
        // m_vhf[14] = 2125997.847;
        // m_vhf[15] = 2418913.10592;
        // m_vhf[16] = 2730726.12348;
        // m_vhf[17] = 3061436.8996800003;
        // m_vhf[18] = 3411045.43452;
        // m_vhf[19] = 3779551.728;
        // m_vhf[20] = 4166955.7801200002;
        // m_vhf[21] = 4573257.59088;
        // m_vhf[22] = 4998457.16028;
        // m_vhf[23] = 5442554.4883200005;
        // m_vhf[24] = 5905549.575;
        // m_vhf[25] = 6387442.4203200005;
        // m_vhf[26] = 6888233.0242800005;
        // m_vhf[27] = 7407921.38688;
        // m_vhf[28] = 7946507.50812;
        // m_vhf[29] = 8503991.388;
        // m_vhf[30] = 9080373.02652;
        // m_vhf[31] = 9675652.42368;
        // m_vhf[32] = 10289829.57948;
        // m_vhf[33] = 10922904.49392;
        // m_vhf[34] = 11574877.167;
        // m_vhf[35] = 12245747.59872;
        // m_vhf[36] = 12935515.78908;
        // m_vhf[37] = 13644181.73808;
        // m_vhf[38] = 14371745.44572;
        // m_vhf[39] = 15118206.912;
        // m_vhf[40] = 15883566.136920001;
        // m_vhf[41] = 16667823.120480001;
        // m_vhf[42] = 17470977.86268;
        // m_vhf[43] = 18293030.36352;
        // m_vhf[44] = 19133980.623;
        // m_vhf[45] = 19993828.64112;
        // m_vhf[46] = 20872574.417880002;



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

    // Store here only the necessary data for calculating species thermodynamics
    const int m_ns = 11; // need to see how to recognize number of states from M++
    const int m_na = 1; // need to see how to recognize number of states from M++
    const int m_nm = 10; // need to see how to recognize number of states from M++
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
        LOOP(op(cp[i], 1.5));
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
        // int ilevel = 0;
        // double sum, fac1, fac2;
        op(cp[0], 0.0);
        LOOP_ATOMS(op(cp[j], 0.0));
        LOOP_MOLECULES(op(cp[j], 0.0));
//!        LOOP_MOLECULES(
//!            sum = 0.0;
//!            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
//!                fac1 = mp_vib_temps[ilevel] / Tv;
//!                fac2 = std::exp(fac1);
//!                fac1 *= fac1*fac2;
//!                fac2 -= 1.0;
//!                sum += fac1/(fac2*fac2);
//!            }
//!            op(cp[j], sum);
//!        )
    }
    /**
     * Computes the translational enthalpy of each species in K.
     */
    // template <typename OP>
    // void hT(double T, double Te, double* const h, const OP& op) {
    //     if (m_has_electron)
    //         op(h[0], 2.5 * Te);
    //     LOOP_HEAVY(op(h[j], 2.5 * T))
    // }

    /**
     * Computes the rotational enthalpy of each species in K.
     */
    template <typename OP>
    void hR(double T, double* const h, const OP& op) {
//!        LOOP_MOLECULES(op(h[j], mp_rot_data[i].linearity * T))
    }

    /**
     * Computes the vibrational enthalpy of each species in K.
     */
    template <typename OP>
    void hV(double T, double* const h, const OP& op) {
        if (T < 10.0) {
            LOOP_MOLECULES(op(h[j], 0.0));
        } else {
//!            int ilevel = 0;
//!            double sum;
//!            LOOP_MOLECULES(
//!                sum = 0.0;
//!                for (int k = 0; k < mp_nvib[i]; ++k, ilevel++)
//!                    sum += mp_vib_temps[ilevel] /
//!                        (std::exp(mp_vib_temps[ilevel] / T) - 1.0);
//!                op(h[j], sum);
//!            )
        }
    }
    /**
     * Computes the unitless translational entropy of each species.
     */
    // template <typename OP>
    // void sT(double Th, double Te, double P, double* const s, const OP& op) {
    //     double fac = 2.5 * (1.0 + std::log(Th)) - std::log(P);
    //     if (m_has_electron)
    //         op(s[0], 2.5 * std::log(Te / Th) + fac + mp_lnqtmw[0]);
    //     for (int i = (m_has_electron ? 1 : 0); i < m_ns; ++i)
    //         op(s[i], fac + mp_lnqtmw[i]);
    // }

    /**
     * Computes the unitless rotational entropy of each species.
     */
    template <typename OP>
    void sR(double T, double* const s, const OP& op) {
        const double onelnT = 1.0 + std::log(T);
//!        LOOP_MOLECULES(
//!            op(s[j], mp_rot_data[i].linearity * (onelnT -
//!                mp_rot_data[i].ln_omega_t));
//!        )
    }

    /**
     * Computes the unitless vibrational entropy of each species.
     */
    template <typename OP>
    void sV(double T, double* const s, const OP& op) {
        int ilevel = 0;
        double fac, sum1, sum2;
//!        LOOP_MOLECULES(
//!            sum1 = sum2 = 0.0;
//!            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
//!                fac  =  std::exp(mp_vib_temps[ilevel] / T);
//!                sum1 += mp_vib_temps[ilevel] / (fac - 1.0);
//!                sum2 += std::log(1.0 - 1.0 / fac);
//!            }
//!            op(s[j], (sum1 / T - sum2));
//!        )
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
    std::vector<double> m_sv {}; //should this be in private??
    std::vector<double> m_st {}; //should this be in private??
    std::vector<double> m_sr {}; //should this be in private??
    std::vector<double> m_sel {}; //should this be in private??



    //     spnm  spwt(g/mol) ih  ie
    //   O2   32.0000     2   0
    //  h0sp (kJ/mol)
    //   0.0
    //  factr(homogeneous)
    //   0.5d0
    //  dissociation energy (cm-1)
    //   41280.0
    //  spectral data
    //   1 ! number of electronic excited states
    //   n  state     Te        re      g  dzero     we        wexe     weye
    //      weze       be       alphae     de         betae    spn-orb lambda spin
    //   1  X3SIGg-   0.00      1.2075  3  41280.00  1580.193  11.9808  4.747E-02
    //     -1.273E-03  1.43768  1.593E-02  4.839E-06  0.000E+00 -8.400E-03  0  3
    // 36
    //    0   222   7.87380953594E+02
    //    1   218   2.34376026609E+03
    //    2   213   3.87656829159E+03
    //    3   209   5.38602874609E+03
    //    4   205   6.87233479359E+03
    //    5   200   8.33564904609E+03
    //    6   196   9.77610356359E+03
    //    7   192   1.11937998541E+04
    //    8   188   1.25888088736E+04
    //    9   183   1.39611710261E+04
    //   10   179   1.53108961636E+04
    //   11   175   1.66379635861E+04
    //   12   170   1.79423220416E+04
    //   13   166   1.92238897261E+04
    //   14   162   2.04825542836E+04
    //   15   157   2.17181728061E+04
    //   16   153   2.29305718336E+04
    //   17   148   2.41195473541E+04
    //   18   143   2.52848648036E+04
    //   19   138   2.64262590661E+04
    //   20   134   2.75434344736E+04
    //   21   129   2.86360648061E+04
    //   22   124   2.97037932916E+04
    //   23   118   3.07462326061E+04
    //   24   113   3.17629648736E+04
    //   25   108   3.27535416661E+04
    //   26   102   3.37174840036E+04
    //   27   96    3.46542823541E+04
    //   28   90    3.55633966336E+04
    //   29   84    3.64442562061E+04
    //   30   77    3.72962598836E+04
    //   31   70    3.81187759261E+04
    //   32   63    3.89111420416E+04
    //   33   55    3.96726653861E+04
    //   34   45    4.04026225636E+04
    //   35   34    4.11002596261E+04
    //   36   19    4.17647920736E+04


}; // class STS10DB
#undef LOOP
#undef LOOP_HEAVY
#undef LOOP_MOLECULES


// Register the STS10DB model with the other thermodynamic databases
Utilities::Config::ObjectProvider<STS10DB, ThermoDB> sts10DB("STS10");

    } // namespace Thermodynamics
} // namespace Mutation


