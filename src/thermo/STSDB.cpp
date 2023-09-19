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
#include <fstream>


using namespace std;
using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

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

    STSDB(int arg) : ThermoDB(298.15, 101325.0){}
        //   m_has_electron(false),
        //   m_use_tables(false),
        //   m_last_bfacs_T(0.0) {} //

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
            double* const cpv, double* const cpel) {
        
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

     // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
     if (cpt != NULL) {
         for (int i = 0; i < m_ns; i++){
             cpt[i] += 2.5; // Cv = 3/2 R; Cp = Cv + R
             cp[i] += cpt[i];
         }

     } else {
         for (int i = 0; i < m_ns; i++){
             cp[i] += 2.5;
         }
     }

     // Rotation. Assuming fulling active rotational mode
     if (cpr != NULL) {
         for (int i = 0; i < m_ns; i++){
             if (i == 0) {
//                 cp[i] = 0.0; // Ground state
//                 cp[i] += cpr[i];
                 continue; } // Ground state
             cpr[i] += 2.0; //iCv = R; Cp = Cv + R
             cp[i] += cpr[i];
         }

     } else {
         for (int i = 0; i < m_ns; i++){
             if (i == 0) {
//                 cp[i] = 0.0; // Ground state
                 continue; } // Ground state
             cp[i] += 2.0;
         }
     }

     // etc...

     // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
     if (cpv != NULL) {
         for (int i = 0; i < m_ns; i++){
             cpv[i] += 0.0; // Setting as zero for now. Need to think
             cpv[i] += cpv[i];
         }

     } else {
         for (int i = 0; i < m_ns; i++){
             cpv[i] += 0.0;
         }
     }

     // Electronic. For now setting as zero
     // cpel[0] = 0.0;
     // cpel[1] = 0.0;
     // cpel[2] = 0.0;

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
        double Th, double Te, double Tr, double Tv, double Tel,
        double* const h, double* const ht, double* const hr,
        double* const hv, double* const hel, double* const hf)
    {
        // Given Ts calculate h, ht, hr, hv, hel, hf
        // Old equations, before generalize
    //    h[0] = 2.5 + m_vhf[0];
    //    h[1] = m_vh[1]*Th + m_vhf[1];
    //    h[2] = m_vh[2]*Th + m_vhf[2];

        const int SIZE = 37;
        double state[SIZE];
        double blank[SIZE];
        double energy[SIZE];
        string inFileName = "../../data/thermo/oxygen_energy.txt";
        ifstream inFile;
        inFile.open(inFileName.c_str());
        if (inFile.is_open())
        {
            for (int i = 0; i < SIZE; i++)
            {
                inFile >> state[i];
                inFile >> blank[i];
                inFile >> energy[i];
            }

            inFile.close(); // CLose input file
        }
        
        
        
        
        
        // Setting to zero
        for (int i = 0; i < m_ns; i++){
            h[i] = 0.;
        }

        // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        if (ht != NULL) {
            for (int i = 0; i < m_ns; i++){
                ht[i] += 2.5; // Is this non-dimensional too? Taking in work flow too. Otherwise it would be 1.5.
                h[i] += ht[i];
            }

        } else {
            for (int i = 0; i < m_ns; i++){
                h[i] += 2.5;
            }
        }

        // Rotation. Assuming fulling active rotational mode
        if (hr != NULL) {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
//                    h[i] = 0.0; // Ground state
//                    h[i] += hr[i];
                    continue; }
                hr[i] += 1.0;
                h[i] += hr[i];
            }

        } else {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
//                    h[i] = 0.0; // Ground state
                    continue; }
                h[i] += 1.0;
            }
        }

        // etc...

        // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
        if (hv != NULL) {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    hv[i] = 0.0;
//                    h[i] = 0.0; // Ground state
                    continue; }
                hv[i] = energy[i] * 1.42879 / Th * exp(-1*energy[i] * 1.42879 / Th); // See KMH notes
                h[i] += hv[i];
            }

        } else {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    h[i] = 0.0;
//                    h[i] = 0.0; // Ground state
                    continue; }
                h[i] = energy[i] * 1.42879 / Th * exp(-1*energy[i] * 1.42879 / Th); // See KMH notes
            }
        }

        // Electronic. For now setting as zero
        // hel[0] = 0.0;
        // hel[1] = 0.0;
        // hel[2] = 0.0;

        for (int i = 0; i < m_ns; i++){
            h[i] += m_vhf[i];
        }

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
        for (int i = 0; i < m_ns; i++){
            s[i] = 0.;
        }


        // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        // Will need to upload masses of each species
        if (st != NULL) {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*15.999 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // EQ 3.90 of Boyd. // Ground state
                    s[i] += st[i];
                    continue; }
                st[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*31.998 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
                s[i] += st[i];
            }

        } else {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    s[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*15.999 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5; // EQ 3.90 of Boyd. // Ground state
                    s[i] += st[i];
                    continue; }
                s[i] += 2.5 * log(Th) - log(P) + log(pow((2*PI*31.998 / NA / pow(HP,2.0)),1.5) * pow(KB,2.5)) + 2.5;
            }
        }

        // Rotation. Assuming fulling active rotational mode
        if (sr != NULL) {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    sr[i] = 0.0; // Ground state
                    s[i] += sr[i];
                    continue; }
                sr[i] = hr[i]/Th + log(0.5 * Th / 2.1); // From slide 20 of Magin, need to check units
                s[i] += sr[i];
            }
            
            // Old equations, before generalize
            // sr[1] = 1.0 + log((0.5 * Th / 2.1) / N ) + 1.0; // Eq. 3.78 of Boyd. Need to define N or substitute
            // sr[2] = 1.0 + log((0.5 * Th / 2.1) / N ) + 1.0; // Eq. 3.78 of Boyd. Need to define N or substitute

        } else {
            for (int i = 0; i < m_ns; i++){
                if (i == 0) {
                    s[i] += 0.0;
                    continue; }
                s[i] += hr[i]/Th + log(0.5 * Th / 2.1); // From Magin above;
            }

        }

        // etc...

        // Vibration. Assuming the characteristic vib temperature is the vib energy level of that state.
        if (sv != NULL) {
            for (int i = 0; i < m_ns; i++){
                sv[i] = 0.0; // Setting to 0 based on discussion with George -- no degeneracy, don't lose any info since sts
                s[i] += sv[i];
            }
            // Old equations, before generalize
            // sv[1] = 1.0 + log(exp(-7.87380953594E+02 * 1.42879 / Th) / N ) + 7.87380953594E+02 * 1.42879 / Th; // Eq. 3.78 of Boyd. Need to define N or substitute
            // sv[2] =  1.0 + log(exp(-2.34376026609E+03 * 1.42879 / Th) / N ) + 2.34376026609E+03 * 1.42879 / Th;
        } else {
            for (int i = 0; i < m_ns; i++){
                s[i] += 0.0;
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
        // Given Ts calculate g, gt, gr, gv, gel
        // Note: Check if NULL

                   // Following similar approach as enthalpy
            // Setting to zero
        for (int i = 0; i < m_ns; i++){
            g[i] = 0.;
        }

        // Eventually, replace this with a loop over all species as they should have equal translational enthalpy
        // Will need to upload masses of each species
        if (gt != NULL) {
            for (int i = 0; i < m_ns; i++){
                gt[i] += ht[i] - Th * st[i]; // G = H - TS
            }

            for (int i = 0; i < m_ns; i++){
                g[i] += gt[i];
            } // Is there a way to combine this with previous loop?
  
        } else {
            for (int i = 0; i < m_ns; i++){
                gt[i] += ht[i] - Th * st[i]; // G = H - TS
            }
        }

        if (gr != NULL) {
            for (int i = 0; i < m_ns; i++){
                gr[i] += hr[i] - Th * sr[i]; // G = H - TS
            }

            for (int i = 0; i < m_ns; i++){
                g[i] += gr[i];
            } // Is there a way to combine this with previous loop?
  
        } else {
            for (int i = 0; i < m_ns; i++){
                gr[i] += hr[i] - Th * sr[i]; // G = H - TS
            }
        }

        if (gv != NULL) {
            for (int i = 0; i < m_ns; i++){
                gv[i] += hv[i] - Th * st[i]; // G = H - TS // Tv?
            }

            for (int i = 0; i < m_ns; i++){
                g[i] += gv[i];
            } // Is there a way to combine this with previous loop?
  
        } else {
            for (int i = 0; i < m_ns; i++){
                gv[i] += hv[i] - Th * sv[i]; // G = H - TS // Tv?
            }
        }
        // Electronic. For now setting as zero
        // sel[0] = 0.0;
        // sel[1] = 0.0;
        // sel[2] = 0.0;

        //h[0] += m_vhf[0];
        //h[1] += m_vhf[1];
        // h[2] += m_vhf[2];
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
        // Add ht, hr, hv...

        m_vh[0] = 0.; // Atomic oxygen
        m_vh[1] = 7.87380953594E+02;
        m_vh[2] = 1.4E+03;
        // m_vh[3] = 7.87380953594E+02;

        // I think these are enthlapy of formation
        m_vhf[0] = 100.59; // Atomic oxygen
        m_vhf[1] = 7.87380953594E+02;
        m_vhf[2] = 1.0+02;
        // m_vh[3] = 7.87380953594E+02;

        // @todo: 1/18/2023
        // Load the necessary thermodynamic data
    }

private:
    // Store here only the necessary data for calculating species thermodynamics
    const int m_ns = 3; // need to see how to recognize number of states from M++
    std::vector<double> m_vh {};
    std::vector<double> m_vhf {};
    std::vector<double> hv {}; //should this be in private?? //m_ for private
    std::vector<double> ht {}; //should this be in private??
    std::vector<double> hr {}; //should this be in private??
    std::vector<double> sv {}; //should this be in private??
    std::vector<double> st {}; //should this be in private??
    std::vector<double> sr {}; //should this be in private??



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


}; // class STSDB

// Register the STSDB model with the other thermodynamic databases
Utilities::Config::ObjectProvider<STSDB, ThermoDB> stsDB("STS");

    } // namespace Thermodynamics
} // namespace Mutation


