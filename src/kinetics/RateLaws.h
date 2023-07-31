/**
 * @file RateLaws.h
 *
 * @brief Declaration of various RateLaw classes.
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

#ifndef RATELAW_H
#define RATELAW_H

#include <vector>
#include <cmath>
#include <cstdlib>

#include "Utilities.h"

namespace Mutation {
    namespace Kinetics {

/**
 * Abstract base class for all rate laws which allows owners such as class
 * Reaction to store any rate law polymorphically.
 */
class RateLaw
{
public:

    virtual ~RateLaw() { };
    virtual RateLaw* clone() const = 0;
};

/**
 * Arrhenius rate law \f$ k_f(T) = A T^\eta exp(-E_a / (R_u T)) \f$.
 */
class Arrhenius : public RateLaw
{
public:

    static void setUnits(const Mutation::Utilities::IO::XmlElement& node);
    
    Arrhenius(const Mutation::Utilities::IO::XmlElement& node, const int order);
    
    Arrhenius(const Arrhenius& to_copy)
        : m_lnA(to_copy.m_lnA), m_n(to_copy.m_n), m_temp(to_copy.m_temp)
    { }
    
    virtual ~Arrhenius() { };
    
    Arrhenius* clone() const {
        return new Arrhenius(*this);
    }
    
    inline double getLnRate(const double lnT, const double invT) const {
        return (m_lnA + m_n * lnT - m_temp * invT);
    }
    
    inline double derivative(const double k, const double lnT, const double invT) const {
        return (k*invT*(m_n + m_temp*invT));
    }

    double A() const {
        return std::exp(m_lnA);
    }
    
    double n() const {
        return m_n;
    }
    
    double T() const {
        return m_temp;
    }
    
private:

    static std::vector<Mutation::Utilities::Units> sm_aunits;
    static std::vector<Mutation::Utilities::Units> sm_eunits;

    double m_lnA;
    double m_n;
    double m_temp;
};
    
/**
 * Arrhenius rate law \f$ k_f(T) = A T^\eta exp(-E_a / (R_u T)) \f$.
 */
class MMT : public RateLaw
{
public:

    static void setUnits(const Mutation::Utilities::IO::XmlElement& node);
    
    MMT(const Mutation::Utilities::IO::XmlElement& node, const int order);
    
    MMT(const MMT& to_copy)
        : m_lnA(to_copy.m_lnA), m_n(to_copy.m_n), m_temp(to_copy.m_temp), m_temp(to_copy.m_a), m_temp(to_copy.m_U_s)
    { }
    
    virtual ~MMT() { };
    
    MMT* clone() const {
        return new MMT(*this);
    }
    
    inline double getLnRate(const double lnTtr, const double invTtr, const double invTv) const {
//         Possible log approximation Taylor Series?
//        val1 = m_lnA + m_n * lnT - m_temp * invT;
//        invT == 1/Ttr? if so, we can simplify division in U,TF,lnQTR
        U = (1/invTtr * m_U_s) / \
                (1/invTtr + m_a * m_U_s);
        TF = -1 * (1/invTtr * 1/invTv * U) \
                / (1/invTtr * 1/invTv - 1/invTtr * U
                   + 1/invTtr * U);
        lnQTr = std::log(1 - std::exp(-m_temp*invTtr)) - \
                std::log(1 - std::exp(-m_theta_v*invTtr));
        lnQTF = std::log(1 - std::exp(-m_temp/TF)) - \
                std::log(1 - std::exp(-m_theta_v/TF));
        lnQTv = std::log(1 - std::exp(-m_temp*invTv)) - \
                std::log(1 - std::exp(-m_theta_v*invTv));
        lnQU =  std::log(1 - std::exp(m_temp/U)) - \
                std::log(1 - std::exp(m_theta_v/U));
//        lnZ = lnQTr + lnQTF - lnQTv - lnQU
        return (m_lnA + m_n * lnTtr - m_temp * invT + lnQTr + lnQTF - lnQTv - lnQU);
    }
    
    // Can I change function arguments?
    inline double derivative(const double k, const double lnTtr, const double invTtr, const double invTv)) const {
        // k must be the rate value --> derivative with respect to T
        // Jacobian? [d/dTtr, d/dTv] or just d/dTtr
        //        invT == 1/Ttr? if so, we can simplify division in all variables
        val1 = m_temp * (std::exp(m_temp*invTtr) - 2) / (std::exp(m_temp*invTtr) - 1)
        val2 = m_n / invTtr
        val3 = m_theta_v / (std::exp(m_theta_v*invTtr))
        return (k*pow(invT,2)*(val1 + val2 + val3));
    }

    double A() const {
        return std::exp(m_lnA);
    }
    
    double n() const {
        return m_n;
    }
    
    double T() const {
        return m_temp;
    }
    
    double a() const {
        return m_a;
    }
    
    double U_s() const {
        return m_U_s;
    }
    
    double thetaV() const {
        return m_theta_v;
    }
    
    
private:

    static std::vector<Mutation::Utilities::Units> sm_aunits;
    static std::vector<Mutation::Utilities::Units> sm_eunits;

    double m_lnA;
    double m_n;
    double m_temp; // TD?
    double m_theta_v; // Maybe referenced from other part of M++?
//    double m_temp_Tv; // Maybe referenced from other part of M++?
//    double m_temp_Ttr; // Maybe referenced from other part of M++?
    double m_a;
    double m_U_s;
};


    } // namespace Kinetics
} // namespace Mutation



#endif // RATELAW_HPP
