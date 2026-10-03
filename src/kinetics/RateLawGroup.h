/**
 * @file RateLawGroup.h
 *
 * Defines the abstract class RateLawGroup which is used by the RateManager
 * class to combine rate laws which are evaluated at the same temperature in
 * order to compute the rate coefficient for these reactions efficiently.
 *
 * @see class RateLawGroup
 * @see class RateLawGroup1T
 * @see class RateLawGroupCollection
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

#ifndef KINETICS_RATE_LAW_GROUP_H
#define KINETICS_RATE_LAW_GROUP_H

#include <map>
#include <typeinfo>
#include <vector>

#include "RateLaws.h"
#include "Reaction.h"
//#include "StateModel.h"
#include "StoichiometryManager.h"

class StateModel;

namespace Mutation {
    namespace Kinetics {

/**
 * Abstract base class which defines the interface for all RateLawGroup objects
 * which evaluate like rate raws to increase efficiency.
 */
class RateLawGroup
{
public:

    /**
     * Constructor.
     */
    RateLawGroup() : m_last_t(-1.0) {}

    /**
     * Destructor.
     */
    virtual ~RateLawGroup() { };

    /**
     * Adds a new rate to evaluate with this group.
     */
    virtual void addRateCoefficient(
        const size_t rxn, const RateLaw* const p_rate) = 0;
    
    /**
     * Adds a reaction which uses the backward temperature represented by this
     * group.
     */
    void addReaction(const size_t rxn, const Reaction& reaction) {
        m_reacs.addReaction(rxn, reaction.reactants());
        m_prods.addReaction(rxn, reaction.products());
    }
    
    /**
     * Returns the temperature used in the last evaluation of the rate
     * coefficients.
     */
    double getT() const { return m_t; }
    
    /**
     * Evaluates all of the rates in the group and stores in the given vector.
     */
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk) = 0;
        
    /**
     * Evaluates all of the temperature derivative of rates in the group and stores in the given vector.
     */
    virtual void invkdkdT(
        const Thermodynamics::StateModel* const p_state, double* const p_dkdT) = 0;
  
    /**
     * Computes \Delta G / RT for this rate law group and subtracts these values
     * for each of the reactions in this group.
     */
    void subtractLnKeq(size_t ns, double* const p_g, double* const p_r) const
    {        
        // Compute G_i/RT - ln(Patm/RT)
        const double val = std::log(ONEATM / (RU * m_t));
        for (size_t i = 0; i < ns; ++i)
            p_g[i] -= val;
            
        // Now subtract \Delta[G_i/RT - ln(Patm/RT)]_j
        m_reacs.decrReactions(p_g, p_r);
        m_prods.incrReactions(p_g, p_r);
    }
    
    /**
     * Computes 1 / Keq \frac{\partial Keq}{\partial T} for this rate law group 
     * and subtracts these values for each of the reactions in this group.
     */
    void derivativeKeq(double* const p_dKeqdT, double* const p_dkdT) const
    {        
        m_reacs.decrReactions(p_dKeqdT, p_dkdT);
        m_prods.incrReactions(p_dKeqdT, p_dkdT);
    }

protected:

    /// This is the temperature computed to evaluate the rate law (should be set
    /// in the lnk() function)
    double m_t;
    double m_last_t;
    
    /// Stores the reactants for reactions that will use this rate law for the
    /// reverse direction
    StoichiometryManager m_reacs;
    
    /// Stores the products for reactions that will use this rate law for the
    /// reverse direction
    StoichiometryManager m_prods;
};


/**
 * Groups reaction rate laws based on a single temperature that are the same
 * kind and evaluated at the same temperature together so that they may be 
 * evaluated efficiently.
 */
template <typename RateLawType, typename TSelectorType>
class RateLawGroup1T : public RateLawGroup
{
public:

    typedef TSelectorType Selector;

    /**
     * Constructor.  Selectors which carry parameters (ie: TaTvSelector) are
     * passed in, otherwise the default constructed selector is used.
     */
    explicit RateLawGroup1T(const TSelectorType& selector = TSelectorType())
        : m_selector(selector)
    { }

    /**
     * Adds a new rate to evaluate with this group.
     */
    virtual void addRateCoefficient(
        const size_t rxn, const RateLaw* const p_rate)
    {
        m_rates.push_back(
            std::make_pair(rxn, dynamic_cast<const RateLawType&>(*p_rate))
        );
    }

    /**
     * Evaluates all of the rates in the group and stores in the given vector.
     */
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk)
    {
        // Determine the reaction temperature for this group
        m_t = m_selector.getT(p_state);

        // Update only if the temperature has changed
        //if (std::abs(m_t - m_last_t) > 1.0e-10) {
            const double lnT  = std::log(m_t);
            const double invT = 1.0 / m_t;

            for (int i = 0; i < m_rates.size(); ++i) {
                const std::pair<size_t, RateLawType>& rate = m_rates[i];
                p_lnk[rate.first] = rate.second.getLnRate(lnT, invT);
            }
        //}

        // Save this temperature
        m_last_t = m_t;
    }

    /**
     * Evaluates all of the temperature derivative rates in the group and stores in the given vector.
     */
    virtual void invkdkdT(
        const Thermodynamics::StateModel* const p_state, double* const p_dkdT)
    {
        // Determine the reaction temperature for this group
        m_t = m_selector.getT(p_state);

        // Update only if the temperature has changed
        //if (std::abs(m_t - m_last_t) > 1.0e-10) {
            const double invT = 1.0 / m_t;

            for (size_t i = 0; i < m_rates.size(); ++i) {
                const std::pair<size_t, RateLawType>& rate = m_rates[i];
                p_dkdT[rate.first] = rate.second.derivativebykf(invT);
            }
        //}

        // Save this temperature
        m_last_t = m_t;
    }

private:

    /// Selects the temperature at which the rates in this group are evaluated
    TSelectorType m_selector;

    /// vector of rates to evaluate
    std::vector< std::pair<size_t, RateLawType> > m_rates;
};


/**
 * Small helper class which provides comparison between two std::type_info
 * pointers.
 */
struct CompareTypeInfo {
    bool operator ()(const std::type_info* a, const std::type_info* b) const {
        return a->before(*b);
    }
};


/**
 * Manages a collection of RateLawGroup objects such that only one object of any
 * RateLawGroup type is ever created in each collection.
 */
class RateLawGroupCollection
{
public:

    typedef std::map<const std::type_info*, RateLawGroup*, CompareTypeInfo>
        GroupMap;

    /// Groups evaluated at T^a * Tv^b, keyed by (a, b)
    typedef std::map<std::pair<double, double>, RateLawGroup*> TaTvGroupMap;

    /**
     * Destructor.
     */
    ~RateLawGroupCollection()
    {
        for (size_t i = 0; i < m_groups.size(); ++i)
            delete m_groups[i];
    }
    
    /**
     * Returns the number of different rate law groups in this collection.
     */
    size_t nGroups() const { return m_groups.size(); }
    
    /**
     * Returns all of the RateLawGroup objects in this collection.
     */
    const std::vector<RateLawGroup*>& groups() const { return m_groups; }

    /**
     * Adds a new rate law to be managed by this collection of rate law groups.
     */
    template <typename GroupType>
    void addRateCoefficient(const size_t rxn, const RateLaw* const p_rate)
    {
        getGroup<GroupType>()->addRateCoefficient(rxn, p_rate);
    }

    /**
     * Adds a new rate law to be evaluated at T^a * Tv^b.  One group is created
     * for each unique (a, b) pair.  GroupType::Selector must be constructible
     * from (a, b).
     */
    template <typename GroupType>
    void addRateCoefficient(
        const size_t rxn, const RateLaw* const p_rate,
        const double a, const double b)
    {
        RateLawGroup*& p_group = m_tatv_map[std::make_pair(a, b)];
        if (p_group == NULL) {
            p_group = new GroupType(typename GroupType::Selector(a, b));
            m_groups.push_back(p_group);
        }
        p_group->addRateCoefficient(rxn, p_rate);
    }
    
    /**
     * Adds a reaction to the manager which allows for the calculation of the 
     * \Delta G / RT term for reverse rate coefficients.
     */
    template <typename GroupType>
    void addReaction(const size_t rxn, const Reaction& reaction)
    {
        getGroup<GroupType>()->addReaction(rxn, reaction);
    }

    /**
     * Computes the rate coefficients in this collection and stores them in
     * the vector at the index corresponding to their respective reaction.
     */
    void logOfRateCoefficients(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk)
    {
        // Compute the forward rate constants
        for (size_t i = 0; i < m_groups.size(); ++i)
            m_groups[i]->lnk(p_state, p_lnk);
    }
    
    /**
     * Computes the derivative of rate coefficients 
     * \f$ \frac{1}{k_{f,j}} \frac{dk_{f,j}}{dT_{reac}} \f$
     * in this collection and stores them in the vector at the index corresponding 
     * to their respective reaction.
     */
    void derivativeOfRateCoefficients(
        const Thermodynamics::StateModel* const p_state, double* const p_dkdT)
    {
        // Compute the forward rate constants
        for (size_t i = 0; i < m_groups.size(); ++i)
            m_groups[i]->invkdkdT(p_state, p_dkdT);
    }

    /**
     * Subtracts ln(keq) from the provided rate coefficients.
     */
    void subtractLnKeq(
        const Thermodynamics::Thermodynamics& thermo, double* const p_g,
        double* const p_lnk)
    {
        const size_t ns = thermo.nSpecies();
        for (size_t i = 0; i < m_groups.size(); ++i) {
            const RateLawGroup* p_group = m_groups[i];
            thermo.speciesSTGOverRT(p_group->getT(), p_g);
            p_group->subtractLnKeq(ns, p_g, p_lnk);
        }
    }

    /**
     * Computes the temperature derivative of keq as follows:
     * \f$ \frac{1}{keq} \frac{\partial keq}{\partial T} \f$
     */
    void derivativeKeq(
        const Thermodynamics::Thermodynamics& thermo, double* const p_dKeqdT, double* const p_dkdT)
    {
        const size_t ns = thermo.nSpecies();
        for (size_t g = 0; g < m_groups.size(); ++g) {
            const RateLawGroup* p_group = m_groups[g];
            thermo.speciesSTdGOverRT(p_group->getT(), p_dKeqdT);
	    for(size_t i = 0; i < ns; ++i)
                p_dKeqdT[i] += 1./p_group->getT();
            p_group->derivativeKeq(p_dKeqdT, p_dkdT);
        }
    }
  
private:

    /**
     * Returns the group of the given type, creating it if necessary.
     */
    template <typename GroupType>
    RateLawGroup* getGroup()
    {
        RateLawGroup*& p_group = m_group_map[&typeid(GroupType)];
        if (p_group == NULL) {
            p_group = new GroupType();
            m_groups.push_back(p_group);
        }
        return p_group;
    }
    
    /// Lookup of RateLawGroup objects with compile-time temperature selectors
    GroupMap m_group_map;

    /// Lookup of RateLawGroup objects evaluated at T^a * Tv^b
    TaTvGroupMap m_tatv_map;

    /// All RateLawGroup objects in this collection (owns the pointers)
    std::vector<RateLawGroup*> m_groups;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // KINETICS_RATE_LAW_GROUP_H

