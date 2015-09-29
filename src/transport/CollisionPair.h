/**
 * @file CollisionPair.h
 *
 * Provides the CollisionPair type.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2015 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef TRANSPORT_COLLISION_PAIR_H
#define TRANSPORT_COLLISION_PAIR_H

#include "CollisionIntegral.h"
#include "SharedPtr.h"
#include "XMLite.h"
#include <string>

// Forward declarations
namespace Mutation { namespace Thermodynamics { class Species; }}

namespace Mutation {
    namespace Transport {

/**
 * Enumerates the different type of collision pairs.
 */
enum CollisionType {
    NEUTRAL_NEUTRAL,
    ELECTRON_NEUTRAL,
    ION_NEUTRAL,
    ATTRACTIVE,
    REPULSIVE
};

/**
 * Encapsulates data corresponding to a particular collision pair.
 */
class CollisionPairNew
{
public:
    /**
     * Loads the collision pair information provided the species in the pair and
     * the root node of the XML collision database.
     */
    CollisionPairNew(
        const Mutation::Thermodynamics::Species& s1,
        const Mutation::Thermodynamics::Species& s2,
        const Utilities::IO::XmlElement& xml);

    // Getter functions
    const std::string& species1() const { return m_sp1; }
    const std::string& species2() const { return m_sp2; }

    SharedPtr<CollisionIntegral> Q11() const { return mp_Q11; }
    SharedPtr<CollisionIntegral> Q22() const { return mp_Q22; }
    SharedPtr<CollisionIntegral> Bst() const { return mp_Bst; }
    SharedPtr<CollisionIntegral> Cst() const { return mp_Cst; }

private:

    /**
     * Initializes the species names and the type of collision integral this is
     * from the two species objects.
     */
    void initSpeciesData(
        const Mutation::Thermodynamics::Species& s1,
        const Mutation::Thermodynamics::Species& s2);

    /**
     * Returns the iterator pointing to the XmlElement which holds the collision
     * integral of type kind for this collision pair.  If this doesn't exist,
     * then database.end() is returned.
     */
    Mutation::Utilities::IO::XmlElement::const_iterator
    findXmlElementWithIntegralType(
        const std::string& kind,
        const Mutation::Utilities::IO::XmlElement& database) const;

    /**
     * Loads a particular collision integral from the database.
     */
    SharedPtr<CollisionIntegral> loadIntegral(
        const std::string& type,
        const Mutation::Utilities::IO::XmlElement& database) const;

private:

    CollisionType m_type;
    std::string   m_sp1;
    std::string   m_sp2;

    // collision integrals
    SharedPtr<CollisionIntegral> mp_Q11;
    SharedPtr<CollisionIntegral> mp_Q22;
    SharedPtr<CollisionIntegral> mp_Bst;
    SharedPtr<CollisionIntegral> mp_Cst;
};

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_COLLISION_PAIR_H
