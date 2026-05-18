//
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-.
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| |
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | |
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  |
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)|
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_)
//      (__)              (_)      (__)     (__)     (__)
//      Official webSite: https://code-mphi.github.io/ECOGEN/
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names
//  are listed in the copyright file included with this source
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published
//  by the Free Software Foundation, either version 3 of the License,
//  or (at your option) any later version.
//
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef GRADPHASESHALLOWWATER_H
#define GRADPHASESHALLOWWATER_H

#include "../GradPhase.h"
#include "PhaseShallowWater.h"

class Phase;

//! \class     GradPhaseShallowWater
//! \brief     Phase variable gradients for ShallowWater model. Stored for 2nd-order computation on unstructured mesh (O2 NS)
class GradPhaseShallowWater : public GradPhase
{
  public:
    GradPhaseShallowWater();
    ~GradPhaseShallowWater() override;

    void initializeGradientVectors() override;

    void computeDistanceGradientScalarProduct(Coord const& distance, Phase* phase) const override;
    void limitGradients(const Phase& gradientLimiter) override;

    int numberOfTransmittedGradients() const override;

  protected:
    //! \brief     Enumeration for the phase flow variables, specific to ShallowWater
    enum VarLocal
    {
      height,
      velocityU,
      velocityV
    };
};

#endif
