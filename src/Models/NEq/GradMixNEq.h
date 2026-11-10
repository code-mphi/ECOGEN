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

#ifndef GRADMIXNEQ_H
#define GRADMIXNEQ_H

#include "../GradMixture.h"
#include "MixNEq.h"

class Mixture;

//! \class     GradMixNEq
//! \brief     Mixture variable gradients for NEq model. Stored for 2nd-order computation on unstructured mesh (O2 NS)

class GradMixNEq : public GradMixture
{
  public:
    GradMixNEq() {};
    ~GradMixNEq() override {};

    void initializeGradientVectors() override {};

    void computeDistanceGradientScalarProduct(Coord const& /*distance*/, Mixture* /*mixture*/) const override {};
    void limitGradients(const Mixture& /*gradientLimiter*/) override {};

    // O2 parallel
    int numberOfTransmittedGradients() const override { return 0; };
    void getBufferGradients(double* /*buffer*/, int& /*counter*/) override {};
    void fillBufferGradients(double* /*buffer*/, int& /*counter*/) override {};
};

#endif //GRADMIXNEQ_H
