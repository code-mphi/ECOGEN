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

#ifndef MIXSHALLOWWATER_H
#define MIXSHALLOWWATER_H

#include "../Mixture.h"

//! \class     MixShallowWater
//! \brief     Mixture variables for ShallowWater equations (single phase)
class MixShallowWater : public Mixture
{
  public:
    MixShallowWater() {};
    ~MixShallowWater() override {};

    void allocateAndCopyMixture(Mixture** mixture) override { *mixture = new MixShallowWater(*this); };
    void copyMixture(Mixture& /*mixture*/) override {};

    void computeTotalEnergy(std::vector<QuantitiesAddPhys*>& /*vecGPA*/) override {};                                 // mandatory
    void localProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override {};   // mandatory
    void reverseProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override {}; // mandatory

    //Specific methods for data printing
    //----------------------------------

    //Specific methods for parallel computing
    //---------------------------------------
    int numberOfTransmittedVariables() const override { return 0; };
    void fillBuffer(double* /*buffer*/, int& /*counter*/) const override {};
    void fillBuffer(std::vector<double>& /*dataToSend*/) const override {};
    void getBuffer(double* /*buffer*/, int& /*counter*/) override {};
    void getBuffer(std::vector<double>& /*dataToReceive*/, int& /*counter*/) override {};

    //Accessors
    //---------

    //Operators
    //---------

  protected:
  private:
};

#endif // MIXSHALLOWWATER_H
