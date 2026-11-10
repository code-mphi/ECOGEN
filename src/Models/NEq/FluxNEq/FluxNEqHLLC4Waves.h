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

#ifndef FLUXNEQHLLC4WAVES_H
#define FLUXNEQHLLC4WAVES_H

#include "FluxNEq.h"

//! \class     FluxNEqHLLC4Waves
//! \brief     Flux class for the non-equilibrium system of equations using HLLC 4-Waves solver (Furfaro & Saurel, 2015)
class FluxNEqHLLC4Waves : public FluxNEq
{
  public:
    FluxNEqHLLC4Waves(const int& numbPhases);
    ~FluxNEqHLLC4Waves() override;

    void setToZero() override;
    void addNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal) override;
    void subtractNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal) override;
};

#endif // FLUXNEQHLLC4WAVES_H
