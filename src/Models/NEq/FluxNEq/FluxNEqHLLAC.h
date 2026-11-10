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

#ifndef FLUXNEQHLLAC_H
#define FLUXNEQHLLAC_H

#include "FluxNEq.h"

//! \class     FluxNEqHLLAC
//! \brief     Flux class for the non-equilibrium system of equations using HLLAC solver (Lochon et al. 2016)
class FluxNEqHLLAC : public FluxNEq
{
  public:
    FluxNEqHLLAC(const int& numbPhases);
    ~FluxNEqHLLAC() override;

    void setToZero() override;
    void addNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal) override;
    void subtractNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal) override;

    double* m_alphaL; //!< array of left volume fraction, for non-conservative terms
    double* m_alphaR; //!< array of right volume fraction, for non-conservative terms
    double m_pStarL;  //!< left star pressure of phase 1, for non-conservative terms
    double m_pStarR;  //!< right star pressure of phase 1, for non-conservative terms
};

#endif // FLUXNEQHLLAC_H
