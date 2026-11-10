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

#ifndef MODNEQHLLAC_H
#define MODNEQHLLAC_H

#include "ModNEq.h"

class ModNEqHLLAC;

#include "../FluxNEq/FluxNEqHLLAC.h"

//! \class     ModNEqHLLAC
//! \brief     Model class for the non-equilibrium system of equations using the HLLAC (Lochon et al. 2016)
class ModNEqHLLAC : public ModNEq
{
  public:
    //! \brief     NEq model constructor with HLLAC solver
    //! \param     numbTransports    number of additional transport equations
    //! \param     numbPhases        number of phases
    ModNEqHLLAC(const int& numbTransports, const int& numbPhases);
    ~ModNEqHLLAC() override;

    void allocateCons(Flux** cons) override;

    //! \brief     HLLAC Riemann solver of Lochon et al. (2016)
    void solver(Cell& cellLeft,
                Cell& cellRight,
                const double& dxLeft,
                const double& dxRight,
                double& dtMax,
                std::vector<double>& /*boundData*/ = DEFAULT_VEC_INTERFACE_DATA) const override; // Riemann between two computed cells

  private:
    static const std::string NAME;

    friend class FluxNEqHLLAC;
    friend class Input;
};

#endif // MODNEQHLLAC_H
