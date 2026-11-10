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

#ifndef FLUXNEQDEM_H
#define FLUXNEQDEM_H

#include "../NEq/FluxNEq/FluxNEq.h"
#include "ContactAreas/HeaderContactAreas.h"

class FluxNEqDEM;

//! \class     FluxNEqDEM
//! \brief     Flux class for the non-equilibrium system of equations using the DEM
class FluxNEqDEM : public FluxNEq
{
  public:
    FluxNEqDEM(const int& numbPhases);
    ~FluxNEqDEM() override;

    void setToZero() override;
    void addNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal) override;
    void subtractNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal) override;

    double** m_uI;    //!< interfacial velocity [m/s] //KS//Shouldn't m_uI and m_pI be extern variables?
    double** m_pI;    //!< interfacial pressure [Pa]
    double** m_pI_02; //!< interfacial pressure for 2nd order
    double** m_uI_02; //!< interfacial velocity for 2nd order

    friend class ModNEqDEM; // !< ModNEqDEM need to access to contact areas and define interfacial variables of interface fluxes
};

extern ContactAreas* contactAreas;    //!< Contact area pointer to either symmetric or non-symmetric contact areas algorithms.
extern ContactAreas* contactAreas_02; //!< Contact area 02 pointer to either symmetric or non-symmetric contact areas algorithms.

#endif // FLUXNEQDEM_H
