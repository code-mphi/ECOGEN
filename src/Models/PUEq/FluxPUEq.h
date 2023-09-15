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

#ifndef FLUXPUEQ_H
#define FLUXPUEQ_H

#include "../UEq/FluxUEq.h"

class FluxPUEq;

#include "ModPUEq.h"

//! \class     FluxPUEq
//! \brief     Flux class for the pressure-velocity-equilibrium (mechanical equilibrium) system of equations (Kapila)
class FluxPUEq : public FluxUEq
{
  public:
    FluxPUEq(const int& numbPhases);
    virtual ~FluxPUEq();

    virtual void addNonCons(double coefA, const Cell* cell, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/);
    virtual void subtractNonCons(double coefA, const Cell* cell, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/);
    virtual void schemeCorrection(Cell& /*cell*/) const {};
    virtual void correctionEnergy(Cell* cell, Prim type = vecPhases) const;

  private:
    friend class ModPUEq;
};

#endif // FLUXPUEQ_H