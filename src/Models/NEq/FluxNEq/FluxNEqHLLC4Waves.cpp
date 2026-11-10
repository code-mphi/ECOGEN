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

#include "FluxNEqHLLC4Waves.h"

//***********************************************************************

FluxNEqHLLC4Waves::FluxNEqHLLC4Waves(const int& numbPhases) : FluxNEq(numbPhases) {}

//***********************************************************************

FluxNEqHLLC4Waves::~FluxNEqHLLC4Waves() {}

//***********************************************************************

void FluxNEqHLLC4Waves::setToZero()
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    = 0.;
    m_mass[k]     = 0.;
    m_momentum[k] = 0.;
    m_energ[k]    = 0.;
  }
  m_uI = 0.;
  m_pI = 0.;
}

//***********************************************************************

void FluxNEqHLLC4Waves::addNonCons(double coefA, const Cell* cell, const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  // Saurel & Abgrall closure law for non-conservative terms of Furfaro & Saurel solver
  // ----------------------------------------------------------------------------------
  Phase* phase;
  for (int k = 0; k < numberPhases; k++) {
    phase       = cell->getPhase(k);
    m_alpha[k] += -coefA * phase->getAlpha() * static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_uI; // Topological equation choice in ECOGEN
    m_momentum[k].localProjection(normal, tangent, binormal);
    m_momentum[k].setX(m_momentum[k].getX() + coefA * phase->getAlpha() * static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_pI);
    m_momentum[k].reverseProjection(normal, tangent, binormal);
    m_energ[k] += coefA * phase->getAlpha() * static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_pI * static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_uI;
  }
}

//***********************************************************************

void FluxNEqHLLC4Waves::subtractNonCons(double coefA, const Cell* cell, const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  // Saurel & Abgrall closure law for non-conservative terms of Furfaro & Saurel solver
  // ----------------------------------------------------------------------------------
  Phase* phase;
  for (int k = 0; k < numberPhases; k++) {
    phase       = cell->getPhase(k);
    m_alpha[k] -= -coefA * phase->getAlpha() * static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_uI;
    m_momentum[k].localProjection(normal, tangent, binormal);
    m_momentum[k].setX(m_momentum[k].getX() - coefA * phase->getAlpha() * static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_pI);
    m_momentum[k].reverseProjection(normal, tangent, binormal);
    m_energ[k] -= coefA * phase->getAlpha() * static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_pI * static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_uI;
  }
}

//***********************************************************************
