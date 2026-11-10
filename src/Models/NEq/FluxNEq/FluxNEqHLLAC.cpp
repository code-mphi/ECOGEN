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

#include "FluxNEqHLLAC.h"

//***********************************************************************

FluxNEqHLLAC::FluxNEqHLLAC(const int& numbPhases) : FluxNEq(numbPhases)
{
  m_alphaL = new double[numbPhases];
  m_alphaR = new double[numbPhases];
}

//***********************************************************************

FluxNEqHLLAC::~FluxNEqHLLAC()
{
  delete[] m_alphaL;
  delete[] m_alphaR;
}

//***********************************************************************

void FluxNEqHLLAC::setToZero()
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    = 0.;
    m_mass[k]     = 0.;
    m_momentum[k] = 0.;
    m_energ[k]    = 0.;
    m_alphaL[k]   = 0.;
    m_alphaR[k]   = 0.;
  }
  m_pStarL = 0.;
  m_pStarR = 0.;
  m_uI     = 0.;
  m_pI     = 0.;
}

//***********************************************************************

void FluxNEqHLLAC::addNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  // Non-conservative terms of HLLAC
  // -------------------------------
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += -coefA * 0.5 * (1. + sign(static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI)) * std::pow(-1., k) *
                  static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI *
                  (static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaR[0] - static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaL[0]);
    m_momentum[k].localProjection(normal, tangent, binormal);
    m_momentum[k].setX(m_momentum[k].getX() - coefA * 0.5 * (1. + sign(static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI)) * std::pow(-1., k) *
                                                (static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaR[1] * static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarR -
                                                 static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaL[1] * static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarL));
    m_momentum[k].reverseProjection(normal, tangent, binormal);
    m_energ[k] += -coefA * 0.5 * (1. + sign(static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI)) * std::pow(-1., k) *
                  static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI *
                  (static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaR[1] * static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarR -
                   static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaL[1] * static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarL);
  }
}

//***********************************************************************

void FluxNEqHLLAC::subtractNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  // Non-conservative terms of HLLAC
  // -------------------------------
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] -= -coefA * 0.5 * (1. - sign(static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI)) * std::pow(-1., k) *
                  static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI *
                  (static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaR[0] - static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaL[0]);
    m_momentum[k].localProjection(normal, tangent, binormal);
    m_momentum[k].setX(m_momentum[k].getX() + coefA * 0.5 * (1. - sign(static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI)) * std::pow(-1., k) *
                                                (static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaR[1] * static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarR -
                                                 static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaL[1] * static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarL));
    m_momentum[k].reverseProjection(normal, tangent, binormal);
    m_energ[k] -= -coefA * 0.5 * (1. - sign(static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI)) * std::pow(-1., k) *
                  static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI *
                  (static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaR[1] * static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarR -
                   static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaL[1] * static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarL);
  }
}

//***********************************************************************
