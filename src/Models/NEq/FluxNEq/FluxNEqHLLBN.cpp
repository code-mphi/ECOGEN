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

#include "FluxNEqHLLBN.h"

//***********************************************************************

FluxNEqHLLBN::FluxNEqHLLBN(const int& numbPhases) : FluxNEq(numbPhases)
{
  m_alphaL = new double[numbPhases];
  m_alphaR = new double[numbPhases];
  m_sL     = new double[numbPhases];
  m_sR     = new double[numbPhases];
}

//***********************************************************************

FluxNEqHLLBN::~FluxNEqHLLBN()
{
  delete[] m_alphaL;
  delete[] m_alphaR;
  delete[] m_sL;
  delete[] m_sR;
}

//***********************************************************************

void FluxNEqHLLBN::setToZero()
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    = 0.;
    m_mass[k]     = 0.;
    m_momentum[k] = 0.;
    m_energ[k]    = 0.;
    m_alphaL[k]   = 0.;
    m_alphaR[k]   = 0.;
    m_sL[k]       = 0.;
    m_sR[k]       = 0.;
  }
  m_uI = 0.;
  m_pI = 0.;
}

//***********************************************************************

void FluxNEqHLLBN::addNonCons(double coefA, const Cell* cell, const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  // HLL/Rusanov non-conservative terms of Lochon (2016) with HLL with B&N closure law
  // ----------------------------------------------------------------------------------
  double uI = cell->getPhase(1)->getVelocity().getX();
  double vI = cell->getPhase(1)->getVelocity().getY();
  double wI = cell->getPhase(1)->getVelocity().getZ();
  velocity.setXYZ(uI, vI, wI);
  velocity.localProjection(normal, tangent, binormal);

  double pI = cell->getPhase(0)->getPressure();
  Phase* phase;
  for (int k = 0; k < numberPhases; k++) {
    phase = cell->getPhase(k);
    phase->getVelocity().localProjection(normal, tangent, binormal);
    m_alpha[k] += coefA *
                  (velocity.getX() * (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] * static_cast<FluxNEqHLLBN*>(fluxBuff)->m_alphaL[k] -
                                      static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k] * phase->getAlpha()) +
                   static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] * static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k] *
                     (phase->getAlpha() - static_cast<FluxNEqHLLBN*>(fluxBuff)->m_alphaL[k])) /
                  (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] - static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k]);
    m_momentum[k].localProjection(normal, tangent, binormal);
    m_momentum[k].setX(m_momentum[k].getX() - coefA * pI *
                                                ((static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] * static_cast<FluxNEqHLLBN*>(fluxBuff)->m_alphaL[k] -
                                                  static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k] * phase->getAlpha()) /
                                                 (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] - static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k])));
    m_momentum[k].reverseProjection(normal, tangent, binormal);
    m_energ[k] += -coefA * pI * velocity.getX() * // velocity.getX() = uI projected on the normal of the face
                  ((static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] * static_cast<FluxNEqHLLBN*>(fluxBuff)->m_alphaL[k] -
                    static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k] * phase->getAlpha()) /
                   (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] - static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k]));
    phase->getVelocity().reverseProjection(normal, tangent, binormal);
  }
  velocity.reverseProjection(normal, tangent, binormal);
}

//***********************************************************************

void FluxNEqHLLBN::subtractNonCons(double coefA, const Cell* cell, const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  // HLL/Rusanov non-conservative terms of Lochon (2016) with HLL with B&N closure law
  // ----------------------------------------------------------------------------------
  double uI = cell->getPhase(1)->getVelocity().getX();
  double vI = cell->getPhase(1)->getVelocity().getY();
  double wI = cell->getPhase(1)->getVelocity().getZ();

  velocity.setXYZ(uI, vI, wI);
  velocity.localProjection(normal, tangent, binormal);

  double pI = cell->getPhase(0)->getPressure();

  Phase* phase;
  for (int k = 0; k < numberPhases; k++) {
    phase = cell->getPhase(k);
    phase->getVelocity().localProjection(normal, tangent, binormal);
    m_alpha[k] -= coefA *
                  (velocity.getX() * (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] * phase->getAlpha() -
                                      static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k] * static_cast<FluxNEqHLLBN*>(fluxBuff)->m_alphaR[k]) +
                   static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] * static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k] *
                     (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_alphaR[k] - phase->getAlpha())) /
                  (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] - static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k]);
    m_momentum[k].localProjection(normal, tangent, binormal);
    m_momentum[k].setX(m_momentum[k].getX() + coefA * pI *
                                                ((static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] * phase->getAlpha() -
                                                  static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k] * static_cast<FluxNEqHLLBN*>(fluxBuff)->m_alphaR[k]) /
                                                 (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] - static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k])));
    m_momentum[k].reverseProjection(normal, tangent, binormal);
    m_energ[k] -= -coefA * pI * velocity.getX() * // velocity.getX() = uI projected on the normal of the face
                  ((static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] * phase->getAlpha() -
                    static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k] * static_cast<FluxNEqHLLBN*>(fluxBuff)->m_alphaR[k]) /
                   (static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sR[k] - static_cast<FluxNEqHLLBN*>(fluxBuff)->m_sL[k]));
    phase->getVelocity().reverseProjection(normal, tangent, binormal);
  }
  velocity.reverseProjection(normal, tangent, binormal);
}

//***********************************************************************
