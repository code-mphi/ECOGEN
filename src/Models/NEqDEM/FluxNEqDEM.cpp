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

#include "FluxNEqDEM.h"

ContactAreas* contactAreas;
ContactAreas* contactAreas_02;

//***********************************************************************

FluxNEqDEM::FluxNEqDEM(const int& numbPhases) : FluxNEq(numbPhases)
{
  m_pI    = new double*[numbPhases];
  m_uI    = new double*[numbPhases];
  m_pI_02 = new double*[numbPhases];
  m_uI_02 = new double*[numbPhases];
  for (int k = 0; k < numbPhases; k++) {
    m_pI[k]    = new double[numbPhases];
    m_uI[k]    = new double[numbPhases];
    m_pI_02[k] = new double[numbPhases];
    m_uI_02[k] = new double[numbPhases];
  }
}

//***********************************************************************

FluxNEqDEM::~FluxNEqDEM()
{
  for (int k = 0; k < numberPhases; k++) {
    delete[] m_pI[k];
    delete[] m_uI[k];
    delete[] m_pI_02[k];
    delete[] m_uI_02[k];
  }
  delete[] m_pI;
  delete[] m_uI;
  delete[] m_pI_02;
  delete[] m_uI_02;
}

//***********************************************************************

void FluxNEqDEM::setToZero()
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    = 0.;
    m_mass[k]     = 0.;
    m_momentum[k] = 0.;
    m_energ[k]    = 0.;
    for (int i = 0; i < numberPhases; i++) {
      m_uI[k][i]    = 0.;
      m_pI[k][i]    = 0.;
      m_uI_02[k][i] = 0.;
      m_pI_02[k][i] = 0.;
    }
  }
}

//***********************************************************************

void FluxNEqDEM::addNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  double epsilon(1.e-12);
  double coefB(0.);
  for (int i = 0; i < numberPhases; i++) {
    for (int j = 0; j < numberPhases; j++) {
      if (contactAreas->m_contacts[i][j] < epsilon) continue;
      if (i != j) {
        coefB = coefA * contactAreas->m_contacts[i][j];
        if (static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j] >= 0.) {
          m_alpha[i] -= -coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j];
          m_momentum[i].localProjection(normal, tangent, binormal);
          m_momentum[i].setX(m_momentum[i].getX() - coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j]);
          m_momentum[i].reverseProjection(normal, tangent, binormal);
          m_energ[i] -= coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j] * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j];

          m_alpha[j] += -coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j];
          m_momentum[j].localProjection(normal, tangent, binormal);
          m_momentum[j].setX(m_momentum[j].getX() + coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j]);
          m_momentum[j].reverseProjection(normal, tangent, binormal);
          m_energ[j] += coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j] * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j];
        }
      }
    }
  }
}

//***********************************************************************

void FluxNEqDEM::subtractNonCons(double coefA, const Cell* /*cell*/, const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  // First order non-conservative flux (at the interface)
  double epsilon(1.e-12);
  double coefB(0.);
  for (int i = 0; i < numberPhases; i++) {
    for (int j = 0; j < numberPhases; j++) {
      if (contactAreas->m_contacts[i][j] < epsilon) continue;
      if (i != j) {
        coefB = coefA * contactAreas->m_contacts[i][j];
        if (static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j] < 0.) {
          m_alpha[i] -= -coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j];
          m_momentum[i].localProjection(normal, tangent, binormal);
          m_momentum[i].setX(m_momentum[i].getX() - coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j]);
          m_momentum[i].reverseProjection(normal, tangent, binormal);
          m_energ[i] -= coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j] * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j];

          m_alpha[j] += -coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j];
          m_momentum[j].localProjection(normal, tangent, binormal);
          m_momentum[j].setX(m_momentum[j].getX() + coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j]);
          m_momentum[j].reverseProjection(normal, tangent, binormal);
          m_energ[j] += coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j] * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j];
        }
      }
    }
  }

  // Second-order non-conservative flux (intern to the cell)
  if (simulationOrder == SimulationOrder::second) {
    for (int i = 0; i < numberPhases; i++) {
      for (int j = 0; j < numberPhases; j++) {
        if (i != j && contactAreas_02->m_contacts[i][j] > epsilon) {
          coefB       = coefA * contactAreas_02->m_contacts[i][j];
          m_alpha[i] -= -coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI_02[i][j];
          m_momentum[i].localProjection(normal, tangent, binormal);
          m_momentum[i].setX(m_momentum[i].getX() - coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI_02[i][j]);
          m_momentum[i].reverseProjection(normal, tangent, binormal);
          m_energ[i] -= coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI_02[i][j] * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI_02[i][j];

          m_alpha[j] += -coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI_02[i][j];
          m_momentum[j].localProjection(normal, tangent, binormal);
          m_momentum[j].setX(m_momentum[j].getX() + coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI_02[i][j]);
          m_momentum[j].reverseProjection(normal, tangent, binormal);
          m_energ[j] += coefB * static_cast<FluxNEqDEM*>(fluxBuff)->m_pI_02[i][j] * static_cast<FluxNEqDEM*>(fluxBuff)->m_uI_02[i][j];
        }
      }
    }
  }
}

//***********************************************************************
