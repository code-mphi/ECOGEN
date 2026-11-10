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

#include "SymmetricContactAreas.h"

//***********************************************************************

SymmetricContactAreas::SymmetricContactAreas(const int& numbPhases) : ContactAreas(numbPhases)
{
  m_deltaTotal = 0.;
  m_delta.resize(numbPhases);
  m_deltaPos.resize(numbPhases);
  m_deltaNeg.resize(numbPhases);
}

//***********************************************************************

SymmetricContactAreas::~SymmetricContactAreas() {}

//***********************************************************************

void SymmetricContactAreas::computeContactAreas(const std::vector<double>& alphaLeft, const std::vector<double>& alphaRight)
{
  // Contact areas for N-phases with Paula's method

  // Contact areas between same phases
  //----------------------------------
  for (int i = 0; i < numberPhases; i++) {
    m_contacts[i][i] = std::min(alphaLeft[i], alphaRight[i]);
  }

  // Contact areas between different phases
  //---------------------------------------
  m_deltaTotal = 0.;
  for (int i = 0; i < numberPhases; i++) {
    m_delta[i]     = alphaLeft[i] - alphaRight[i];
    m_deltaPos[i]  = std::max(0.0, m_delta[i]);
    m_deltaNeg[i]  = std::max(0.0, -m_delta[i]);
    m_deltaTotal  += std::abs(m_delta[i]);
  }
  m_deltaTotal *= 0.5;

  double epsilon(1.e-15);
  if (m_deltaTotal > epsilon) {
    for (int i = 0; i < numberPhases; i++) {
      for (int j = 0; j < numberPhases; j++) {
        if (i != j) {
          m_contacts[i][j] = (m_deltaPos[i] * m_deltaNeg[j]) / m_deltaTotal;
        }
      }
    }
  }
  else {
    for (int i = 0; i < numberPhases; i++)
      for (int j = 0; j < numberPhases; j++)
        if (i != j) m_contacts[i][j] = 0.;
  }

  // Recalibrate contact areas
  //--------------------------
  this->recalibrateContactAreas();
}

//***********************************************************************
