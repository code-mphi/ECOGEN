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

#include "NonSymmetricContactAreas.h"

//***********************************************************************

NonSymmetricContactAreas::NonSymmetricContactAreas(const int& numbPhases) : ContactAreas(numbPhases)
{
  m_contactSurfacesRight.resize(numbPhases);
  m_contactSurfacesLeft.resize(numbPhases);
}

//***********************************************************************

NonSymmetricContactAreas::~NonSymmetricContactAreas() {}

//***********************************************************************

void NonSymmetricContactAreas::computeContactAreas(const std::vector<double>& alphaLeft, const std::vector<double>& alphaRight)
{
  // Contact areas for N-phases with Chinnayya's method

  // Initialisation
  // --------------
  for (int i = 0; i < numberPhases; i++) {
    m_contactSurfacesLeft[i]  = alphaLeft[i];
    m_contactSurfacesRight[i] = alphaRight[i];
  }

  // Contact intra-phase (same phase)
  // --------------------------------
  for (int i = 0; i < numberPhases; i++) {
    double sii                 = std::min(m_contactSurfacesLeft[i], m_contactSurfacesRight[i]);
    m_contacts[i][i]           = sii;
    m_contactSurfacesLeft[i]  -= sii;
    m_contactSurfacesRight[i] -= sii;
  }

  // Contact between phase 1 and 2 (or 2 and 1)
  // ------------------------------------------
  for (int i = 0; i < numberPhases; i++) {
    for (int j = 0; j < numberPhases; j++) {
      if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
        double sij                 = std::min(m_contactSurfacesLeft[i], m_contactSurfacesRight[j]);
        m_contacts[i][j]           = sij;
        m_contactSurfacesLeft[i]  -= sij;
        m_contactSurfacesRight[j] -= sij;
      }
    }
  }

  // Contact between others phases
  //------------------------------
  for (int i = 0; i < numberPhases; i++) {
    for (int j = 0; j < numberPhases; j++) {
      if (i != j && !(i == 1 && j == 2) && !(i == 2 && j == 1)) {
        double sij                 = std::min(m_contactSurfacesLeft[i], m_contactSurfacesRight[j]);
        m_contacts[i][j]           = sij;
        m_contactSurfacesLeft[i]  -= sij;
        m_contactSurfacesRight[j] -= sij;
      }
    }
  }

  // Recalibrate contact areas
  //--------------------------
  this->recalibrateContactAreas();
}

//***********************************************************************
