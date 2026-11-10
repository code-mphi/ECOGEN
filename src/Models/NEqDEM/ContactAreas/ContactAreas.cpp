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

#include "ContactAreas.h"

//***********************************************************************

ContactAreas::ContactAreas(const int& numbPhases)
{
  m_contacts.resize(numbPhases);
  for (int k = 0; k < numbPhases; k++) {
    m_contacts[k].resize(numbPhases);
    for (int i = 0; i < numbPhases; i++) {
      m_contacts[k][i] = 0.;
    }
  }
}

//***********************************************************************

ContactAreas::~ContactAreas() {}

//***********************************************************************

void ContactAreas::recalibrateContactAreas()
{
  //Recalibrate contact areas
  //-------------------------
  //To make sure the sum of all contact areas equals 1
  double one(0.);
  double epsilon(1.e-12);
  for (int i = 0; i < numberPhases; i++) {
    for (int j = 0; j < numberPhases; j++) {
      if (m_contacts[i][j] < epsilon) m_contacts[i][j] = 0.;
      one += m_contacts[i][j];
    }
  }
  for (int i = 0; i < numberPhases; i++) {
    for (int j = 0; j < numberPhases; j++) {
      m_contacts[i][j] /= one;
    }
  }
}
