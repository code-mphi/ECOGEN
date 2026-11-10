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

#ifndef NONSYMMETRICCONTACTAREAS_H
#define NONSYMMETRICCONTACTAREAS_H

#include "ContactAreas.h"

//! \class     NonSymmetricContactAreas
//! \brief     NonSymmetricContactAreas class for the non-equilibrium system of equations, algorithm from the PhD thesis of Ashwinn Chinnayya (2004)
class NonSymmetricContactAreas : public ContactAreas
{
  public:
    //! \brief  Constructor of the non-symmetric contact areas algorithm
    NonSymmetricContactAreas(const int& numbPhases);
    //! \brief  Destructor of the non-symmetric contact areas algorithm
    ~NonSymmetricContactAreas() override;
    //! \brief  Compute the contact areas (non-symmetric)
    //! \param  alphaLeft    left volume fraction
    //! \param  alphaRight   right volume fraction
    void computeContactAreas(const std::vector<double>& alphaLeft, const std::vector<double>& alphaRight) override;

    std::vector<double> m_contactSurfacesLeft;  //!< left contact areas, from Chinnayya thesis (2004) algorithm
    std::vector<double> m_contactSurfacesRight; //!< right contact areas, from Chinnayya thesis (2004) algorithm
};

#endif // NONSYMMETRICCONTACTAREAS_H
