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

#ifndef CONTACTAREAS_H
#define CONTACTAREAS_H

#include "../../../Order1/Cell.h"

//! \class     Contact areas
//! \brief     Contact areas abstract class for the non-equilibrium system of equations
class ContactAreas
{
  public:
    //! \brief  Constructor of the contact areas
    ContactAreas(const int& numbPhases);
    //! \brief  Virtual destructor of the contact areas
    virtual ~ContactAreas();
    //! \brief  Compute the contact areas, either symmetric or non-symmetric
    //! \param  alphaLeft    left volume fraction
    //! \param  alphaRight   right volume fraction
    virtual void computeContactAreas(const std::vector<double>& /*alphaLeft*/, const std::vector<double>& /*alphaRight*/)
    {
      Errors::errorMessage("computeContactAreas not available for required ContactArea");
    };
    //! \brief  Recalibration of contact areas
    void recalibrateContactAreas();

    std::vector<std::vector<double>> m_contacts; //!< contact areas [] (non-dimensional, since based on sum of volume fraction = 1)
};

#endif // CONTACTAREAS_H
