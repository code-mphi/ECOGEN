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

#ifndef SYMMETRICCONTACTAREAS_H
#define SYMMETRICCONTACTAREAS_H

#include "ContactAreas.h"

//! \class     SymmetricContactAreas
//! \brief     SymmetricContactAreas class for the non-equilibrium system of equations, algorithm from Paula et al. (2023)
class SymmetricContactAreas : public ContactAreas
{
  public:
    //! \brief  Constructor of the symmetric contact areas algorithm
    SymmetricContactAreas(const int& numbPhases);
    //! \brief  Destructor of the symmetric contact areas algorithm
    ~SymmetricContactAreas() override;
    //! \brief  Compute the contact areas (symmetric) at the interface
    //! \param  alphaLeft    left volume fraction
    //! \param  alphaRight   right volume fraction
    void computeContactAreas(const std::vector<double>& alphaLeft, const std::vector<double>& alphaRight) override;

    std::vector<double> m_delta;    //!< contact areas between different phases, from Paula et al. (2023) algorithm
    std::vector<double> m_deltaPos; //!< positive contact areas between different phases, from Paula et al. (2023) algorithm
    std::vector<double> m_deltaNeg; //!< negative contact areas between different phases, from Paula et al. (2023) algorithm

    double m_deltaTotal; //!< total contact areas between different phases, from Paula et al. (2023) algorithm
};

#endif // SYMMETRICCONTACTAREAS_H
