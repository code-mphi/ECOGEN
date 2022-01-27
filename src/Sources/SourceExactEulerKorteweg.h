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

#ifndef SOURCEEXACTEULERKORTEWEG_H
#define SOURCEEXACTEULERKORTEWEG_H

#include "SourceExact.h"
#include "../Models/EulerKorteweg/FluxEulerKorteweg.h"

//! \class     SourceExactEulerKorteweg
//! \brief     Class for Euler--Korteweg source terms
class SourceExactEulerKorteweg : public SourceExact
{
  public:
    //! \brief    SourceExactEulerKorteweg constructor depending on physical entity to apply source
    //! \param    physicalEntity  the entity to which the source term is applied (default whole domain)
    SourceExactEulerKorteweg(int physicalEntity = 0);
    virtual ~SourceExactEulerKorteweg();

    //! \brief     Integration via an exact solution
    //! \param     cell           cell for source term integration
    //! \param     dt             explicit integration time step
    virtual void integrationExactSolution(Cell* cell, const double& dt);
};

#endif // SOURCEEXACTEULERKORTEWEG_H
