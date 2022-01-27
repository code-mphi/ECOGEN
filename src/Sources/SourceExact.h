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

#ifndef SOURCEEXACT_H
#define SOURCEEXACT_H

#include "Source.h"

//! \class     SourceExact
//! \brief     Abstract class for source terms solved by an exact solution
class SourceExact : public Source
{
  public:
    //! \brief    SourceExact constructor depending on physical entity to apply source
    //! \param    physicalEntity  the entity to which the source term is applied (default whole domain)
    SourceExact(int physicalEntity = 0);
    virtual ~SourceExact();

    //! \brief     Source terms integration on conservative quantities
    //! \param     cell           cell for source term integration
    //! \param     dt             integration time step
    virtual void integrateSourceTerms(Cell* cell, const double& dt);

    //! \brief     Integration via an exact solution
    //! \param     cell           cell for source term integration
    //! \param     dt             explicit integration time step
    virtual void integrationExactSolution(Cell* /*cell*/, const double& /*dt*/) { Errors::errorMessage("integrationExactSolution not available for required source"); };
};

#endif // SOURCEEXACT_H
