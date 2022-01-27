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

#ifndef SOURCE_H
#define SOURCE_H

#include "../libTierces/tinyxml2.h"
#include "../Errors.h"
#include "../Tools.h"
#include "../Order1/Cell.h"
#include "../Models/Flux.h"

//! \class     Source
//! \brief     Abstract class for source terms
class Source
{
  public:
    //! \brief    Source constructor depending on integration order and physical entity to apply source
    //! \param    physicalEntity  the entity to which the source term is applied (default whole domain)
    Source(int physicalEntity = 0);
    virtual ~Source();

    //! \brief     Source terms preparation for integration
    //! \param     cell           cell for source term integration
    virtual void prepSourceTerms(Cell* /*cell*/, const int& /*i*/ = 0) { Errors::errorMessage("prepSourceTerms not available for required source"); };

    //! \brief     Source terms integration on conservative quantities
    //! \param     cell           cell for source term integration
    //! \param     dt             integration time step
    virtual void integrateSourceTerms(Cell* /*cell*/, const double& /*dt*/) { Errors::errorMessage("integrateSourceTerms not available for required source"); };

    //! \brief     Allows to modifiy the source term along time
    //! \param     time            physical time of the computation
    virtual void sourceEvolution(const double& /*time*/) {};

    //! \brief     Compute the absolute velocity in the fixed coordinate system
    //! \param     relVelocity     velocity in the moving coordinate system
    //! \param     position        position vector in the fixed coordinate system
    virtual Coord computeAbsVelocity(const Coord& /*relVelocity*/, const Coord& /*position*/) { Errors::errorMessage("computeAbsVelocity not available for required source"); return 0.; };

    int getPhysicalEntity() { return m_physicalEntity; }

  protected:
    int m_physicalEntity;
};

#endif // SOURCE_H
