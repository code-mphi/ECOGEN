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

#ifndef SOURCENUM_H
#define SOURCENUM_H

#include "Source.h"

//! \class     SourceNum
//! \brief     Abstract class for source terms solved by a numerical scheme
class SourceNum : public Source
{
  public:
    //! \brief    SourceNum constructor depending on integration order and physical entity to apply source
    //! \param    order           integration order (could EULER, RK2 or RK4 scheme)
    //! \param    physicalEntity  the entity to which the source term is applied (default whole domain)
    SourceNum(int order, int physicalEntity = 0);
    virtual ~SourceNum();

    //! \brief     Source terms preparation for integration
    //! \param     cell           cell for source term integration
    virtual void prepSourceTerms(Cell* /*cell*/, const int& /*i*/ = 0) { Errors::errorMessage("prepSourceTerms not available for required source"); };

    //! \brief     Source terms integration on conservative quantities
    //! \param     cell           cell for source term integration
    //! \param     dt             integration time step
    virtual void integrateSourceTerms(Cell* cell, const double& dt);

    //! \brief     Euler explicite integration (order 1)
    //! \param     cell           cell for source term integration
    //! \param     dt             explicit integration time step
    void integrationEuler(Cell* cell, const double& dt);

    //! \brief     Runge-Kutta integration (order 2)
    //! \param     cell           cell for source term integration
    //! \param     dt             explicit integration time step
    void integrationRK2(Cell* cell, const double& dt);

    //! \brief     Runge-Kutta integration (order 4)
    //! \param     cell           cell for source term integration
    //! \param     dt             explicit integration time step
    void integrationRK4(Cell* cell, const double& dt);

    //! \brief     Allows to modifiy the source term along time
    //! \param     time            physical time of the computation
    virtual void sourceEvolution(const double& /*time*/) {};

    //! \brief     Compute the absolute velocity in the fixed coordinate system
    //! \param     relVelocity     velocity in the moving coordinate system
    //! \param     position        position vector in the fixed coordinate system
    virtual Coord computeAbsVelocity(const Coord& /*relVelocity*/, const Coord& /*position*/) { Errors::errorMessage("computeAbsVelocity not available for required source"); return 0.; };

  protected:
    int m_order;
};

#endif // SOURCENUM_H
