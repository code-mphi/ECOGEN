//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

//! \file      Source.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      January 10 2018

#include <string>
#include "../libTierces/tinyxml2.h"
#include "../Errors.h"
#include "../Tools.h"
#include "../Cell.h"

//! \class     Source
//! \brief     Abstract class for source terms
class Source
{
  public:
    Source();
    virtual ~Source();

    //! \brief     Source terms integration on conservative quantities
    //! \param     cell           cell for source term integration
    //! \param     numberPhases   number of phases
    //! \param     dt             explicit integration time step
    virtual void integrateSourceTerms(Cell *cell, const int &numberPhases, const double &dt){ Errors::errorMessage("integrateSourceTerms not available for required source"); };
    virtual void sourceEvolution(const double &time) {};

    virtual Coord computeAbsVelocity(const Coord relVelocity, const Coord position) { Errors::errorMessage("computeAbsVelocity not available for required source"); return 0.; };
};

#endif // SOURCE_H