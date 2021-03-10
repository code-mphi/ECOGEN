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

#ifndef RELAXATION_H
#define RELAXATION_H

class Relaxation; //Predeclaration of class Relaxation to include Cell.h

#include <string>
#include "../libTierces/tinyxml2.h"
#include "../Errors.h"
#include "../Tools.h"
#include "../Order1/Cell.h"

//! \class     Relaxation
//! \brief     Abstract class for Relaxations
class Relaxation
{
public:
  Relaxation();
  virtual ~Relaxation();

  //! \brief     Stiff relaxation methods
  //! \details   Call for this method computes the one of the stiff relaxation methods.
  //! \param     cell           cell to relax
  //! \param     numberPhases   number of phases
  //! \param     type           enumeration allowing to relax either state in the cell or second order half time step state
  virtual void relaxation(Cell* /*cell*/, const int& /*numberPhases*/, const double& /*dt*/, Prim /*type*/ = vecPhases) { Errors::errorMessage("relaxation not available for required relaxation"); };

  //! \brief     Return the type of the corresponding relaxation method
  virtual int getType() const { Errors::errorMessage("getType not available for required relaxation"); return 0; }

private:
};

#endif // RELAXATION_H