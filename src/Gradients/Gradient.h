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

#ifndef GRADIENT_H
#define GRADIENT_H

class Gradient; //Predeclaration of class to include following .h

#include <vector>
#include "../Maths/Coord.h"
#include "../Order1/Cell.h"

//! \class     Gradient
//! \brief     Base class for the gradient method
class Gradient
{
  public:

    Gradient();
    virtual ~Gradient();
    
    //! \brief  Compute gradient coord (only used to compute density gradient)
    //! \param  cell           Cell whose gradient must be calculated
    //! \param  nameVariables  Name of the variable for which the gradient is calculated
    //! \param  numPhase       Phase number
    virtual Coord computeGradient(Cell* /*cell*/, Variable /*nameVariable*/, int /*numPhase*/ = -1) { 
      Errors::errorMessage("computeGradient not available for required gradient method");
      return Coord::defaultCoord;
    };

    //! \brief  Compute gradients (temperature, velocity, density) of a cell
    //! \param  cell           Cell whose gradient must be calculated
    //! \param  grads          Array of desired gradients, e.g. for temperature each component represents phase temperature and for velocity each component represents the gradient of a velocity component (grad(u), grad(v), grad(w))
    //! \param  nameVariables  Name of the variable for which the gradient is calculated
    //! \param  numPhases      Phase numbers
    virtual void computeGradient(Cell* /*cell*/, std::vector<Coord>& /*grads*/, std::vector<Variable>& /*nameVariables*/, std::vector<int>& /*numPhases*/) { 
      Errors::errorMessage("computeGradient not available for required gradient method"); 
    };

    virtual const TypeGrad& getType() const { return m_type; };

  protected:
    TypeGrad m_type;
};

#endif // GRADIENT_H