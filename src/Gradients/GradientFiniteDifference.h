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

#ifndef GRADIENTFINITEDIFFERENCE_H
#define GRADIENTFINITEDIFFERENCE_H

#include "Gradient.h"

//! \class     GradientFiniteDifference
//! \brief     Class for the finite difference like gradient method (only working on Cartesian mesh with/without AMR)

class GradientFiniteDifference : public Gradient
{
  public:

    GradientFiniteDifference();
    virtual ~GradientFiniteDifference();

    //! \brief  Compute gradients (temperature, velocity, density) of a cell
    //! \param  grads          Array of desired gradients, e.g. for temperature each component represents phase temperature and for velocity each component represents the gradient of a velocity component (grad(u), grad(v), grad(w))
    //! \param  nameVariables  Name of the variable for which the gradient is calculated
    //! \param  numPhases      Phases numbers
    void computeGradients(Cell* cell, std::vector<Coord>& grads, const std::vector<Variable>& nameVariables, const std::vector<int>& numPhases);

  protected:
};

#endif // GRADIENTFINITEDIFFERENCE_H