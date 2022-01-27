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

#ifndef GRADIENTGREENGAUSS_H
#define GRADIENTGREENGAUSS_H

#include "Gradient.h"

//! \class     GradientGreenGauss
//! \brief     Class for the Green-Gauss gradient method working on all mesh types
class GradientGreenGauss : public Gradient
{
  public:

    GradientGreenGauss();
    virtual ~GradientGreenGauss();

    //! \brief Add the contribution of the interface to build the cell gradient
    //! \param grad           Gradient to add the interface contribution
    //! \param cellInterface  Cell interface of interest used to set outward normal
    //! \param cell           Cell whose gradient must be calculated (here to access its element)
    //! \param faceValue      Value of the face interpolated using barycenter
    void addGradInterface(Coord &grad, CellInterface &cellInterface, Cell &cell, double const& faceValue);

    //! \brief  Compute gradient coord (only used to compute density gradient)
    //! \param  nameVariables  Name of the variable for which the gradient is calculated
    //! \param  numPhases      Phases number's
    Coord computeGradient(Cell* cell, Variable nameVariable, int numPhase = -1);

    //! \brief  Compute gradients (temperature, velocity, density) of a cell
    //! \param  grads          Array of desired gradients, e.g. for temperature each component represents phase temperature and for velocity each component represents the gradient of a velocity component (grad(u), grad(v), grad(w))
    //! \param  nameVariables  Name of the variable for which the gradient is calculated
    //! \param  numPhases      Phases number's
    virtual void computeGradient(Cell* cell, std::vector<Coord>& grads, std::vector<Variable>& nameVariables, std::vector<int>& numPhases);

  protected:
};

#endif // GRADIENTGRADIENTGREENGAUSS_H