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

#ifndef RELAXATIONP_H
#define RELAXATIONP_H

#include "Relaxation.h"

//! \class     RelaxationP
//! \brief     Pressure relaxation
class RelaxationP : public Relaxation
{
  public:
    RelaxationP();
    virtual ~RelaxationP();

    //! \brief     Return the pressure-relaxation type
    virtual int getType() const { return P; }

    //! \brief     Newton-Raphson method for the infinite pressure relaxation
    //! \details   Call of this method computes the totally relaxed pressure in a given cell.
    //! \param     pStar          initial and final pressure value
    //! \param     iteration      number of iterations for convergence of the method
    void NewtonRaphson(double& pStar, int& iteration);

    //! \brief     Compute interface pressure
    //! \details   Call for this method computes the interface pressure in a cell.
    //! \param     cell           cell
    //! \param     type           enumeration allowing to relax either state in the cell or second order half time step state
    //! \return    interface pressure
    double computeInterfacePressure(Cell* cell, Prim type = vecPhases);
};

#endif // RELAXATIONP_H


