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

#ifndef MODEULERHOMOGENEOUSOUS_H
#define MODEULERHOMOGENEOUSOUS_H

//! \file      ModEulerHomogeneousous.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 15 2018

#include "../Model.h"
#include "../../Cell.h"
#include "MixEulerHomogeneous.h"

class ModEulerHomogeneous;

#include "FluxEulerHomogeneous.h"

//! \class     ModEulerHomogeneous
//! \brief     Model class for Homogeneous Euler mathematical system of equations (velocity and thermodynamical equilibrium)
class ModEulerHomogeneous : public Model
{
  public:
    //! \brief     Homogeneous Euler model constructor
    //! \param     numberTransports    number of additional transport equations
    //! \param     liquid              fluid number for liquid phase
    //! \param     vapor               fluid number for vapor phase
    ModEulerHomogeneous(const int &numberTransports, const int liquid = 0, const int vapor = 1);
    virtual ~ModEulerHomogeneous();

    virtual void allocateCons(Flux **cons, const int &numberPhases);
    virtual void allocatePhase(Phase **phase);
    virtual void allocateMixture(Mixture **mixture);

    virtual void fulfillState(Phase **phases, Mixture *mixture, const int &numberPhases, Prim type = vecPhases);

    //Hydrodynamic Riemann solvers
    //----------------------------
    virtual void solveRiemannIntern(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const;
    virtual void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const;

    //Accessors
    //---------
    virtual double getSM();
    virtual Coord getVelocity(Cell *cell) const;
    int getLiq();
    int getVap();
    virtual std::string whoAmI() const;

  protected:

  private:
    static const std::string NAME;
	int m_liq;                                 //!< Liquid phase number for phase change
	int m_vap;                                 //!< Vapor phase number for phase change

    friend class FluxEulerHomogeneous;
};

#endif // MODEULERHOMOGENEOUSOUS_H
