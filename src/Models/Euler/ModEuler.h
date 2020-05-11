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

#ifndef MODEULER_H
#define MODEULER_H

//! \file      ModEuler.h
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.1
//! \date      June 5 2019

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "FluxEuler.h"
#include "MixEuler.h"

//! \class     ModEuler
//! \brief     Model class for Euler mathematical system of equations (single phase)
class ModEuler : public Model
{
  public:
    //! \brief     Euler model constructor
    //! \param     numberTransports    number of additional transport equations
    ModEuler(const int &numberTransports);
    virtual ~ModEuler();

    virtual void allocateCons(Flux **cons, const int &numberPhases);
    virtual void allocatePhase(Phase **phase);
    virtual void allocateMixture(Mixture **mixture);

    //! \details    Complete single fluid state from pressure, density and velocity
    virtual void fulfillState(Phase **phases, Mixture *mixture, const int &numberPhases, Prim type = vecPhases);

    //Hydrodynamic Riemann solvers
    //----------------------------
    virtual void solveRiemannIntern(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const; 
    virtual void solveRiemannWall(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax) const; 
    virtual void solveRiemannInflow(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax, const double m0, const double *ak0, const double *rhok0, const double *pk0) const;
    virtual void solveRiemannTank(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax, const double *ak0, const double *rhok0, const double &p0, const double &T0) const;
    virtual void solveRiemannOutflow(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax, const double p0, double *debitSurf) const; 

    virtual void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const;

    //Accessors
    //---------
    virtual const double& getSM();
    virtual const Coord& getVelocity(const Cell *cell) const { return cell->getPhase(0)->getVelocity(); };
    virtual Coord& getVelocity(Cell *cell) { return cell->getPhase(0)->getVelocity(); };

    virtual const std::string& whoAmI() const { return m_name; };
  
  protected:

  private:
    static const std::string NAME;
};

#endif // MODEULER_H
