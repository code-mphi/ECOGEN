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

#ifndef MODKAPILA_H
#define MODKAPILA_H

//! \file      ModKapila.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "MixKapila.h"

class ModKapila;

#include "FluxKapila.h"

//! \class     ModKapila
//! \brief     Model class for mechanical equilibrium multiphase flows
class ModKapila : public Model
{
  public:
    //! \brief     Kapila model constructor
    //! \param     numberTransports    number of additional transport equations
    //! \param     numberPhases        number of phases
    ModKapila(int &numberTransports, const int &numberPhases);
    virtual ~ModKapila();

    virtual void allocateCons(Flux **cons, const int &numberPhases);
    virtual void allocatePhase(Phase **phase);
    virtual void allocateMixture(Mixture **mixture);

    //! \details    Complete multiphase mechanical equilibrium state from volume fractions, pressure, densities, velocity
    virtual void fulfillState(Phase **phases, Mixture *mixture, const int &numberPhases, Prim type = vecPhases);

    //Hydrodynamic Riemann solvers
    //----------------------------
    virtual void solveRiemannIntern(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const; // Riemann between two computed cells
    virtual void solveRiemannWall(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax) const; // Riemann between left cell and wall
    virtual void solveRiemannInflow(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax, const double m0, const double *ak0, const double *rhok0, const double *pk0) const; // Riemann for inflow (injection)
    virtual void solveRiemannTank(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax, const double *ak0, const double *rhok0, const double &p0, const double &T0) const; // Riemann for tank
    virtual void solveRiemannOutflow(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax, const double p0, double *debitSurf) const; // Riemann for outflow with imposed pressure

    //Transports Riemann solvers
    //--------------------------
    virtual void solveRiemannTransportIntern(Cell &cellLeft, Cell &cellRight, const int &numberTransports);
    virtual void solveRiemannTransportWall(const int &numberTransports);
    virtual void solveRiemannTransportInflow(Cell &cellLeft, const int &numberTransports, double *valueTransports);
    virtual void solveRiemannTransportTank(Cell &cellLeft, const int &numberTransports, double *valueTransports);
    virtual void solveRiemannTransportOutflow(Cell &cellLeft, const int &numberTransports, double *valueTransport);

    virtual void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const;

    //Accessors
    //---------
    virtual const double& getSM();
    virtual const Coord& getVelocity(const Cell *cell) const { return cell->getMixture()->getVelocity(); };
    virtual Coord& getVelocity(Cell *cell) { return cell->getMixture()->getVelocity(); };

    virtual const std::string& whoAmI() const { return m_name; };

  protected:
  
  private:
    static const std::string NAME;

    friend class FluxKapila;
};

#endif // MODKAPILA_H
