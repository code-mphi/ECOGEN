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

#ifndef MODUEQ_H
#define MODUEQ_H

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "MixUEq.h"

class ModUEq;

#include "FluxUEq.h"

//! \class     ModUEq
//! \brief     Model class for the velocity-equilibrium system of equations
class ModUEq : public Model
{
  public:
    //! \brief     UEq model constructor
    //! \param     numberTransports    number of additional transport equations
    //! \param     numberPhases        number of phases
    ModUEq(int& numberTransports, const int& numberPhases);
    //! \brief     Generic model constructor (used by derived classes)
    //! \param     name                 model name
    //! \param     numberTransports     number of additional transport equations
    ModUEq(const std::string& name, const int& numberTransports);
    virtual ~ModUEq();

    virtual void allocateCons(Flux** cons, const int& numberPhases);
    virtual void allocatePhase(Phase** phase);
    virtual void allocateMixture(Mixture** mixture);

    //! \details    Complete multiphase state from volume fractions, pressures, densities and velocity
    virtual void fulfillState(Phase** phases, Mixture* mixture, const int& numberPhases, Prim type = vecPhases);

    //Hydrodynamic Riemann solvers
    //----------------------------
    virtual void solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const int& numberPhases, const double& dxLeft, const double& dxRight, double& dtMax, double& massflow, double& powerFlux) const; // Riemann between two computed cells
    virtual void solveRiemannWall(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax) const; // Riemann between left cell and wall
    virtual void solveRiemannInflow(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax, const double m0, const double* ak0, const double* rhok0, const double* pk0, double& massflow, double& powerFlux) const; // Riemann for inflow (injection)
    virtual void solveRiemannTank(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax, const double* ak0, const double* rhok0, const double& p0, const double& /*T0*/, double& massflow, double& powerFlux) const; // Riemann for tank
    virtual void solveRiemannOutflow(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax, const double p0, double& massflow, double& powerFlux) const; // Riemann for outflow with imposed pressure

    //Transports Riemann solvers
    //--------------------------
    virtual void solveRiemannTransportIntern(Cell& cellLeft, Cell& cellRight, const int& numberTransports);
    virtual void solveRiemannTransportWall(const int& numberTransports);
    virtual void solveRiemannTransportInflow(Cell& cellLeft, const int& numberTransports, double* valueTransports);
    virtual void solveRiemannTransportTank(Cell& cellLeft, const int& numberTransports, double* valueTransports);
    virtual void solveRiemannTransportOutflow(Cell& cellLeft, const int& numberTransports, double* valueTransport);

    virtual void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const;

    //Accessors
    //---------
    virtual const double& getSM();
    virtual const Coord& getVelocity(const Cell* cell) const { return cell->getMixture()->getVelocity(); };
    virtual Coord& getVelocity(Cell* cell) { return cell->getMixture()->getVelocity(); };

    virtual const std::string& whoAmI() const { return m_name; };

  private:
    static const std::string NAME;

    friend class FluxUEq;
};

#endif // MODUEQ_H
