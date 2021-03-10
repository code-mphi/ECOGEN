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

#ifndef MODELE_H
#define MODELE_H

class Model; //Predeclaration of class Model to include Flux.h

#include "Flux.h"
#include "../Maths/Coord.h"
#include "../Errors.h"
#include "../Relaxations/Relaxation.h"

//! \class     Model
//! \brief     Abstract class for mathematical flow models
class Model
{
  public:
    //! \brief     Generic model constructor
    //! \param     name                 model name
    //! \param     numberTransports     number of additional transport equations
    Model(const std::string& name, const int& numberTransports);
    virtual ~Model();

    //! \brief     Allocate conservative variable arrays
    //! \param     cons           conservative variable array to allocate
    //! \param     numberPhases   number of phases
    virtual void allocateCons(Flux** /*cons*/, const int& /*numberPhases*/) { Errors::errorMessage("allocateCons not available for required model"); };
    //! \brief     Instanciate phase variable
    //! \param     phase          phase to instanciate
    virtual void allocatePhase(Phase** /*phase*/) { Errors::errorMessage("allocatePhase not available for required model"); };
    //! \brief     Instanciate mixture variable
    //! \param     mixture        mixture to instanciate
    virtual void allocateMixture(Mixture** /*mixture*/) { Errors::errorMessage("allocateMixture not available for required model"); };
    //! \brief     Associate equations of state
    //! \param     cell           original cell for equation of state linking
    //! \param     numberPhases   number of phases
    void allocateEos(Cell& cell, const int& numberPhases) const;

    //! \brief     Complete a thermodynamics state frome minimum variables depending on the model
    //! \param     phases         phases array variables
    //! \param     mixture        mixture variables
    //! \param     numberPhases   number of phases
    virtual void fulfillState(Phase** /*phases*/, Mixture* /*mixture*/, const int& /*numberPhases*/, Prim /*type*/ = vecPhases) { Errors::errorMessage("fulfillState not available for required model"); };

    //Hydrodynamic Riemann solvers
    //----------------------------
    //! \brief     Cell to cell Riemann solver 
    //! \param     cellLeft          left cell
    //! \param     cellRight         right cell
    //! \param     numberPhases      number of phases
    //! \param     dxLeft            left characteristic lenght
    //! \param     dxRight           right characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     massflow          massflow through cell interface
    //! \param     powerFlux         power flux through cell interface
    virtual void solveRiemannIntern(Cell& /*cellLeft*/, Cell& /*cellRight*/, const int& /*numberPhases*/, const double& /*dxLeft*/, const double& /*dxRight*/, double& /*dtMax*/, double& /*massflow*/, double& /*powerFlux*/) const { Errors::errorMessage("solveRiemannIntern not available for required model"); };
    //! \brief     Wall half Riemann solver 
    //! \param     cellLeft          left cell
    //! \param     numberPhases      number of phases
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    virtual void solveRiemannWall(Cell& /*cellLeft*/, const int& /*numberPhases*/, const double& /*dxLeft*/, double& /*dtMax*/) const { Errors::errorMessage("solveRiemannWall not available for required model"); };
    //! \brief     Inflow (injection) half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     numberPhases      number of phases
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     m0                specific mass flow rate (kg/s/m�)
    //! \param     ak0               volume fraction array of injected fluids
    //! \param     rhok0             density array of injected fluids
    //! \param     pk0               pressure array of injected fluids
    //! \param     massflow          massflow through cell interface
    //! \param     powerFlux         power flux through cell interface
    virtual void solveRiemannInflow(Cell& /*cellLeft*/, const int& /*numberPhases*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*m0*/, const double* /*ak0*/, const double* /*rhok0*/, const double* /*pk0*/, double& /*massflow*/, double& /*powerFlux*/) const { Errors::errorMessage("solveRiemannInflow not available for required model"); };
    //! \brief     Subsonic injection half Riemann solver (Euler specific)
    //! \param     cellLeft          left cell
    //! \param     numberPhases      number of phases
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     m0                specific mass flow rate (kg.s-1.m-2)
    //! \param     T0                temperature of injected fluid
    //! \param     massflow          massflow through cell interface
    //! \param     powerFlux         power flux through cell interface
    virtual void solveRiemannSubInj(Cell& /*cellLeft*/, const int& /*numberPhases*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*m0*/, const double /*T0*/, double& /*massflow*/, double& /*powerFlux*/) const { Errors::errorMessage("solveRiemannSubInj not available for required model"); };
    //! \brief     Tank half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     numberPhases      number of phases
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     ak0               volume fraction array of fluids in tank
    //! \param     rhok0             density array of fluids in tank
    //! \param     pk0               pressure array of fluids in tank
    //! \param     massflow          massflow through cell interface
    //! \param     powerFlux         power flux through cell interface
    virtual void solveRiemannTank(Cell& /*cellLeft*/, const int& /*numberPhases*/, const double& /*dxLeft*/, double& /*dtMax*/, const double* /*ak0*/, const double* /*rhok0*/, const double& /*p0*/, const double& /*T0*/, double& /*massflow*/, double& /*powerFlux*/) const { Errors::errorMessage("solveRiemannTank not available for required model"); };
    //! \brief     Outflow half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     numberPhases      number of phases
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     p0                external pressure
    //! \param     massflow          massflow through cell interface
    //! \param     powerFlux         power flux through cell interface
    virtual void solveRiemannOutflow(Cell& /*cellLeft*/, const int& /*numberPhases*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*p0*/, double& /*massflow*/, double& /*powerFlux*/) const { Errors::errorMessage("solveRiemannOutflow not available for required model"); };

    //Transports Riemann solvers
    //--------------------------
    //! \brief     Cell to cell Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     cellRight         right cell
    //! \param     numberTransports  number of transports
    virtual void solveRiemannTransportIntern(Cell& /*cellLeft*/, Cell& /*cellRight*/, const int& /*numberTransports*/) { Errors::errorMessage("solveRiemannTransportIntern not available for required model"); };
    //! \brief     Wall half Riemann solver for transport equations
    //! \param     numberTransports  number of transports
    virtual void solveRiemannTransportWall(const int& /*numberTransports*/) { Errors::errorMessage("solveRiemannTransportWall not available for required model"); };
    //! \brief     Flow injection half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     numberTransports  number of transports
    //! \param     valueTransports   array of transport quantities injected
    virtual void solveRiemannTransportInflow(Cell& /*cellLeft*/, const int& /*numberTransports*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportInflow not available for required model"); };
    //! \brief     Tank half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     numberTransports  number of transports
    //! \param     valueTransports   array of transport quantities in tank
    virtual void solveRiemannTransportTank(Cell& /*cellLeft*/, const int& /*numberTransports*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportTank not available for required model"); };
    //! \brief     Outflow half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     numberTransports  number of transports
    //! \param     valueTransports   array of external transport quantities
    virtual void solveRiemannTransportOutflow(Cell& /*cellLeft*/, const int& /*numberTransports*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportOutflow not available for required model"); };

    //! \brief     Flux reverse projection in the absolute cartesian coordinate system
    //! \param     normal            normal vector associated to the cell interface
    //! \param     tangent           tangent vector associated to the cell interface
    //! \param     binormal          binormal vector associated to the cell interface
    virtual void reverseProjection(const Coord /*normal*/, const Coord /*tangent*/, const Coord /*binormal*/) const { Errors::errorMessage("reverseProjection not available for required model"); };

	//Relaxations
	//-----------
	void relaxations(Cell* cell, const int& numberPhases, const double& dt, Prim type = vecPhases) const;

    //Accessors
    //---------
    //! \brief     Return the local fluid velocity
    //! \return    the velocity solution of the local Riemann problem
    virtual const double& getSM() { Errors::errorMessage("getSM not available for required model"); return Errors::defaultDouble; };
    //! \brief     Return the fluid velocity of the corresponding cell
    //! \param     cell       pointer to corresponding cell
    //! \return    velocity
    virtual const Coord& getVelocity(const Cell* /*cell*/) const { Errors::errorMessage("getVelocity not available for required model"); return Coord::defaultCoord; };
    virtual Coord& getVelocity(Cell* /*cell*/) { Errors::errorMessage("getVelocity not available for required model"); return Coord::defaultCoordNonConst; };

    std::vector<Relaxation*>* getRelaxations() { return &m_relaxations; };
    
    void printInfo() const;
    virtual const std::string& whoAmI() const { return Errors::defaultString; };

  protected:
    std::string m_name;                        //!< Name of the required model
    std::vector<Relaxation *> m_relaxations;   //!< Vector of relaxation procedure
    
  private:
};

#endif // MODELE_H
