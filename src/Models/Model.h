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

class Model; //Predeclaration of class to include following .h

#include "Flux.h"
#include "../Maths/Coord.h"
#include "../Errors.h"
#include "../Relaxations/Relaxation.h"

//Only boundary Riemann solvers require to extract interface data for output
//This vector allows to use default value for intern Riemann solvers
static std::vector<double> DEFAULT_VEC_INTERFACE_DATA(VarBoundary::SIZE, 0.);

//! \class     Model
//! \brief     Abstract class for mathematical flow models
class Model
{
  public:
    //! \brief     Generic model constructor
    //! \param     name                 model name
    //! \param     numbTransports       number of additional transport equations
    Model(const std::string& name, const int& numbTransports);
    virtual ~Model();

    //! \brief     Allocate conservative variable arrays
    //! \param     cons           conservative variable array to allocate
    virtual void allocateCons(Flux** /*cons*/) { Errors::errorMessage("allocateCons not available for required model"); };
    //! \brief     Instanciate phase variable
    //! \param     phase          phase to instanciate
    virtual void allocatePhase(Phase** /*phase*/) { Errors::errorMessage("allocatePhase not available for required model"); };
    //! \brief     Instanciate mixture variable
    //! \param     mixture        mixture to instanciate
    virtual void allocateMixture(Mixture** /*mixture*/) { Errors::errorMessage("allocateMixture not available for required model"); };
    //! \brief     Associate equations of state
    //! \param     cell           original cell for equation of state linking
    void allocateEos(Cell& cell) const;

    //! \brief     Complete a thermodynamics state frome minimum variables depending on the model
    //! \param     phases         phases array variables
    //! \param     mixture        mixture variables
    virtual void fulfillState(Phase** /*phases*/, Mixture* /*mixture*/) { Errors::errorMessage("fulfillState not available for required model"); };

    //! \brief     Complete some variables if necessary when restarting a simulation
    //! \param     phases         phases array variables
    //! \param     mixture        mixture variables
    virtual void fulfillStateRestart(Phase** /*phases*/, Mixture* /*mixture*/) { Errors::errorMessage("fulfillStateRestart not available for required model"); };

    //! \brief     Complete the augmented variables (such as the ones of Euler-Korteweg model)
    //! \param     cell           cell
    virtual void initializeAugmentedVariables(Cell* /*cell*/) { Errors::errorMessage("initializeAugmentedVariables not available for required model"); };

    //Hydrodynamic Riemann solvers
    //----------------------------
    //! \brief     Cell to cell Riemann solver 
    //! \param     cellLeft          left cell
    //! \param     cellRight         right cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dxRight           right characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannIntern(Cell& /*cellLeft*/, Cell& /*cellRight*/, const double& /*dxLeft*/, const double& /*dxRight*/, double& /*dtMax*/, std::vector<double>& /*boundData*/ = DEFAULT_VEC_INTERFACE_DATA) const { Errors::errorMessage("solveRiemannIntern not available for required model"); };
    //! \brief     Wall half Riemann solver 
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannWall(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannWall not available for required model"); };
    //! \brief     Inflow (injection) half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     m0                specific mass flow rate (kg/s/m�)
    //! \param     ak0               volume fraction array of injected fluids
    //! \param     rhok0             density array of injected fluids
    //! \param     pk0               pressure array of injected fluids
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannInflow(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*m0*/, const double* /*ak0*/, const double* /*rhok0*/, const double* /*pk0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannInflow not available for required model"); };
    //! \brief     Subsonic injection half Riemann solver (Euler specific)
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     m0                specific mass flow rate (kg.s-1.m-2)
    //! \param     Tk0               temperature of injected fluid (same for both phases)
    //! \param     ak0               volume fraction of each injected phase
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannSubInj(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*m0*/, const double* /*Tk0*/, const double* /*ak0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannSubInj not available for required model"); };
    //! \brief     Tank half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     ak0               volume fraction array of fluids in tank
    //! \param     rhok0             density array of fluids in tank
    //! \param     pk0               pressure array of fluids in tank
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannTank(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double* /*ak0*/, const double* /*rhok0*/, const double& /*p0*/, const double& /*T0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannTank not available for required model"); };
    //! \brief     Outflow half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     p0                external pressure
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannOutflow(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*p0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannOutflow not available for required model"); };
    //! \brief     No flux half Riemann solver (return null flux to use with 1D geometry with smooth varying cross section)
    virtual void solveRiemannNullFlux() const { Errors::errorMessage("solveRiemannNullFlux not available for required model"); };

    //Transports Riemann solvers
    //--------------------------
    //! \brief     Cell to cell Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     cellRight         right cell
    virtual void solveRiemannTransportIntern(Cell& /*cellLeft*/, Cell& /*cellRight*/) { Errors::errorMessage("solveRiemannTransportIntern not available for required model"); };
    //! \brief     Wall half Riemann solver for transport equations
    virtual void solveRiemannTransportWall() { Errors::errorMessage("solveRiemannTransportWall not available for required model"); };
    //! \brief     Flow injection half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     valueTransports   array of transport quantities injected
    virtual void solveRiemannTransportInflow(Cell& /*cellLeft*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportInflow not available for required model"); };
    //! \brief     Tank half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     valueTransports   array of transport quantities in tank
    virtual void solveRiemannTransportTank(Cell& /*cellLeft*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportTank not available for required model"); };
    //! \brief     Outflow half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     valueTransports   array of external transport quantities
    virtual void solveRiemannTransportOutflow(Cell& /*cellLeft*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportOutflow not available for required model"); };

    //! \brief     Flux reverse projection in the absolute Cartesian coordinate system
    //! \param     normal            normal vector associated to the cell interface
    //! \param     tangent           tangent vector associated to the cell interface
    //! \param     binormal          binormal vector associated to the cell interface
    virtual void reverseProjection(const Coord /*normal*/, const Coord /*tangent*/, const Coord /*binormal*/) const { Errors::errorMessage("reverseProjection not available for required model"); };

    //Relaxations
    //-----------
    void relaxations(Cell* cell, const double& dt, Prim type = vecPhases) const;

    //Low-Mach preconditioning
    //------------------------
    virtual void lowMachSoundSpeed(double& /*machRef*/, const double& /*uL*/, double& /*cL*/, const double& /*uR*/ = Errors::defaultDouble, double& /*cR*/ = Tools::uselessDouble) const { Errors::errorMessage("lowMachSoundSpeed not available for required model"); };

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

    virtual const Coord& getVectorP(const Cell* /*cell*/) const { Errors::errorMessage("getVectorP not available for required model"); return Coord::defaultCoord; };
    virtual Coord& getVectorP(Cell* /*cell*/) { Errors::errorMessage("getVectorP not available for required model"); return Coord::defaultCoordNonConst; };

    std::vector<Relaxation*>* getRelaxations() { return &m_relaxations; };
    
    void printInfo() const;
    virtual const std::string& whoAmI() const { return Errors::defaultString; };

    virtual void setSmoothCrossSection1d(const bool& /*applySmooth*/) { Errors::errorMessage("setSmoothCrossSection1d not available for required model"); };
    virtual const bool& isSmoothCrossSection1d() const { return m_smoothCrossSection1d; };

    virtual void setLowMach(const bool& /*lowMach*/) { Errors::errorMessage("setLowMach not available for required model"); };

  protected:
    std::string m_name;                        //!< Name of the required model
    std::vector<Relaxation *> m_relaxations;   //!< Vector of relaxation procedure

    bool m_lowMach;                             //!< Low-Mach preconditioning (default: false)
    bool m_smoothCrossSection1d;                //!< 1D geometry with smooth cross section variation (default: false)
  private:
};

#endif // MODELE_H
