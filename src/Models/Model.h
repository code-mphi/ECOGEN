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
#include "GradPhase.h"
#include "GradMixture.h"

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
    //! \brief     Instanciate fluid phase variable
    //! \param     phase          phase to instanciate
    virtual void allocatePhase(Phase** /*phase*/) { Errors::errorMessage("allocatePhase not available for required model"); };
    //! \brief     Instanciate solid phase variable
    //! \param     phase          phase to instanciate
    virtual void allocatePhaseSolid(Phase** /*phase*/) {};
    //! \brief     Instanciate gradient phase variable
    //! \param     phase          phase to instanciate
    virtual void allocatePhaseGradient(GradPhase** /*phase*/) { Errors::errorMessage("allocatePhaseGradient not available for required model"); };
    //! \brief     Instanciate gradient solid-phase variable
    //! \param     phase          phase to instanciate
    virtual void allocatePhaseSolidGradient(GradPhase** /*phase*/) { Errors::errorMessage("allocatePhaseSolidGradient not available for required model"); };
    //! \brief     Instanciate mixture variable
    //! \param     mixture        mixture to instanciate
    virtual void allocateMixture(Mixture** /*mixture*/) { Errors::errorMessage("allocateMixture not available for required model"); };
    //! \brief     Instanciate gradient mixture variable
    //! \param     mixture        mixture to instanciate
    virtual void allocateMixtureGradient(GradMixture** /*mixture*/) { Errors::errorMessage("allocateMixtureGradient not available for required model"); };
    //! \brief     Associate equations of state
    //! \param     cell           original cell for equation of state linking
    void allocateEos(Cell& cell) const;
    //! \brief     Initialize the theoritical critical pressure of the fluid (only required for PTMu relax)
    //! \param     cell           cell to get the eos
    void initializeRelaxation(Cell* cell) const;

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
    //! \brief     Cell to cell Riemann solver + compute fluxBuffMRF for MRF interface
    //! \param     cellLeft          left cell
    //! \param     cellRight         right cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dxRight           right characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     omega             rotating velocity
    //! \param     normal            face normal
    //! \param     tangent           face tangent
    //! \param     binormal          face binormal
    //! \param     position          face position
    virtual void solveRiemannInternMRF(Cell& /*cellLeft*/, Cell& /*cellRight*/, const double& /*dxLeft*/, const double& /*dxRight*/, double& /*dtMax*/, const Coord& /*omega*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/, const Coord& /*position*/) const { Errors::errorMessage("solveRiemannInternMRF not available for required model"); };
    //! \brief     Wall half Riemann solver 
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannWall(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannWall not available for required model"); };
    //! \brief     Piston half Riemann solver 
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     boundData         boundary dataset used for output
    //! \param     uPiston           piston velocity
    virtual void solveRiemannPiston(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, std::vector<double>& /*boundData*/, const double& /*uPiston*/) const { Errors::errorMessage("solveRiemannPiston not available for required model"); };
    //! \brief     Inlet tank half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     ak0               volume fraction array of fluids in tank
    //! \param     rhok0             density array of fluids in tank
    //! \param     pk0               pressure array of fluids in tank
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannInletTank(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double* /*ak0*/, const double* /*rhok0*/, const double& /*p0*/, const double& /*T0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannInletTank not available for required model"); };
    //! \brief     Inlet injection using stagnation state half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     m0                specific mass flow rate (kg.s-1.m-2)
    //! \param     ak0               volume fraction array of injected fluids
    //! \param     rhok0             density array of injected fluids
    //! \param     pk0               pressure array of injected fluids
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannInletInjStagState(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*m0*/, const double* /*ak0*/, const double* /*rhok0*/, const double* /*pk0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannInletInjStagState not available for required model"); };
    //! \brief     Inlet injection using temperature half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     m0                specific mass flow rate (kg.s-1.m-2)
    //! \param     Tk0               temperature of injected fluid (same for both phases)
    //! \param     ak0               volume fraction of each injected phase
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannInletInjTemp(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*m0*/, const double* /*Tk0*/, const double* /*ak0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannInletInjTemp not available for required model"); };
    //! \brief     Outlet at imposed pressure half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     p0                external pressure
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannOutletPressure(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*p0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannOutletPressure not available for required model"); };
    //! \brief     Outlet at imposed massflow half Riemann solver
    //! \param     cellLeft          left cell
    //! \param     dxLeft            left characteristic lenght
    //! \param     dtMax             maximum explicit time step
    //! \param     m0                specific mass flow rate (kg.s-1.m-2)
    //! \param     boundData         boundary dataset used for output
    virtual void solveRiemannOutletMassflow(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, const double /*m0*/, std::vector<double>& /*boundData*/) const { Errors::errorMessage("solveRiemannOutletMassflow not available for required model"); };
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
    //! \brief     Piston half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     uPiston           piston velocity
    virtual void solveRiemannTransportPiston(Cell& /*cellLeft*/, double /*uPiston*/) { Errors::errorMessage("solveRiemannTransportPiston not available for required model"); };
    //! \brief     Inlet injection using stagnation state half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     valueTransports   array of transport quantities injected
    virtual void solveRiemannTransportInletInjStagState(Cell& /*cellLeft*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportInletInjStagState not available for required model"); };
    //! \brief     Inlet tank half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     valueTransports   array of transport quantities in tank
    virtual void solveRiemannTransportInletTank(Cell& /*cellLeft*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportInletTank not available for required model"); };
    //! \brief     Outlet at imposed pressure half Riemann solver for transport equations
    //! \param     cellLeft          left cell
    //! \param     valueTransports   array of external transport quantities
    virtual void solveRiemannTransportOutletPressure(Cell& /*cellLeft*/, double* /*valueTransports*/) { Errors::errorMessage("solveRiemannTransportOutletPressure not available for required model"); };

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

    //Moving Reference Frame
    //----------------------
    virtual void addNonConsMrfFlux(Phase** /*phases*/) { Errors::errorMessage("addNonConsMrfFlux not available for required model"); };
    virtual void reverseProjectionMrfFlux(const Coord /*normal*/, const Coord /*tangent*/, const Coord /*binormal*/) const { Errors::errorMessage("reverseProjectionMrfFlux not available for required model"); };

    //Compaction
    //----------
    //! \brief    Read compaction parameters
    //! \param    elementCompaction  XML element
    //! \param    fileName           string name of readed XML file
    //! \param    k                  phase number
    virtual void readCompactionParameters(tinyxml2::XMLElement* /*elementCompaction*/, std::string /*fileName*/, const int& /*k*/) { Errors::errorMessage("readCompactionParameters not available for required model"); };

    //! \brief    Return if compaction is activated or not
    //! \param    k               phase number
    virtual bool compaction(const int& /*k*/) const { Errors::errorMessage("compaction not available for required model"); return false; };

    //! \brief    Return the compaction reference density of phase k
    //! \param    k               phase number
    virtual const double& getDensityReference(const int& /*k*/) const { Errors::errorMessage("getDensityReference not available for required model"); return Errors::defaultDouble; };

    //! \brief    Return the compaction energy of the fluid/solid
    //! \param    k               phase number
    //! \param    alpha           volume fraction (\f$ \alpha \f$)
    //! \param    density         density (\f$ \rho \f$)
    //! \param    lambda          compaction variable (\f$ \lambda \f$)
    //!return    \f$ \kappa \f$ (compaction energy)
    virtual double computeEnergyCompaction(const int& /*k*/, const double& /*alpha*/, const double& /*density*/, const double& /*lambda*/) const { Errors::errorMessage("computeEnergyCompaction not available for required model"); return Errors::defaultDouble; };

    //! \brief  Compute the compaction pressure of the fluid/solid for relaxation
    //! \param    k               phase number
    //! \param    alpha           volume fraction (\f$ \alpha \f$)
    //! \param    density         density (\f$ \rho \f$)
    //! \param    lambda          compaction variable (\f$ \lambda \f$)
    //!return    \f$ p_c \f$ (compaction pressure)
    virtual double computeCompactionPressure(const int& /*k*/, const double& /*alpha*/, const double& /*density*/, const double& /*lambda*/) const { Errors::errorMessage("computeCompactionPressure not available for required model"); return Errors::defaultDouble; };

    //! \brief  Compute the derivative of the compaction function f with respect to xi
    //! \param    k               phase number
    //! \param    density         density (\f$ \rho \f$)
    //! \param    lambda          compaction variable (\f$ \lambda \f$)
    //!return    \f$ \partial f / \partial \xi \f$
    virtual double computeDerivativeCompactionFunctionF(const int& /*k*/, const double& /*density*/, const double& /*lambda*/) const { Errors::errorMessage("computeDerivativeCompactionFunctionF not available for required model"); return Errors::defaultDouble; };

    //! \brief  Compute the compaction plasticity term of the fluid/solid for relaxation
    //! \param    k               phase number
    //! \param    alpha           volume fraction (\f$ \alpha \f$)
    //! \param    density         density (\f$ \rho \f$)
    //! \param    lambda          compaction variable (\f$ \lambda \f$)
    //! \param    pc              compaction pressure (\f$ p_c \f$)
    //!return    \f$ \dlambda \f$ (compaction plasticity)
    virtual double computeCompactionPlasticity(const int& /*k*/, const double& /*alpha*/, const double& /*density*/, const double& /*lambda*/, const double& /*pc*/) const { Errors::errorMessage("computeCompactionPlasticity not available for required model"); return Errors::defaultDouble; };

    //Solid elasticity and plasticity
    //-------------------------------
    //! \brief    Read solid parameters
    //! \param    elementSolid       XML element
    //! \param    fileName           string name of readed XML file
    //! \param    k                  phase number
    virtual void readSolidParameters(tinyxml2::XMLElement* /*elementSolid*/, std::string /*fileName*/, const int& /*k*/) { Errors::errorMessage("readSolidParameters not available for required model"); };

    //! \brief    Return the reference density for each phase
    virtual const double* getDensityReference() const { Errors::errorMessage("getDensityReference not available for required model"); return 0; };

    //! \brief    Return the elastic parameter a of the one-parameter model for each phase
    virtual const double* getElasticParameterA() const { Errors::errorMessage("getElasticParameterA not available for required model"); return 0; };

    //! \brief    Return the shear modulus for each phase
    virtual const double* getShearModulus() const { Errors::errorMessage("getShearModulus not available for required model"); return 0; };

    //! \brief    Return the limit of elasticity for each phase
    virtual const double* getElasticityLimit() const { Errors::errorMessage("getElasticityLimit not available for required model"); return 0; };

    //! \brief    Compute the elastic energy of the solid
    //! \param    k               phase number
    //! \param    cobase          cobase tensor (\f$ e^\beta \f$)
    virtual double computeElasticEnergy(const int& /*k*/, const Tensor& /*cobase*/) const { Errors::errorMessage("computeElasticEnergy not available for required model"); return Errors::defaultDouble; };

    //! \brief    Compute the elastic energy and stress tensor of the solid
    //! \param    k               phase number
    //! \param    cobase          cobase tensor (\f$ e^\beta \f$)
    //! \param    pressure        pressure (\f$ p \f$)
    //! \param    density         density (\f$ \rho \f$)
    //! \param    elasticEnergy   elastic energy (\f$ e^e \f$)
    //! \param    stressTensor    stress tensor (\f$ \sigma \f$)
    virtual void computeElasticEnergyAndStressTensor(const int& /*k*/, const Tensor& /*cobase*/, const double& /*pressure*/, const double& /*density*/,
                                                       double& /*elasticEnergy*/, Tensor& /*stressTensor*/) const { Errors::errorMessage("computeElasticEnergyAndStressTensor not available for required model"); };

    //! \brief    Compute the square longitudinal wave speed of the solid
    //! \param    k               phase number
    //! \param    phase           phase object
    //! \return   square longitudinal wave speed
    virtual double computeSquareLongitudinalWaveSpeed(const int& /*k*/, const Phase& /*phase*/) { Errors::errorMessage("computeSquareLongitudinalWaveSpeed not available for required model"); return Errors::defaultDouble; };

    //Accessors
    //---------
    //! \brief  Select a specific scalar variable
    //! \param  phases         phases array variables
    //! \param  mixture        mixture variables
    //! \param  vecTransports  vector of transport variables
    //! \param  nameVariables  Name of the variable to select
    //! \param  numPhases      Phases number's
    virtual double selectScalar(Phase** /*phases*/, Mixture* /*mixture*/, Transport* /*transports*/, Variable /*nameVariable*/, int /*num*/ = 0) const { Errors::errorMessage("selectScalar not available for required model"); return Errors::defaultDouble; };
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
    virtual void setMachRefMin(const double& /*machRefMin*/) { Errors::errorMessage("setMachRefMin not available for required model"); };

  protected:
    std::string m_name;                        //!< Name of the required model
    std::vector<Relaxation *> m_relaxations;   //!< Vector of relaxation procedure

    bool m_lowMach;                             //!< Low-Mach preconditioning (default: false)
    double m_machRefMin;                        //!< Minimum Mach number limit for L-M preconditionning used when local Mach number is below this value (default: 0.01) 
    bool m_smoothCrossSection1d;                //!< 1D geometry with smooth cross section variation (default: false)
  private:
};

#endif // MODELE_H
