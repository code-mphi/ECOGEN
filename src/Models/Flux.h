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

#ifndef FLUX_H
#define FLUX_H

class Flux; //Predeclaration of class Flux to include Cell.h

#include "Phase.h"
#include "../Order1/Cell.h"
#include "../Tools.h"

//! \class     Flux
//! \brief     Abstract class for conservative variables and fluxes
class Flux
{
  public:
    Flux();
    virtual ~Flux();

    virtual void printFlux() const { Errors::errorMessage("printFlux not available for required model"); };

    //! \brief     Add flux passed in parameter to the correspond model flux
    //! \param     flux           flux to add to the current one
    //! \param     numberPhases   number of phases
    virtual void addFlux(Flux* /*flux*/, const int& /*numberPhases*/) { Errors::errorMessage("addFlux not available for required model"); };
    //! \brief     Add flux to the corresponding model flux
    //! \param     coefA          possibility to multiply the flux before adding (set 1.d0 if not needed)
    //! \param     numberPhases   number of phases
    virtual void addFlux(double /*coefA*/, const int& /*numberPhases*/){ Errors::errorMessage("addFlux not available for required model"); };
    //! \brief     Subtract flux to the corresponding model buffer flux
    //! \param     coefA          possibility to multiply the flux before subtraction (set 1.d0 if not needed)
    //! \param     numberPhases   number of phases
    virtual void subtractFlux(double /*coefA*/, const int& /*numberPhases*/){ Errors::errorMessage("subtractFlux not available for required model"); };
    //! \brief     multiply the flux of the corresponding model by a constant
    //! \param     scalar       constant
    //! \param     numberPhases   number of phases
    virtual void multiply(double /*scalar*/, const int& /*numberPhases*/){ Errors::errorMessage("multiply not available for required model"); };
    //! \brief     Temporary store the conservative variables of a given cell
    //! \details   The conservatvie variables are temporary stored in the corresponding model buffer flux
    //! \param     cell           cell used for conservative variables calculus
    //! \param     numberPhases   number of phases
    virtual void setBufferFlux(Cell& /*cell*/, const int& /*numberPhases*/){ Errors::errorMessage("setBufferFlux not available for required model"); };
    // virtual void setBufferFlux(Cell* cell, const int& numberPhases){ Errors::errorMessage("setBufferFlux not available for required model"); };
    //! \brief     Build the conservative variables for a given cell from primitive one
    //! \param     phases         Phases array used for conservative variables calculus
    //! \param     mixture        Mixture used for conservative variables calculus
    //! \param     numberPhases   number of phases
    virtual void buildCons(Phase** /*phases*/, const int& /*numberPhases*/, Mixture* /*mixture*/) { Errors::errorMessage("buildCons not available for required model"); };
    //! \brief     Build the primitive variables for a given cell from conservative one
    //! \param     phases         Phases array to fill
    //! \param     mixture        Mixture to fill
    //! \param     numberPhases   number of phases
    virtual void buildPrim(Phase** /*phases*/, Mixture* /*mixture*/, const int& /*numberPhases*/) { Errors::errorMessage("buildPrim not available for required model"); };
    //! \brief     set each attribute of the flux to zero
    //! \param     numberPhases   number of phases
    virtual void setToZero(const int& /*numberPhases*/){ Errors::errorMessage("setToZero not available for required model"); };
    //! \brief     set each attribute of the corresponding buffer flux to zero
    //! \param     numberPhases   number of phases
    virtual void setToZeroBufferFlux(const int& /*numberPhases*/) { Errors::errorMessage("setToZero not available for required model"); };
    
    //! \brief     Add non conservative term to the flux
    //! \param     coefA          possibility to multiply the non conservative term before adding (set 1.d0 if not needed)
    //! \param     cell           reference cell used to approximate the non conservative term
    //! \param     numberPhases   number of phases
    virtual void addNonCons(double /*coefA*/, const Cell* /*cell*/, const int& /*numberPhases*/){ Errors::errorMessage("addNonCons not available for required model"); };
    //! \brief     Subtract non conservative term to the flux
    //! \param     coefA          possibility to multiply the non conservative term before subtraction (set 1.d0 if not needed)
    //! \param     cell           reference cell used to approximate the non conservative term
    //! \param     numberPhases   number of phases
    virtual void subtractNonCons(double /*coefA*/, const Cell* /*cell*/, const int& /*numberPhases*/){ Errors::errorMessage("subtractNonCons not available for required model"); };
    //! \brief     Method to correct energy in non conservative models using total energy conservation
    //! \param     cell           cell to correct
    //! \param     numberPhases   number of phases
    //! \param     type           enumeration allowing to correct either state in the cell or second order half time step state
    virtual void correctionEnergy(Cell* /*cell*/, const int& /*numberPhases*/, Prim /*type*/ = vecPhases) const {};

    virtual void schemeCorrection(const int& /*numberPhases*/) const {};

    //! \brief     Add symetric terms 
    //! \param     r   radial distance of the cell from the axis of symmetry
    //! \param     v   velocity in the radial direction
    virtual void addSymmetricTerms(Phase** /*phases*/, Mixture* /*mixture*/, const int& /*numberPhases*/, const double& /*r*/, const double& /*v*/) { Errors::errorMessage("addSymmetricTerms not implemented for used model"); };

    //! \brief     Gravity source term
    virtual void prepSourceTermsGravity(const int& /*numberPhases*/, const Coord& /*g*/) { Errors::errorMessage("prepSourceTermsGravity not available for required model"); };
    //! \brief     Heating source term
    virtual void prepSourceTermsHeating(const int& /*numberPhases*/, const double& /*q*/) { Errors::errorMessage("prepSourceTermsHeating not available for required model"); };
    //! \brief     MRF source term
    virtual void prepSourceTermsMRF(Cell* /*cell*/, const int& /*numberPhases*/, const Coord& /*omega*/) { Errors::errorMessage("prepSourceTermsMRF not available for required model"); };

    virtual void addTuyere1D(const Coord& /*normal*/, const double& /*surface*/, Cell* /*cell*/, const int& /*numberPhases*/) {};
    virtual void subtractTuyere1D(const Coord& /*normal*/, const double& /*surface*/, Cell* /*cell*/, const int& /*numberPhases*/) {};

    // Accessors
    //----------
    virtual const double& getAlpha(const int& /*numPhase*/) const { return Errors::defaultDouble; };
    virtual const double& getMasse(const int& /*numPhase*/) const { return Errors::defaultDouble; };
    virtual const double& getEnergy(const int& /*numPhase*/) const { return Errors::defaultDouble; };
    virtual const Coord& getQdm() const { return Coord::defaultCoord; };
    virtual const double& getMasseMix() const { return Errors::defaultDouble; };
    virtual const double& getEnergyMix() const { return Errors::defaultDouble; };
    virtual void setCons(const Flux* /*cons*/, const int& /*numberPhases*/) { Errors::errorMessage("setCons not available for required model"); };

  protected:
    double  m_sM;     //!< Fluid velocity for intercell interfaces
  private:
};

extern std::vector<Flux*> sourceCons;
extern Flux* fluxBuff;

#endif // FLUX_H
