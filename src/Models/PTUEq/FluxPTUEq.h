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

#ifndef FLUXPTUEQ_H
#define FLUXPTUEQ_H

#include "../Flux.h"
#include <iostream>

//! \class     FluxPTUEq
//! \brief     Flux class for pressure-temperature-velocity (mechanical and thermal equilibrium) system of equations
class FluxPTUEq : public Flux
{
  public:
    FluxPTUEq(const int& numberPhases);
    virtual ~FluxPTUEq();

    virtual void printFlux() const;
    virtual void addFlux(double coefA, const int& numberPhases);
    virtual void addFlux(Flux* flux, const int& numberPhases);
    virtual void subtractFlux(double coefA, const int& numberPhases);
    virtual void multiply(double scalar, const int& numberPhases);
    virtual void setBufferFlux(Cell& cell, const int& numberPhases);
    virtual void buildCons(Phase** phases, const int& numberPhases, Mixture* mixture);
    virtual void buildPrim(Phase** phases, Mixture* mixture, const int& numberPhases);
    virtual void setToZero(const int& numberPhases);
    virtual void setToZeroBufferFlux(const int& numberPhases);
    virtual void addNonCons(double /*coefA*/, const Cell* /*cell*/, const int& /*numberPhases*/) {};
    virtual void subtractNonCons(double /*coefA*/, const Cell* /*cell*/, const int& /*numberPhases*/) {};

    virtual void prepSourceTermsHeating(const int& numberPhases, const double& q);

    // Accessors
    //----------
    virtual const double& getMasse(const int& numPhase) const { return m_masse[numPhase]; };
    virtual const Coord& getQdm() const { return m_qdm; };
    virtual const double& getEnergyMix() const { return m_energMixture; };
    virtual void setCons(const Flux* cons, const int& numberPhases);

  protected:
    double* m_masse;          //!< mass array
    Coord m_qdm;              //!< momentum array
    double m_energMixture;    //!< mixture energy

  private:
    friend class ModPTUEq;
    // To modify if needed, example: To add a class APPTUEqViscosity, add friend class APPTUEqViscosity.
};

#endif // FLUXPTUEQ_H
