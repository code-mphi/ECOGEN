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

#ifndef FLUXUEQTOTENERGY_H
#define FLUXUEQTOTENERGY_H

#include <iostream>
#include "../Flux.h"

class FluxUEqTotEnergy;

#include "ModUEqTotEnergy.h"

//! \class     FluxUEqTotEnergy
//! \brief     Flux class for the velocity-equilibrium system of equations
class FluxUEqTotEnergy : public Flux
{
  public:
    FluxUEqTotEnergy(const int& numbPhases);
    virtual ~FluxUEqTotEnergy();

    virtual void printFlux() const;
    virtual void addFlux(double coefA);
    virtual void addFlux(Flux* flux);
    virtual void subtractFlux(double coefA);
    virtual void multiply(double scalar);
    virtual void setBufferFlux(Cell& cell);
    virtual void buildCons(Phase** phases, Mixture* mixture);
    virtual void buildPrim(Phase** phases, Mixture* mixture);
    virtual void setToZero();
    virtual void addNonCons(double coefA, const Cell* cell, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/);
    virtual void subtractNonCons(double coefA, const Cell* cell, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/);

    // Accessors
    //----------
    virtual const double& getAlpha(const int& numPhase) const { return m_alpha[numPhase]; };
    virtual const double& getMass(const int& numPhase) const { return m_mass[numPhase]; };
    virtual const double& getTotEnergy(const int& numPhase) const { return m_totEnerg[numPhase]; };
    virtual const Coord& getMomentum() const { return m_momentum; };
    virtual void setCons(const Flux* cons);

  protected:
    double* m_alpha;          //!< Volume fraction array
    double* m_mass;           //!< Mass array
    double* m_totEnerg;       //!< Specific total energy array
    Coord m_momentum;         //!< Momentum array
    double* m_alphap;         //!< One of the non-conservative flux for phasic total energy

  private:
    friend class ModUEqTotEnergy;
};

#endif // FLUXUEQTOTENERGY_H
