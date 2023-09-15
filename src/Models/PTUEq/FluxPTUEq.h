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

//! \class     FluxPTUEq
//! \brief     Flux class for pressure-temperature-velocity (mechanical and thermal equilibrium) system of equations
class FluxPTUEq : public Flux
{
  public:
    FluxPTUEq(const int& numbPhases);
    virtual ~FluxPTUEq();

    virtual void printFlux() const;
    virtual void addFlux(double coefA);
    virtual void addFlux(Flux* flux);
    virtual void subtractFlux(double coefA);
    virtual void multiply(double scalar);
    virtual void setBufferFlux(Cell& cell);
    virtual void buildCons(Phase** phases, Mixture* mixture);
    virtual void buildPrim(Phase** phases, Mixture* mixture);
    virtual void setToZero();
    virtual void addNonCons(double /*coefA*/, const Cell* /*cell*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) {};
    virtual void subtractNonCons(double /*coefA*/, const Cell* /*cell*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) {};

    virtual void prepSourceTermsHeating(const double& q);

    // Accessors
    //----------
    virtual const double& getMass(const int& numPhase) const { return m_mass[numPhase]; };
    virtual const Coord& getMomentum() const { return m_momentum; };
    virtual const double& getEnergyMix() const { return m_energMixture; };
    virtual void setCons(const Flux* cons);

  protected:
    double* m_mass;          //!< mass array
    Coord m_momentum;        //!< momentum array
    double m_energMixture;   //!< mixture energy

  private:
    friend class ModPTUEq;
    // To modify if needed, example: To add a class APPTUEqViscosity, add friend class APPTUEqViscosity.
};

#endif // FLUXPTUEQ_H
