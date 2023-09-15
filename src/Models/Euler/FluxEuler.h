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

#ifndef FLUXEULER_H
#define FLUXEULER_H

#include "../Flux.h"

//! \class     FluxEuler
//! \brief     Model class for Euler Flux (single phase)
class FluxEuler : public Flux
{
  public:
    FluxEuler();
    virtual ~FluxEuler();

    virtual void printFlux() const;
    virtual void addFlux(double coefA);
    virtual void addFlux(Flux* flux);
    virtual void subtractFlux(double coefA);
    virtual void addFluxRotatingRegion(double coefA);
    virtual void subtractFluxRotatingRegion(double coefA);
    virtual void multiply(double scalar);
    virtual void setBufferFlux(Cell& cell);
    virtual void buildCons(Phase** phase, Mixture* /*mixture*/);
    virtual void buildPrim(Phase** phase, Mixture* /*mixture*/);
    virtual void setToZero();
    virtual void addNonCons(double /*coefA*/, const Cell* /*cell*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) {};
    virtual void subtractNonCons(double /*coefA*/, const Cell* /*cell*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) {};
    
    virtual void addFluxSmooth1D(double coefA, const Coord& normal, Cell* cell);
    virtual void substractFluxSmooth1D(double coefA, const Coord& normal, Cell* cell);

    virtual void addSymmetricTerms(Phase** phases, Mixture* /*mixture*/, const double& r, const double& v);
    virtual void prepSourceTermsGravity(const Coord& g);
    virtual void prepSourceTermsHeating(const double& q);
    virtual void prepSourceTermsMRF(Cell* cell, const Coord& omega);

    //Moving Reference Frame
    virtual void addNonConsMrfFlux(Phase** /*phases*/){};

    // Accessors
    //----------
    virtual const Coord& getMomentum() const { return m_momentum; };
    virtual const double& getMassMix() const { return m_mass; }; 
    virtual const double& getEnergyMix() const { return m_energ; };
    virtual void setCons(const Flux* cons);

  protected:
    double m_mass;                   //!< mass
    Coord m_momentum;                //!< momentum
    double m_energ;                  //!< total energy

  private:

    friend class ModEuler;
    // To modify if needed, example: to add a class APEViscosity, add friend class APEViscosity.
    friend class APEuler;
    friend class APEViscosity;
    friend class APEConductivity;
};

#endif // FLUXEULER_H


