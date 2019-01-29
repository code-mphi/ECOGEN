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

#ifndef FLUXTHERMALEQ_H
#define FLUXTHERMALEQ_H

//! \file      FluxThermalEq.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      May 04 2018

#include "../Flux.h"
#include <iostream>

//! \class     FluxThermalEq
//! \brief     Model class for for ThermalEq system of equations (mechanical and thermal equilibrium) flux
class FluxThermalEq : public Flux
{
  public:
    FluxThermalEq();
    FluxThermalEq(const int &numberPhases);
    virtual ~FluxThermalEq();

    virtual void printFlux() const;
    virtual void addFlux(double coefA, const int &numberPhases);
    virtual void subtractFlux(double coefA, const int &numberPhases);
    virtual void multiply(double scalar, const int &numberPhases);
    virtual void setBufferFlux(Cell &cell, const int &numberPhases);
    virtual void buildCons(Phase **phases, const int &numberPhases, Mixture *mixture);
    virtual void buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases);
    virtual void setToZero(const int &numberPhases);
    virtual void setToZeroBufferFlux(const int &numberPhases);
    virtual void addNonCons(double coefA, const Cell *cell, const int &numberPhases) {};
    virtual void subtractNonCons(double coefA, const Cell *cell, const int &numberPhases) {};
    virtual void correctionEnergy(Cell *cell, const int &numberPhases, Prim type = vecPhases) const {};

    virtual void integrateSourceTermsHeating(Cell *cell, const double &dt, const int &numberPhases, const double &q);

    virtual void addTuyere1D(Coord &normal, double const &surface, Cell *cell, const int &numberPhases){};
    virtual void subtractTuyere1D(Coord &normal, double const &surface, Cell *cell, const int &numberPhases){};

    // Accessors
    //----------
    virtual double getMasse(const int &numPhase) const;
    virtual Coord getQdm() const;
    virtual double getEnergyMix() const;
    virtual void setCons(const Flux *cons, const int &numberPhases);

protected:
    double *m_masse;          //!< mass array
    Coord m_qdm;              //!< momentum array
    double m_energMixture;    //!< mixture energy

  private:

    friend class ModThermalEq;
    // To modify if needed, example: to add a class PATViscosity, add friend class PATViscosity.

};

extern FluxThermalEq *fluxBufferThermalEq;
extern FluxThermalEq *sourceConsThermEq;

#endif // FLUXTHERMALEQ_H
