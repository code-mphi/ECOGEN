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

#ifndef FLUXEULERHOMOGENEOUS_H
#define FLUXEULERHOMOGENEOUS_H

//! \file      FluxEulerHomogeneous.h
//! \author    F. Petitpas, K. Schmidmayer, J. Caze
//! \version   1.1
//! \date      November 18 2019

#include <iostream>
#include "../Flux.h"

class FluxEulerHomogeneous;

#include "ModEulerHomogeneous.h"

//! \class     FluxEulerHomogeneous
//! \brief     Model class for Homogeneous Euler Flux (liquid-vapor in thermodynamical equilibrium)
class FluxEulerHomogeneous : public Flux
{
  public:
    FluxEulerHomogeneous();
    FluxEulerHomogeneous(ModEulerHomogeneous *model);
    virtual ~FluxEulerHomogeneous();

    virtual void printFlux() const;
    virtual void addFlux(double coefA, const int &numberPhases);
    virtual void subtractFlux(double coefA, const int &numberPhases);
    virtual void multiply(double scalar, const int &numberPhases);
    virtual void setBufferFlux(Cell &cell, const int &numberPhases);
    virtual void buildCons(Phase **phase, const int &numberPhases, Mixture *mixture);
    virtual void buildPrim(Phase **phase, Mixture *mixture, const int &numberPhases);
    virtual void setToZero(const int &numberPhases);
    virtual void addNonCons(double coefA, const Cell *cell, const int &numberPhases){};
    virtual void subtractNonCons(double coefA, const Cell *cell, const int &numberPhases){};
    virtual void correctionEnergy(Cell *cell, const int &numberPhases, Prim type = vecPhases) const{};
    
    virtual void addTuyere1D(const Coord normal, const double surface, Cell *cell, const int &numberPhases);
    virtual void subtractTuyere1D(const Coord normal, const double surface, Cell *cell, const int &numberPhases);
    
    // Accessors
    //----------
    virtual const Coord& getQdm() const { return m_qdm; };
    virtual const double& getMasseMix() const { return m_masse; };
    virtual const double& getEnergyMix() const { return m_energ; };
    virtual void setCons(const Flux *cons, const int &numberPhases);

  protected:
    double m_masse;                   //!< mass
    Coord m_qdm;                      //!< momentum
    double m_energ;                   //!< total energy
    ModEulerHomogeneous *m_model;    //!< associated model

  private:

    friend class ModEulerHomogeneous;
    //friend class PAEHuler;  // To modify if needed, example: to add a class PAEHViscosity, add friend class PAEHViscosity.

};

#endif // FLUXEULERHOMOGENEOUS_H


