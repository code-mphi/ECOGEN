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

#ifndef FLUXKAPILA_H
#define FLUXKAPILA_H

//! \file      FluxKapila.h
//! \author    F. Petitpas, K. Schmidmayer, J. Caze
//! \version   1.1
//! \date      November 19 2019

#include <iostream>
#include "../Flux.h"

class FluxKapila;

#include "ModKapila.h"

//! \class     FluxKapila
//! \brief     Model class for for Kapila system of equations (mechanical equilibrium) flux
class FluxKapila : public Flux
{
  public:
    FluxKapila();
    FluxKapila(ModKapila *model, const int &numberPhases); //JC//Q why some constructors don't need to take into parameter *model ?
    virtual ~FluxKapila();

    virtual void printFlux() const;
    virtual void addFlux(double coefA, const int &numberPhases);
    virtual void addFlux(Flux* flux, const int& numberPhases);
    virtual void subtractFlux(double coefA, const int &numberPhases);
    virtual void multiply(double scalar, const int &numberPhases);
    virtual void setBufferFlux(Cell &cell, const int &numberPhases);
    virtual void buildCons(Phase **phases, const int &numberPhases, Mixture *mixture);
    virtual void buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases);
    virtual void setToZero(const int &numberPhases);
    virtual void setToZeroBufferFlux(const int &numberPhases);
    virtual void addNonCons(double coefA, const Cell *cell, const int &numberPhases);
    virtual void subtractNonCons(double coefA, const Cell *cell, const int &numberPhases);
    virtual void correctionEnergy(Cell *cell, const int &numberPhases, Prim type = vecPhases) const;

    virtual void addSymmetricTerms(Phase **phases, Mixture *mixture, const int &numberPhases, const double &r, const double &v);
    virtual void prepSourceTermsGravity(Cell *cell, const double &dt, const int &numberPhases, const Coord &g);
    virtual void prepSourceTermsHeating(Cell *cell, const double &dt, const int &numberPhases, const double &q);
    virtual void prepSourceTermsMRF(Cell *cell, const double &dt, const int &numberPhases, const Coord &omega);

    // Accessors
    //----------
    virtual const double& getAlpha(const int &numPhase) const { return m_alpha[numPhase]; };
    virtual const double& getMasse(const int &numPhase) const { return m_masse[numPhase]; };
    virtual const double& getEnergy(const int &numPhase) const { return m_energ[numPhase]; };
    virtual const Coord& getQdm() const { return m_qdm; };
    virtual const double& getEnergyMix() const { return m_energMixture; };
    virtual void setCons(const Flux *cons, const int &numberPhases);

protected:
    double *m_alpha;          //!< volume fraction array
    double *m_masse;          //!< mass array
    double *m_energ;          //!< specific internal energy array
    Coord m_qdm;              //!< momentum array
    double m_energMixture;    //!< mixture energy
    ModKapila *m_model;       //!< associated model

  private:

    friend class ModKapila;
    // To modify if needed, example: to add a class APKViscosity, add friend class APKViscosity.
    friend class APKSurfaceTension;
    friend class APKViscosity;
    friend class APKConductivity;

};

#endif // FLUXKAPILA_H
