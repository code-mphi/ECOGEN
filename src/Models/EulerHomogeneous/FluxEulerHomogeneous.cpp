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

//! \file      FluxEulerHomogeneous.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      December 22 2017

#include <cmath>
#include "FluxEulerHomogeneous.h"

using namespace std;

FluxEulerHomogeneous fluxBufferEulerHomogeneous;

//***********************************************************************

FluxEulerHomogeneous::FluxEulerHomogeneous() {}

//***********************************************************************

FluxEulerHomogeneous::FluxEulerHomogeneous(ModEulerHomogeneous *model) : m_model(model)
{}

//***********************************************************************

FluxEulerHomogeneous::~FluxEulerHomogeneous(){}

//***********************************************************************

void FluxEulerHomogeneous::printFlux() const
{
  cout << m_masse << " " << m_qdm.getX() << " " << m_energ << endl;
}

//***********************************************************************

void FluxEulerHomogeneous::addFlux(double coefA, const int &numberPhases)
{
    m_masse += coefA*fluxBufferEulerHomogeneous.m_masse;
    m_qdm   += coefA*fluxBufferEulerHomogeneous.m_qdm;
    m_energ += coefA*fluxBufferEulerHomogeneous.m_energ;
}

//***********************************************************************

void FluxEulerHomogeneous::subtractFlux(double coefA, const int &numberPhases)
{
    m_masse -= coefA*fluxBufferEulerHomogeneous.m_masse;
    m_qdm   -= coefA*fluxBufferEulerHomogeneous.m_qdm;
    m_energ -= coefA*fluxBufferEulerHomogeneous.m_energ;
}

//***********************************************************************

void FluxEulerHomogeneous::multiply(double scalar, const int &numberPhases)
{
    m_masse *= scalar;
    m_qdm   *= scalar;
    m_energ *= scalar;
}

//***********************************************************************

void FluxEulerHomogeneous::setBufferFlux(Cell &cell, const int &numberPhases)
{
  fluxBufferEulerHomogeneous.buildCons(cell.getPhases(), numberPhases, cell.getMixture());
}

//***********************************************************************

void FluxEulerHomogeneous::buildCons(Phase **phases, const int &numberPhases, Mixture *mixture)
{
  double energieInterne(0.);
  double rhok, alphak, ek;

  Phase *phase(0);
  m_masse = 0.;
  m_energ = 0.;

  for (int k = 0; k < numberPhases; k++)
  {
    phase = phases[k];
    //Mixture density calculus
    alphak = phase->getAlpha();
    rhok = phase->getDensity();
    m_masse += alphak * rhok;
    //Mixture specific internal energy calculus
    ek = phase->getEos()->computeEnergy(rhok, phase->getPressure());
    energieInterne += alphak*rhok*ek;
  }
  m_qdm = m_masse*mixture->getVelocity();
  m_energ = energieInterne + 0.5*m_masse*mixture->getVelocity().squaredNorm();
}

//***********************************************************************

void FluxEulerHomogeneous::buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases)
{
  double pressure, Tsat, internalEnergy;
  Phase* phase(0);
  
  int liq(m_model->m_liq), vap(m_model->m_vap);

  //Simple extractions
  mixture->setVelocity(m_qdm / m_masse);
  //Erasing small velocity variations
  if (abs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (abs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (abs(mixture->getW()) < 1.e-8) mixture->setW(0.);
  internalEnergy = m_energ / m_masse - 0.5*mixture->getVelocity().squaredNorm();
  
  //Pressure determination
  pressure = mixture->computePressure(m_masse, internalEnergy, phases, mixture, numberPhases, liq, vap);
  Tsat = mixture->computeTsat(phases[liq]->getEos(), phases[vap]->getEos(), pressure);
  phases[vap]->setDensity(phases[vap]->getEos()->computeDensitySaturation(pressure, Tsat, Tsat));
  phases[liq]->setDensity(phases[liq]->getEos()->computeDensitySaturation(pressure, Tsat, Tsat));
  phases[vap]->setAlpha((m_masse - phases[liq]->getDensity())/(phases[vap]->getDensity()- phases[liq]->getDensity()));
  phases[liq]->setAlpha(1. - phases[vap]->getAlpha());
  
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->setPressure(pressure);
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  //Mixture variables
  mixture->setPressure(pressure);
  mixture->setTemperature(Tsat);
  mixture->computeMixtureVariables(phases, numberPhases);
  //Reconstruction of total energy from total energy equation
  double totalEnergy(m_energ / m_masse);
  mixture->setTotalEnergy(totalEnergy);
}

//***********************************************************************

void FluxEulerHomogeneous::setToZero(const int &numberPhases)
{
  m_masse = 0.;
  m_qdm   = 0.;
  m_energ = 0.;
}

//***********************************************************************

void FluxEulerHomogeneous::addTuyere1D(const Coord normal, const double surface, Cell *cell, const int &numberPhases)
{
  double coef = normal.getX()*surface / cell->getElement()->getVolume();
  Phase * phase;
  phase = cell->getPhase(0);

  m_qdm.setX(m_qdm.getX() - phase->getPressure()*coef);
}
//***********************************************************************

void FluxEulerHomogeneous::subtractTuyere1D(const Coord normal, const double surface, Cell *cell, const int &numberPhases)
{
  double coef = normal.getX()*surface / cell->getElement()->getVolume();
  Phase * phase;
  phase = cell->getPhase(0);

  m_qdm.setX(m_qdm.getX() + phase->getPressure()*coef);

}

//***********************************************************************

Coord FluxEulerHomogeneous::getQdm() const
{
  return m_qdm;
}

//***********************************************************************

double FluxEulerHomogeneous::getMasseMix() const
{
  return m_masse;
}

//***********************************************************************

double FluxEulerHomogeneous::getEnergyMix() const
{
  return m_energ;
}

//***********************************************************************

void FluxEulerHomogeneous::setCons(const Flux *cons, const int &numberPhases)
{
  m_masse = cons->getMasseMix();
  m_qdm = cons->getQdm();
  m_energ = cons->getEnergyMix();
}

//***********************************************************************