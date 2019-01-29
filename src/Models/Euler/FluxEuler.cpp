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

//! \file      FluxEuler.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      December 22 2017

#include <cmath>
#include "FluxEuler.h"

using namespace std;

FluxEuler fluxBufferEuler;
FluxEuler sourceConsEul;

//***********************************************************************

FluxEuler::FluxEuler(){}

//***********************************************************************

FluxEuler::~FluxEuler(){}

//***********************************************************************

void FluxEuler::printFlux() const
{
  cout << m_masse << " " << m_qdm.getX() << " " << m_energ << endl;
}

//***********************************************************************

void FluxEuler::addFlux(double coefA, const int &numberPhases)
{
    m_masse += coefA*fluxBufferEuler.m_masse;
    m_qdm   += coefA*fluxBufferEuler.m_qdm;
    m_energ += coefA*fluxBufferEuler.m_energ;
}

//***********************************************************************

void FluxEuler::subtractFlux(double coefA, const int &numberPhases)
{
    m_masse -= coefA*fluxBufferEuler.m_masse;
    m_qdm   -= coefA*fluxBufferEuler.m_qdm;
    m_energ -= coefA*fluxBufferEuler.m_energ;
}

//***********************************************************************

void FluxEuler::multiply(double scalar, const int &numberPhases)
{
    m_masse *= scalar;
    m_qdm   *= scalar;
    m_energ *= scalar;
}

//***********************************************************************

void FluxEuler::setBufferFlux(Cell &cell, const int &numberPhases)
{
  fluxBufferEuler.buildCons(cell.getPhases(), numberPhases, cell.getMixture());
}

//***********************************************************************

void FluxEuler::buildCons(Phase **phases, const int &numberPhases, Mixture *mixture)
{
  Phase* phase(phases[0]);

  Eos *eos(phase->getEos());
  
  m_masse = phase->getDensity();
  m_qdm = m_masse*phase->getVelocity();
  m_energ = m_masse*phase->getTotalEnergy();
}

//***********************************************************************

void FluxEuler::buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases)
{
  double pressure(0.), energieInterne(0.), totalEnergy(0.), soundSpeed(0.), temperature(0.);
  Phase* phase(phases[0]);
  Eos *eos(phase->getEos());

  phase->setDensity(m_masse);
  phase->setVelocity(m_qdm.getX() / m_masse, m_qdm.getY() / m_masse, m_qdm.getZ() / m_masse);
  //Erasing small velocity variations
  if (abs(phase->getU()) < 1.e-8) phase->setU(0.);
  if (abs(phase->getV()) < 1.e-8) phase->setV(0.);
  if (abs(phase->getW()) < 1.e-8) phase->setW(0.);

  totalEnergy = m_energ / m_masse;
  phase->setTotalEnergy(totalEnergy);
  energieInterne = totalEnergy - 0.5*phase->getVelocity().squaredNorm();
  phase->setEnergy(energieInterne);
  pressure = eos->computePressure(m_masse, energieInterne);
  phase->setPressure(pressure);
  soundSpeed = eos->computeSoundSpeed(m_masse, pressure);
  phase->setSoundSpeed(soundSpeed);
  
  temperature = eos->computeTemperature(m_masse, pressure);
  phase->setTemperature(temperature);

}

//***********************************************************************

void FluxEuler::setToZero(const int &numberPhases)
{
  m_masse = 0.;
  m_qdm   = 0.;
  m_energ = 0.;
}

//***********************************************************************

void FluxEuler::addTuyere1D(const Coord normal, const double surface, Cell *cell, const int &numberPhases)
{
  double coef = normal.getX()*surface / cell->getElement()->getVolume();
  Phase * phase;
  phase = cell->getPhase(0);

  m_qdm.setX(m_qdm.getX() - phase->getPressure()*coef);
}

//***********************************************************************

void FluxEuler::subtractTuyere1D(const Coord normal, const double surface, Cell *cell, const int &numberPhases)
{
  double coef = normal.getX()*surface / cell->getElement()->getVolume();
  Phase * phase;
  phase = cell->getPhase(0);

  m_qdm.setX(m_qdm.getX() + phase->getPressure()*coef);

}

//***********************************************************************

void FluxEuler::integrateSourceTermsHeating(Cell *cell, const double &dt, const int &numberPhases, const double &q)
{
  sourceConsEul.setToZero(1);
  sourceConsEul.m_energ = q;

  m_energ += dt*sourceConsEul.m_energ;
}

//***********************************************************************

void FluxEuler::integrateSourceTermsMRF(Cell *cell, const double &dt, const int &numberPhases, const Coord &omega)
{
  sourceConsEul.setToZero(1);
  //Mass and velocity extraction
  double rho = cell->getPhase(0)->getDensity();
  Coord u = cell->getPhase(0)->getVelocity();

  //Coriolis acceleration
  sourceConsEul.m_qdm = -2.*rho*Coord::crossProduct(omega, u);
  //Centrifugal acceleration
  sourceConsEul.m_qdm -= rho*Coord::crossProduct(omega, Coord::crossProduct(omega, cell->getPosition())) ;
  //Centrifugal acceleration work
  sourceConsEul.m_energ = Coord::scalarProduct(u, sourceConsEul.m_qdm);

  //Euler integration (order 1)
  m_qdm += dt*sourceConsEul.m_qdm;
  m_energ += dt*sourceConsEul.m_energ;
}

//***********************************************************************

Coord FluxEuler::getQdm() const
{
  return m_qdm;
}

//***********************************************************************

double FluxEuler::getMasseMix() const
{
  return m_masse;
}

//***********************************************************************

double FluxEuler::getEnergyMix() const
{
  return m_energ;
}

//***********************************************************************

void FluxEuler::setCons(const Flux *cons, const int &numberPhases)
{
  m_masse = cons->getMasseMix();
  m_qdm = cons->getQdm();
  m_energ = cons->getEnergyMix();
}

//***********************************************************************