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

//! \file      FluxThermalEq.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      May 04 2018

#include <cmath>
#include <algorithm>
#include "FluxThermalEq.h"
#include "../Mixture.h"

using namespace std;

FluxThermalEq *fluxBufferThermalEq;
FluxThermalEq *sourceConsThermEq;

//***********************************************************************

FluxThermalEq::FluxThermalEq(){}

//***********************************************************************

FluxThermalEq::FluxThermalEq(const int &numberPhases)
{
  m_masse = new double[numberPhases];
}

//***********************************************************************

FluxThermalEq::~FluxThermalEq()
{
  delete[] m_masse;
}

//***********************************************************************

void FluxThermalEq::printFlux() const
{
  cout << m_masse << " " << m_qdm.getX() << endl;
}

//***********************************************************************

void FluxThermalEq::addFlux(double coefA, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_masse[k] += coefA*fluxBufferThermalEq->m_masse[k];
  }
  m_qdm += coefA*fluxBufferThermalEq->m_qdm;
  m_energMixture += coefA*fluxBufferThermalEq->m_energMixture;
}

//***********************************************************************

void FluxThermalEq::subtractFlux(double coefA, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_masse[k] -= coefA*fluxBufferThermalEq->m_masse[k];
  }
  m_qdm -= coefA*fluxBufferThermalEq->m_qdm;
  m_energMixture -= coefA*fluxBufferThermalEq->m_energMixture;
}

//***********************************************************************

void FluxThermalEq::multiply(double scalar, const int &numberPhases)
{
    for(int k=0;k<numberPhases;k++)
    {
      m_masse[k] *= scalar;
    }
    m_qdm  *= scalar;
    m_energMixture *= scalar;
}

//***********************************************************************

void FluxThermalEq::setBufferFlux(Cell &cell, const int &numberPhases)
{
  fluxBufferThermalEq->buildCons(cell.getPhases(), numberPhases, cell.getMixture());
}

//***********************************************************************

void FluxThermalEq::buildCons(Phase **phases, const int &numberPhases, Mixture *mixture)
{
	double energieInterne(0.);
	Phase *phase(0);
  for (int k = 0; k < numberPhases; k++)
	{
		phase = phases[k];
    TB->ak[k] = phase->getAlpha();
    TB->rhok[k] = phase->getDensity();
		m_masse[k] = TB->ak[k] * TB->rhok[k];
	}
	m_qdm = mixture->getDensity()*mixture->getVelocity();
  m_energMixture = mixture->getDensity()*mixture->getTotalEnergy();
}

//***********************************************************************

void FluxThermalEq::buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases)
{
  double pressure(0.), temperature(0.), internalEnergy(0.), rhoMel(0.);

  //Simple extractions
  for (int k = 0; k < numberPhases; k++) {
    rhoMel += m_masse[k];
  }
  mixture->setVelocity(m_qdm.getX() / rhoMel, m_qdm.getY() / rhoMel, m_qdm.getZ() / rhoMel);
  //Erasing small velocity variations
  if (abs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (abs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (abs(mixture->getW()) < 1.e-8) mixture->setW(0.);
  internalEnergy = m_energMixture / rhoMel - 0.5*(mixture->getU()*mixture->getU() + mixture->getV()*mixture->getV() + mixture->getW()*mixture->getW()) ;
  
  //Pressure and temperature determination
  pressure = mixture->computePressure(m_masse, internalEnergy, phases, numberPhases);
  temperature = mixture->computeTemperature(m_masse, pressure, phases, numberPhases);

  //Phasic variables
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->setPressure(pressure);
    phases[k]->setDensity(phases[k]->getEos()->computeDensity(pressure, temperature));
    phases[k]->setAlpha(m_masse[k] / phases[k]->getDensity());
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }

  mixture->computeMixtureVariables(phases, numberPhases);
  //Reconstruction of total energy from total energy equation
  double totalEnergyMixture(0.);
  totalEnergyMixture = m_energMixture / rhoMel;
  mixture->setTotalEnergy(totalEnergyMixture);
}

//***********************************************************************

void FluxThermalEq::setToZero(const int &numberPhases)
{
  for(int k=0;k<numberPhases;k++){
    m_masse[k] = 0.;
  }
  m_qdm = 0.;
  m_energMixture = 0.;
}

//***********************************************************************

void FluxThermalEq::setToZeroBufferFlux(const int &numberPhases)
{
  for (int k = 0; k<numberPhases; k++) {
    fluxBufferThermalEq->m_masse[k] = 0.;
  }
  fluxBufferThermalEq->m_qdm = 0.;
  fluxBufferThermalEq->m_energMixture = 0.;
}

//***********************************************************************

void FluxThermalEq::integrateSourceTermsHeating(Cell *cell, const double &dt, const int &numberPhases, const double &q)
{
  sourceConsThermEq->setToZero(numberPhases);
  sourceConsThermEq->m_energMixture = q;

  m_energMixture += dt*sourceConsThermEq->m_energMixture;
}

//***********************************************************************

double FluxThermalEq::getMasse(const int &numPhase) const
{
  return m_masse[numPhase];
}

//***********************************************************************

Coord FluxThermalEq::getQdm() const
{
  return m_qdm;
}

//***********************************************************************

double FluxThermalEq::getEnergyMix() const
{
	return m_energMixture;
}

//***********************************************************************

void FluxThermalEq::setCons(const Flux *cons, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_masse[k] = cons->getMasse(k);
  }
  m_qdm = cons->getQdm();
  m_energMixture = cons->getEnergyMix();
}

//***********************************************************************