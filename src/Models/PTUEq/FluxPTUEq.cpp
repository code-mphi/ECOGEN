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

#include <cmath>
#include <algorithm>
#include "FluxPTUEq.h"
#include "../Mixture.h"

//***********************************************************************

FluxPTUEq::FluxPTUEq(const int& numberPhases)
{
  m_masse = new double[numberPhases];
}

//***********************************************************************

FluxPTUEq::~FluxPTUEq()
{
  delete[] m_masse;
}

//***********************************************************************

void FluxPTUEq::printFlux() const
{
  std::cout << m_masse << " " << m_qdm.getX() << std::endl;
}

//***********************************************************************

void FluxPTUEq::addFlux(double coefA, const int& numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_masse[k] += coefA*static_cast<FluxPTUEq*> (fluxBuff)->m_masse[k];
  }
  m_qdm += coefA*static_cast<FluxPTUEq*> (fluxBuff)->m_qdm;
  m_energMixture += coefA*static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture;
}

//***********************************************************************

void FluxPTUEq::addFlux(Flux* flux, const int& numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_masse[k] += static_cast<FluxPTUEq*> (flux)->m_masse[k];
  }
  m_qdm += static_cast<FluxPTUEq*> (flux)->m_qdm;
  m_energMixture += static_cast<FluxPTUEq*> (flux)->m_energMixture;
}

//***********************************************************************

void FluxPTUEq::subtractFlux(double coefA, const int& numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_masse[k] -= coefA*static_cast<FluxPTUEq*> (fluxBuff)->m_masse[k];
  }
  m_qdm -= coefA*static_cast<FluxPTUEq*> (fluxBuff)->m_qdm;
  m_energMixture -= coefA*static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture;
}

//***********************************************************************

void FluxPTUEq::multiply(double scalar, const int& numberPhases)
{
    for(int k=0;k<numberPhases;k++)
    {
      m_masse[k] *= scalar;
    }
    m_qdm *= scalar;
    m_energMixture *= scalar;
}

//***********************************************************************

void FluxPTUEq::setBufferFlux(Cell& cell, const int& numberPhases)
{
  static_cast<FluxPTUEq*> (fluxBuff)->buildCons(cell.getPhases(), numberPhases, cell.getMixture());
}

//***********************************************************************

void FluxPTUEq::buildCons(Phase** phases, const int& numberPhases, Mixture* mixture)
{
	Phase* phase(0);
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

void FluxPTUEq::buildPrim(Phase** phases, Mixture* mixture, const int& numberPhases)
{
  double pressure(0.), temperature(0.), internalEnergy(0.), rhoMel(0.);

  //Simple extractions
  for (int k = 0; k < numberPhases; k++) {
    rhoMel += m_masse[k];
  }
  mixture->setVelocity(m_qdm.getX() / rhoMel, m_qdm.getY() / rhoMel, m_qdm.getZ() / rhoMel);
  //Erasing small velocity variations
  if (std::fabs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (std::fabs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (std::fabs(mixture->getW()) < 1.e-8) mixture->setW(0.);
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

void FluxPTUEq::setToZero(const int& numberPhases)
{
  for(int k=0;k<numberPhases;k++){
    m_masse[k] = 0.;
  }
  m_qdm = 0.;
  m_energMixture = 0.;
}

//***********************************************************************

void FluxPTUEq::setToZeroBufferFlux(const int& numberPhases)
{
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxPTUEq*> (fluxBuff)->m_masse[k] = 0.;
  }
  static_cast<FluxPTUEq*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture = 0.;
}

//***********************************************************************

void FluxPTUEq::prepSourceTermsHeating(const int& numberPhases, const double& q)
{
  m_energMixture = q;

  //Null source components 
  m_qdm = 0.;
  for (int k = 0; k < numberPhases; k++) {
	  m_masse[k] = 0.;
  }
}

//***********************************************************************

void FluxPTUEq::setCons(const Flux* cons, const int& numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_masse[k] = cons->getMasse(k);
  }
  m_qdm = cons->getQdm();
  m_energMixture = cons->getEnergyMix();
}

//***********************************************************************