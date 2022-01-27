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

#include "FluxEulerKorteweg.h"

//***********************************************************************
double alphaEK, betaEK, temperatureEK, kappaEK;   

FluxEulerKorteweg::FluxEulerKorteweg() : m_mass(0.), m_eqOmega(0.), m_eqEta(0.), m_momentum(0.), m_eqVectorP(0.) 
{}

//***********************************************************************

FluxEulerKorteweg::~FluxEulerKorteweg(){}

//***********************************************************************

void FluxEulerKorteweg::printFlux() const
{
  std::cout << m_mass << " " << m_eqOmega << " " << m_eqEta << " " << m_momentum.getX() << " " << m_eqVectorP.getX() << std::endl;
}

//***********************************************************************

void FluxEulerKorteweg::addFlux(double coefA)
{
  m_mass += coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_mass;
  m_eqOmega += coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqOmega;
  m_eqEta += coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqEta;
  m_momentum += coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum;
  m_eqVectorP += coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP;
}

//***********************************************************************

void FluxEulerKorteweg::addFlux(Flux* flux)
{
  m_mass += static_cast<FluxEulerKorteweg*> (flux)->m_mass;
  m_eqOmega += static_cast<FluxEulerKorteweg*> (flux)->m_eqOmega;
  m_eqEta += static_cast<FluxEulerKorteweg*> (flux)->m_eqEta;
  m_momentum += static_cast<FluxEulerKorteweg*> (flux)->m_momentum;
  m_eqVectorP += static_cast<FluxEulerKorteweg*> (flux)->m_eqVectorP;
}

//***********************************************************************

void FluxEulerKorteweg::subtractFlux(double coefA)
{
  m_mass -= coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_mass;
  m_eqOmega -= coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqOmega;
  m_eqEta -= coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqEta;
  m_momentum -= coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum;
  m_eqVectorP -= coefA*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP;
}

//***********************************************************************

void FluxEulerKorteweg::multiply(double scalar)
{
  m_mass *= scalar;
  m_eqOmega *= scalar;
  m_eqEta *= scalar;
  m_momentum *= scalar;
  m_eqVectorP *= scalar; 
}

//***********************************************************************

void FluxEulerKorteweg::setBufferFlux(Cell& cell)
{
  static_cast<FluxEulerKorteweg*> (fluxBuff)->buildCons(cell.getPhases(), cell.getMixture());
}

//***********************************************************************

void FluxEulerKorteweg::buildCons(Phase** phases, Mixture* /*mixture*/)
{
  Phase* phase(phases[0]);
  
  m_mass = phase->getDensity();
  m_eqOmega = m_mass*phase->getOmega();
  m_eqEta = m_mass*phase->getEta();
  m_momentum = m_mass*phase->getVelocity();
  m_eqVectorP = phase->getVectorP();
}

//***********************************************************************

void FluxEulerKorteweg::buildPrim(Phase** phases, Mixture* /*mixture*/)
{
  Phase* phase(phases[0]);

  phase->setDensity(m_mass);
  phase->setOmega(m_eqOmega / m_mass);
  phase->setEta(m_eqEta / m_mass);
  if (phase->getEos() != nullptr) phase->setPressure(phase->getEos()->computePressure(m_mass, temperatureEK));

  phase->setVelocity(m_momentum.getX() / m_mass, m_momentum.getY() / m_mass, m_momentum.getZ() / m_mass);
  //Erasing small velocity variations
  if (std::fabs(phase->getU()) < 1.e-8) phase->setU(0.);
  if (std::fabs(phase->getV()) < 1.e-8) phase->setV(0.);
  if (std::fabs(phase->getW()) < 1.e-8) phase->setW(0.);
  
  phase->setVectorP(m_eqVectorP.getX(), m_eqVectorP.getY(), m_eqVectorP.getZ());
}

//***********************************************************************

void FluxEulerKorteweg::setToZero()
{
  m_mass = 0.;
  m_eqOmega = 0.;
  m_eqEta = 0.; 
  m_momentum = 0.;
  m_eqVectorP = 0.;
}

//***********************************************************************

void FluxEulerKorteweg::setCons(const Flux* cons)
{
  m_mass = cons->getMassMix();
  m_eqOmega = cons->getEqOmega(); 
  m_eqEta = cons->getEqEta();
  m_momentum = cons->getMomentum();
  m_eqVectorP = cons->getEqVectorP();
}

//***********************************************************************