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

#include "FluxEuler.h"

//***********************************************************************

FluxEuler::FluxEuler() :
  m_mass(0.), m_momentum(0.), m_energ(0.)
{}

//***********************************************************************

FluxEuler::~FluxEuler(){}

//***********************************************************************

void FluxEuler::printFlux() const
{
  std::cout << m_mass << " " << m_momentum.getX() << " " << m_energ << std::endl;
}

//***********************************************************************

void FluxEuler::addFlux(double coefA)
{
  m_mass += coefA*static_cast<FluxEuler*> (fluxBuff)->m_mass;
  m_momentum   += coefA*static_cast<FluxEuler*> (fluxBuff)->m_momentum;
  m_energ += coefA*static_cast<FluxEuler*> (fluxBuff)->m_energ;
}

//***********************************************************************

void FluxEuler::addFlux(Flux* flux)
{
  m_mass += static_cast<FluxEuler*> (flux)->m_mass;
  m_momentum   += static_cast<FluxEuler*> (flux)->m_momentum;
  m_energ += static_cast<FluxEuler*> (flux)->m_energ;
}

//***********************************************************************

void FluxEuler::subtractFlux(double coefA)
{
    m_mass -= coefA*static_cast<FluxEuler*> (fluxBuff)->m_mass;
    m_momentum   -= coefA*static_cast<FluxEuler*> (fluxBuff)->m_momentum;
    m_energ -= coefA*static_cast<FluxEuler*> (fluxBuff)->m_energ;
}

//***********************************************************************

void FluxEuler::multiply(double scalar)
{
    m_mass *= scalar;
    m_momentum   *= scalar;
    m_energ *= scalar;
}

//***********************************************************************

void FluxEuler::setBufferFlux(Cell& cell)
{
  static_cast<FluxEuler*> (fluxBuff)->buildCons(cell.getPhases(), cell.getMixture());
}

//***********************************************************************

void FluxEuler::buildCons(Phase** phases, Mixture* /*mixture*/)
{
  Phase* phase(phases[0]);
  
  m_mass = phase->getDensity();
  m_momentum = m_mass*phase->getVelocity();
  m_energ = m_mass*phase->getTotalEnergy();
}

//***********************************************************************

void FluxEuler::buildPrim(Phase** phases, Mixture* /*mixture*/)
{
  double pressure(0.), internalEnergy(0.), totalEnergy(0.), soundSpeed(0.), temperature(0.);
  Phase* phase(phases[0]);
  Eos* eos(phase->getEos());

  phase->setDensity(m_mass);
  phase->setVelocity(m_momentum.getX() / m_mass, m_momentum.getY() / m_mass, m_momentum.getZ() / m_mass);
  //Erasing small velocity variations
  if (std::fabs(phase->getU()) < 1.e-8) phase->setU(0.);
  if (std::fabs(phase->getV()) < 1.e-8) phase->setV(0.);
  if (std::fabs(phase->getW()) < 1.e-8) phase->setW(0.);

  totalEnergy = m_energ / m_mass;
  phase->setTotalEnergy(totalEnergy);
  internalEnergy = totalEnergy - 0.5*phase->getVelocity().squaredNorm();
  phase->setEnergy(internalEnergy);
  pressure = eos->computePressure(m_mass, internalEnergy);
  phase->setPressure(pressure);
  soundSpeed = eos->computeSoundSpeed(m_mass, pressure);
  phase->setSoundSpeed(soundSpeed);
  
  temperature = eos->computeTemperature(m_mass, pressure);
  phase->setTemperature(temperature);

}

//***********************************************************************

void FluxEuler::setToZero()
{
  m_mass = 0.;
  m_momentum   = 0.;
  m_energ = 0.;
}

//***********************************************************************

void FluxEuler::addFluxSmooth1D(double coefA, const Coord& normal, Cell* cell)
{
  Phase* phase(cell->getPhase(0));
  coefA *= normal.getX(); // Switch sign for inflow boundary
  // Contribution only on x-direction 
  if (std::fabs(normal.getY()) > 1.e-6 || std::fabs(normal.getZ()) > 1.e-6) coefA = 0.;
  m_momentum.setX(m_momentum.getX() - phase->getPressure()*coefA);
}

//***********************************************************************

void FluxEuler::substractFluxSmooth1D(double coefA, const Coord& normal, Cell* cell)
{
  Phase* phase(cell->getPhase(0));
  coefA *= normal.getX(); // Switch sign for inflow boundary
  // Contribution only on x-direction
  if (std::fabs(normal.getY()) > 1.e-6 || std::fabs(normal.getZ()) > 1.e-6) coefA = 0.;
  m_momentum.setX(m_momentum.getX() + phase->getPressure()*coefA);
}

//***********************************************************************

void FluxEuler::addSymmetricTerms(Phase** phases, Mixture* /*mixture*/, const double& r, const double& v)
{
	Phase* phase(phases[0]);
	m_mass -= phase->getDensity() * v / r;
	m_momentum -= phase->getDensity()* phase->getVelocity() * v / r;
	m_energ -= (phase->getDensity() * phase->getTotalEnergy() + phase->getPressure()) * v / r;
}

//***********************************************************************

void FluxEuler::prepSourceTermsHeating(const double& q)
{
  m_energ = q;

  // Null source components
  m_mass = 0.;
  m_momentum = 0.;
}

//***********************************************************************

void FluxEuler::prepSourceTermsGravity(const Coord& g)
{
  //Mass and velocity extraction
  double rho(m_mass);
  Coord u(m_momentum/rho);

  // Gravity force and work
  m_momentum = rho * g;
  m_energ = rho * Coord::scalarProduct(g, u);

  // Null source component
  m_mass = 0.;
}

//***********************************************************************

void FluxEuler::prepSourceTermsMRF(Cell* cell, const Coord& omega)
{
  //Mass and velocity extraction
  double rho(m_mass);
  Coord u(m_momentum/rho);

  //Coriolis acceleration
  m_momentum = -2.*rho*Coord::crossProduct(omega, u);
  //Centrifugal acceleration
  m_momentum -= rho*Coord::crossProduct(omega, Coord::crossProduct(omega, cell->getPosition())) ;
  //Centrifugal acceleration work
  m_energ = Coord::scalarProduct(u, m_momentum);

  // Null source components
  m_mass = 0.;
}

//***********************************************************************

void FluxEuler::setCons(const Flux* cons)
{
  m_mass = cons->getMassMix();
  m_momentum = cons->getMomentum();
  m_energ = cons->getEnergyMix();
}

//***********************************************************************