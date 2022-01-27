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

#include "FluxEulerHomogeneous.h"

//***********************************************************************

FluxEulerHomogeneous::FluxEulerHomogeneous() {}

//***********************************************************************

FluxEulerHomogeneous::~FluxEulerHomogeneous(){}

//***********************************************************************

void FluxEulerHomogeneous::printFlux() const
{
  std::cout << m_mass << " " << m_momentum.getX() << " " << m_energ << std::endl;
}

//***********************************************************************

void FluxEulerHomogeneous::addFlux(double coefA)
{
    m_mass += coefA*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_mass;
    m_momentum   += coefA*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum;
    m_energ += coefA*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ;
}

//***********************************************************************

void FluxEulerHomogeneous::addFlux(Flux* flux)
{
  m_mass += static_cast<FluxEulerHomogeneous*> (flux)->m_mass;
  m_momentum   += static_cast<FluxEulerHomogeneous*> (flux)->m_momentum;
  m_energ += static_cast<FluxEulerHomogeneous*> (flux)->m_energ;
}

//***********************************************************************

void FluxEulerHomogeneous::subtractFlux(double coefA)
{
    m_mass -= coefA*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_mass;
    m_momentum   -= coefA*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum;
    m_energ -= coefA*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ;
}

//***********************************************************************

void FluxEulerHomogeneous::multiply(double scalar)
{
    m_mass *= scalar;
    m_momentum   *= scalar;
    m_energ *= scalar;
}

//***********************************************************************

void FluxEulerHomogeneous::setBufferFlux(Cell& cell)
{
  static_cast<FluxEulerHomogeneous*> (fluxBuff)->buildCons(cell.getPhases(), cell.getMixture());
}

//***********************************************************************

void FluxEulerHomogeneous::buildCons(Phase** phases, Mixture* mixture)
{
  double internalEnergy(0.);
  double rhok, alphak, ek;

  Phase* phase(0);
  m_mass = 0.;
  m_energ = 0.;

  for (int k = 0; k < numberPhases; k++)
  {
    phase = phases[k];
    //Mixture density calculus
    alphak = phase->getAlpha();
    rhok = phase->getDensity();
    m_mass += alphak * rhok;
    //Mixture specific internal energy calculus
    ek = phase->getEos()->computeEnergy(rhok, phase->getPressure());
    internalEnergy += alphak*rhok*ek;
  }
  m_momentum = m_mass*mixture->getVelocity();
  m_energ = internalEnergy + 0.5*m_mass*mixture->getVelocity().squaredNorm();
}

//***********************************************************************

void FluxEulerHomogeneous::buildPrim(Phase** phases, Mixture* mixture)
{
  double pressure, Tsat, internalEnergy;
  
  int liq(static_cast<ModEulerHomogeneous*> (model)->m_liq), vap(static_cast<ModEulerHomogeneous*> (model)->m_vap);

  //Simple extractions
  mixture->setVelocity(m_momentum / m_mass);
  //Erasing small velocity variations
  if (std::fabs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (std::fabs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (std::fabs(mixture->getW()) < 1.e-8) mixture->setW(0.);
  internalEnergy = m_energ / m_mass - 0.5*mixture->getVelocity().squaredNorm();
  
  //Pressure determination
  pressure = mixture->computePressure(m_mass, internalEnergy, phases, mixture, liq, vap);
  Tsat = mixture->computeTsat(phases[liq]->getEos(), phases[vap]->getEos(), pressure);
  phases[vap]->setDensity(phases[vap]->getEos()->computeDensitySaturation(pressure, Tsat, Tsat));
  phases[liq]->setDensity(phases[liq]->getEos()->computeDensitySaturation(pressure, Tsat, Tsat));
  phases[vap]->setAlpha((m_mass - phases[liq]->getDensity())/(phases[vap]->getDensity()- phases[liq]->getDensity()));
  phases[liq]->setAlpha(1. - phases[vap]->getAlpha());
  
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->setPressure(pressure);
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  //Mixture variables
  mixture->setPressure(pressure);
  mixture->setTemperature(Tsat);
  mixture->computeMixtureVariables(phases);
  //Reconstruction of total energy from total energy equation
  double totalEnergy(m_energ / m_mass);
  mixture->setTotalEnergy(totalEnergy);
}

//***********************************************************************

void FluxEulerHomogeneous::setToZero()
{
  m_mass = 0.;
  m_momentum   = 0.;
  m_energ = 0.;
}

//***********************************************************************

void FluxEulerHomogeneous::addFluxSmooth1D(double coefA, const Coord& normal, Cell* cell)
{
  Phase* phase(cell->getPhase(0));
  coefA *= normal.getX(); // Switch sign for inflow boundary
  // Contribution only on x-direction 
  if (std::fabs(normal.getY()) > 1.e-6 || std::fabs(normal.getZ()) > 1.e-6) coefA = 0.;
  m_momentum.setX(m_momentum.getX() - phase->getPressure()*coefA);
}

//***********************************************************************

void FluxEulerHomogeneous::substractFluxSmooth1D(double coefA, const Coord& normal, Cell* cell)
{
  Phase* phase(cell->getPhase(0));
  coefA *= normal.getX(); // Switch sign for inflow boundary
  // Contribution only on x-direction
  if (std::fabs(normal.getY()) > 1.e-6 || std::fabs(normal.getZ()) > 1.e-6) coefA = 0.;
  m_momentum.setX(m_momentum.getX() + phase->getPressure()*coefA);
}

//***********************************************************************

void FluxEulerHomogeneous::setCons(const Flux* cons)
{
  m_mass = cons->getMassMix();
  m_momentum = cons->getMomentum();
  m_energ = cons->getEnergyMix();
}

//***********************************************************************