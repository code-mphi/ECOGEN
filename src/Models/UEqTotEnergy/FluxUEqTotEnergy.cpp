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
#include "FluxUEqTotEnergy.h"
#include "../Mixture.h"

//***********************************************************************

FluxUEqTotEnergy::FluxUEqTotEnergy(const int& numbPhases)
{
  m_alpha = new double[numbPhases];
  m_mass = new double[numbPhases];
  m_totEnerg = new double[numbPhases];
  m_alphap = new double[numbPhases];
}

//***********************************************************************

FluxUEqTotEnergy::~FluxUEqTotEnergy()
{
  delete[] m_alpha;
  delete[] m_mass;
  delete[] m_totEnerg;
  delete[] m_alphap;
}

//***********************************************************************

void FluxUEqTotEnergy::printFlux() const
{
  std::cout << m_mass << " " << m_momentum.getX() << " " << m_totEnerg << std::endl;
}

//***********************************************************************

void FluxUEqTotEnergy::addFlux(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += coefA*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alpha[k];
    m_mass[k] += coefA*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_mass[k];
    m_totEnerg[k] += coefA*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_totEnerg[k];
  }
  m_momentum += coefA*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum;
}

//***********************************************************************

void FluxUEqTotEnergy::addFlux(Flux* flux)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += static_cast<FluxUEqTotEnergy*> (flux)->m_alpha[k];
    m_mass[k] += static_cast<FluxUEqTotEnergy*> (flux)->m_mass[k];
    m_totEnerg[k] += static_cast<FluxUEqTotEnergy*> (flux)->m_totEnerg[k];
  }
  m_momentum += static_cast<FluxUEqTotEnergy*> (flux)->m_momentum;
}

//***********************************************************************

void FluxUEqTotEnergy::subtractFlux(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] -= coefA*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alpha[k];
    m_mass[k] -= coefA*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_mass[k];
    m_totEnerg[k] -= coefA*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_totEnerg[k];
  }
  m_momentum -= coefA*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum;
}

//***********************************************************************

void FluxUEqTotEnergy::multiply(double scalar)
{
    for(int k=0;k<numberPhases;k++)
    {
      m_alpha[k] *= scalar;
      m_mass[k] *= scalar;
      m_totEnerg[k] *= scalar;
    }
    m_momentum *= scalar;
}

//***********************************************************************

void FluxUEqTotEnergy::setBufferFlux(Cell& cell)
{
  static_cast<FluxUEqTotEnergy*> (fluxBuff)->buildCons(cell.getPhases(), cell.getMixture());
}

//***********************************************************************

void FluxUEqTotEnergy::buildCons(Phase** phases, Mixture* mixture)
{
	double totEnergy(0.);
  for (int k = 0; k < numberPhases; k++)
	{
    m_alpha[k] = phases[k]->getAlpha();
    m_mass[k] = phases[k]->getAlpha() * phases[k]->getDensity();
    totEnergy = phases[k]->getEos()->computeEnergy(phases[k]->getDensity(), phases[k]->getPressure());
    totEnergy += 0.5 * mixture->getVelocity().squaredNorm();
    m_totEnerg[k] = phases[k]->getAlpha() * phases[k]->getDensity() * totEnergy;
	}
	m_momentum = mixture->getDensity()*mixture->getVelocity();
}

//***********************************************************************

void FluxUEqTotEnergy::buildPrim(Phase** phases, Mixture* mixture)
{
  double pressure(0.), rhoMel(0.);
  Coord vel(0.);
  //Verification and correction if needed (alpha and mass for order 2)
  double un(0.);
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    for (int k = 0; k < numberPhases; k++) {
      if (m_alpha[k] < 0.) m_alpha[k] = 0.;
      if (m_alpha[k] > 1.) m_alpha[k] = 1.;
      un += m_alpha[k];
    }
  }
  else { // alpha = 0 is desactivated (alpha != 0)
    for (int k = 0; k < numberPhases; k++) {
      if (m_alpha[k] <= 1.e-15) m_alpha[k] = 1e-15;
      if (m_alpha[k] >= 1.-1.e-15) m_alpha[k] = 1.0 - 1e-15;
      un += m_alpha[k];
    }
  }
  for (int k = 0; k < numberPhases; k++) { m_alpha[k] /= un; }

  //Phases and mixture variables
  for (int k = 0; k < numberPhases; k++) {
      rhoMel += m_mass[k];
      phases[k]->setAlpha(m_alpha[k]);
      phases[k]->setDensity(m_mass[k] / std::max(m_alpha[k], epsilonAlphaNull));
      phases[k]->verifyAndCorrectDensityMax(m_mass[k]);
  }

  vel = m_momentum / rhoMel;
  mixture->setVelocity(vel.getX(), vel.getY(), vel.getZ());
  //Erasing small velocity variations
  if (std::fabs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (std::fabs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (std::fabs(mixture->getW()) < 1.e-8) mixture->setW(0.);

  //Compute pressures
  for (int k = 0; k < numberPhases; k++) {
    TB->Ek[k] = m_totEnerg[k] / std::max(m_mass[k], epsilonAlphaNull);
    TB->ek[k] = TB->Ek[k] - 0.5 * mixture->getVelocity().squaredNorm();
    pressure = TB->eos[k]->computePressure(phases[k]->getDensity(), TB->ek[k]);
    phases[k]->setPressure(pressure);
    phases[k]->verifyAndCorrectPhase();
  }

  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
    phases[k]->setTotalEnergy(TB->Ek[k]);
  }
  mixture->computeMixtureVariables(phases);
}

//***********************************************************************

void FluxUEqTotEnergy::setToZero()
{
  for(int k=0;k<numberPhases;k++){
    m_alpha[k] = 0.;
    m_mass[k] = 0.;
    m_totEnerg[k] = 0.;
  }
  m_momentum = 0.;
}

//***********************************************************************

void FluxUEqTotEnergy::addNonCons(double coefA, const Cell* cell, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/)
{
  Phase* phase;
  double uStar(static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_uStar);

  double pStar(0.);
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    pStar += static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alphap[k];
  }

  double pMix(cell->getMixture()->getPressure()); 
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] += -coefA * phase->getAlpha()*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_sM;

    m_totEnerg[k] += coefA * ( phase->getY() * (pStar*uStar - pMix*uStar)
     - static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alphap[k]*uStar + phase->getAlpha()*phase->getPressure()*uStar);
  }
}

//***********************************************************************

void FluxUEqTotEnergy::subtractNonCons(double coefA, const Cell* cell, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/)
{
  Phase* phase;
  double uStar(static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_uStar);

  double pStar(0.);
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    pStar += static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alphap[k];
  }

  double pMix(cell->getMixture()->getPressure());
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] -= -coefA * phase->getAlpha()*static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_sM;

    m_totEnerg[k] -= coefA * ( phase->getY() * (pStar*uStar - pMix*uStar)
      - static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alphap[k]*uStar + phase->getAlpha()*phase->getPressure()*uStar);
  }
}

//***********************************************************************

void FluxUEqTotEnergy::setCons(const Flux* cons)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] = cons->getAlpha(k);
    m_mass[k] = cons->getMass(k);
    m_totEnerg[k] = cons->getTotEnergy(k);
  }
  m_momentum = cons->getMomentum();
}

//***********************************************************************