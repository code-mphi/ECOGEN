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

//! \file      FluxKapila.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 15 2018

#include <cmath>
#include <algorithm>
#include "FluxKapila.h"
#include "../Mixture.h"

using namespace std;

FluxKapila *fluxBufferKapila;
FluxKapila *sourceConsKap;

//***********************************************************************

FluxKapila::FluxKapila(){}

//***********************************************************************

FluxKapila::FluxKapila(ModKapila *model, const int &numberPhases) : m_model(model)
{
  m_alpha = new double[numberPhases];
  m_masse = new double[numberPhases];
  m_energ = new double[numberPhases];
}

//***********************************************************************

FluxKapila::~FluxKapila()
{
  delete[] m_alpha;
  delete[] m_masse;
  delete[] m_energ;
}

//***********************************************************************

void FluxKapila::printFlux() const
{
  cout << m_masse << " " << m_qdm.getX() << " " << m_energ << endl;
}

//***********************************************************************

void FluxKapila::addFlux(double coefA, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += coefA*fluxBufferKapila->m_alpha[k];
    m_masse[k] += coefA*fluxBufferKapila->m_masse[k];
    m_energ[k] += coefA*fluxBufferKapila->m_energ[k];
  }
  m_qdm += coefA*fluxBufferKapila->m_qdm;
  m_energMixture += coefA*fluxBufferKapila->m_energMixture;
}

//***********************************************************************

void FluxKapila::subtractFlux(double coefA, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] -= coefA*fluxBufferKapila->m_alpha[k];
    m_masse[k] -= coefA*fluxBufferKapila->m_masse[k];
    m_energ[k] -= coefA*fluxBufferKapila->m_energ[k];
  }
  m_qdm -= coefA*fluxBufferKapila->m_qdm;
  m_energMixture -= coefA*fluxBufferKapila->m_energMixture;
}

//***********************************************************************

void FluxKapila::multiply(double scalar, const int &numberPhases)
{
    for(int k=0;k<numberPhases;k++)
    {
      m_alpha[k] *= scalar;
      m_masse[k] *= scalar;
      m_energ[k] *= scalar;
    }
    m_qdm  *= scalar;
    m_energMixture *= scalar;
}

//***********************************************************************

void FluxKapila::setBufferFlux(Cell &cell, const int &numberPhases)
{
  fluxBufferKapila->buildCons(cell.getPhases(), numberPhases, cell.getMixture());
}

//***********************************************************************

void FluxKapila::buildCons(Phase **phases, const int &numberPhases, Mixture *mixture)
{
	double energieInterne(0.);
  for (int k = 0; k < numberPhases; k++)
	{
    m_alpha[k] = phases[k]->getAlpha();
    m_masse[k] = phases[k]->getAlpha() * phases[k]->getDensity();
    energieInterne = phases[k]->getEos()->computeEnergy(phases[k]->getDensity(), phases[k]->getPressure());
    m_energ[k] = phases[k]->getAlpha() * phases[k]->getDensity() * energieInterne;
	}
	m_qdm = mixture->getDensity()*mixture->getVelocity();
  m_energMixture = mixture->getDensity()*mixture->getTotalEnergy();
}

//***********************************************************************

void FluxKapila::buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases)
{
  double pressure(0.), energieInterne(0.), rhoMel(0.);

  //Verification and correction if needed (alpha and mass for order 2)
  double un(0.);
  if (epsilon > 1.e-20) { // alpha = 0 is activated
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
      rhoMel = rhoMel + m_masse[k];
      phases[k]->setAlpha(m_alpha[k]);
      phases[k]->setDensity(m_masse[k] / max(m_alpha[k], epsilon));
      //Calcul Pressure
      energieInterne = m_energ[k] / max(m_masse[k], epsilon);
      pressure = TB->eos[k]->computePressure(phases[k]->getDensity(), energieInterne);
      phases[k]->setPressure(pressure);
      phases[k]->verifyAndCorrectPhase();
  }
  mixture->setVelocity(m_qdm.getX() / rhoMel, m_qdm.getY() / rhoMel, m_qdm.getZ() / rhoMel);
  //Erasing small velocity variations
  if (abs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (abs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (abs(mixture->getW()) < 1.e-8) mixture->setW(0.);
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  mixture->computeMixtureVariables(phases, numberPhases);
  //Reconstruction of total energy from total energy equation
  double totalEnergyMixture(0.);
  totalEnergyMixture = m_energMixture / rhoMel;
  mixture->setTotalEnergy(totalEnergyMixture);
}

//***********************************************************************

void FluxKapila::setToZero(const int &numberPhases)
{
  for(int k=0;k<numberPhases;k++){
    m_alpha[k] = 0.;
    m_masse[k] = 0.;
    m_energ[k] = 0.;
  }
  m_qdm = 0.;
  m_energMixture = 0.;
}

//***********************************************************************

void FluxKapila::setToZeroBufferFlux(const int &numberPhases)
{
  for (int k = 0; k<numberPhases; k++) {
    fluxBufferKapila->m_alpha[k] = 0.;
    fluxBufferKapila->m_masse[k] = 0.;
    fluxBufferKapila->m_energ[k] = 0.;
  }
  fluxBufferKapila->m_qdm = 0.;
  fluxBufferKapila->m_energMixture = 0.;
}

//***********************************************************************

void FluxKapila::addNonCons(double coefA, const Cell *cell, const int &numberPhases)
{
  Phase *phase;
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] += -coefA*phase->getAlpha()*fluxBufferKapila->m_sM;
    m_energ[k] += coefA*phase->getAlpha()*phase->getPressure()*fluxBufferKapila->m_sM;
  }
}

//***********************************************************************

void FluxKapila::subtractNonCons(double coefA, const Cell *cell, const int &numberPhases)
{
  Phase *phase;
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] -= -coefA*phase->getAlpha()*fluxBufferKapila->m_sM;
    m_energ[k] -= coefA*phase->getAlpha()*phase->getPressure()*fluxBufferKapila->m_sM;
  }
}

//***********************************************************************

void FluxKapila::correctionEnergy(Cell *cell, const int &numberPhases, Prim type) const
{
  Phase *phase;

  //Usefull data extraction
  for (int k = 0; k < numberPhases; k++){
    phase = cell->getPhase(k, type);
    TB->ak[k] = phase->getAlpha();
    TB->rhok[k] = phase->getDensity();
  }
  //Mixture pressure calculus from mixture EOS
  double rhoe = cell->getMixture(type)->getDensity() * cell->getMixture(type)->getEnergy();
  double p(rhoe), denom(0.), gamPinfSurGamMoinsUn(0.), eRef(0.), unSurGamMoinsUn(0.);
  for (int k = 0; k < numberPhases; k++) {
    TB->eos[k]->sendSpecialMixtureEos(gamPinfSurGamMoinsUn, eRef, unSurGamMoinsUn);
    p -= TB->ak[k]*(gamPinfSurGamMoinsUn + TB->rhok[k] * eRef);
    denom += TB->ak[k] * unSurGamMoinsUn;
  }
  p /= denom;

  //Cell update
  for (int k = 0; k < numberPhases; k++){
    phase = cell->getPhase(k, type);
    phase->setPressure(p);
    phase->verifyAndCorrectPhase();
  }
  cell->getMixture(type)->setPressure(p);

  cell->fulfillState();
}

//***********************************************************************

void FluxKapila::addSymmetricTerms(Phase **phases, Mixture *mixture, const int &numberPhases, const double &r, const double &v)
{
  for (int k = 0; k<numberPhases; k++)
  {
    //Note: there is no cylindrical or spherical terms in a transport equation
    m_masse[k] += -phases[k]->getAlpha() * phases[k]->getDensity() * v / r;
    m_energ[k] += -(phases[k]->getAlpha() * phases[k]->getDensity() * phases[k]->getEnergy() + phases[k]->getAlpha() * phases[k]->getPressure()) * v / r;
  }
  m_qdm += -v / r * mixture->getDensity() * mixture->getVelocity();
  m_energMixture += -(mixture->getDensity() * mixture->getTotalEnergy() + mixture->getPressure()) * v / r;
}

//***********************************************************************

void FluxKapila::integrateSourceTermsGravity(Cell *cell, const double &dt, const int &numberPhases, const int &axe, const int &direction, const Coord &g)
{
  sourceConsKap->setToZero(numberPhases);
  //Mass and velocity extraction
  double rho = cell->getMixture()->getDensity();
  Coord u = cell->getMixture()->getVelocity();

  //Gravity force and work
  sourceConsKap->m_qdm = rho*g;
  sourceConsKap->m_energMixture = rho*Coord::scalarProduct(g,u);

  //Euler integration (order 1)
  m_qdm += dt*sourceConsKap->m_qdm;
  m_energMixture += dt*sourceConsKap->m_energMixture;
}

//***********************************************************************

void FluxKapila::integrateSourceTermsHeating(Cell *cell, const double &dt, const int &numberPhases, const double &q)
{
  sourceConsKap->setToZero(numberPhases);

  //Version 1 on 2p model, then relax
  sourceConsKap->m_energMixture = q;
  double sumMasses(0.);
  for (int k = 0; k < numberPhases; k++) {
    sumMasses += m_masse[k];
  }
    for (int k = 0; k < numberPhases; k++) {
    sourceConsKap->m_energ[k] = m_masse[k] / sumMasses * q;
    m_energ[k] += dt*sourceConsKap->m_energ[k];
  }
  m_energMixture += dt*sourceConsKap->m_energMixture;


  ////Version 2 directly on Kapila model
  //double p0 = cell->getMixture()->getPressure();
  //double v = 1./cell->getMixture()->getDensity();
  //vector<double> vk(numberPhases);
  //vector<double> vk0(numberPhases);
  //for (int k = 0; k < numberPhases; k++) {
  //  vk0[k] = 1./cell->getPhase(k)->getDensity();
  //}
  //double pStar(p0+1.);
  //double dvkdp(0.);

  //int iteration(0);
  //double f(0.), df(1.);
  //do {
  //  pStar -= f / df; iteration++;
  //  if (iteration > 50) {
  //    errors.push_back(Errors("not converged in integrateSourceTermsHeating", __FILE__, __LINE__));
  //    break;
  //  }
  //  //Physical pressure?
  //  for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(pStar); }
  //  //Specific volumes after heating
  //  f = 1.;
  //  df = 0.;
  //  for (int k = 0; k < numberPhases; k++) {
  //    vk[k] = TB->eos[k]->computeSpecificVolumeQ(p0, vk0[k], pStar, q*v*dt, &dvkdp);
  //    f -= m_masse[k] * vk[k];
  //    df -= m_masse[k] * dvkdp;
  //  }
  //} while (abs(f)>1e-12);

  ////Cell update
  //Phase *phase;
  //double ak;
  //for (int k = 0; k < numberPhases; k++) {
  //  phase = cell->getPhase(k);
  //  ak = m_masse[k] * vk[k];
  //  phase->setAlpha(ak);
  //  phase->setDensity(1./vk[k]);
  //  phase->setPressure(pStar);

  //  //m_alpha[k] = m_masse[k] * vk[k];
  //  //m_energ[k] = m_masse[k] * TB->eos[k]->computeEnergy(1. / vk[k], pStar);
  //}

  //
  ////m_energMixture += dt*q;

  ////cell->getMixture()->setPressure(pStar);
  //cell->fulfillState();
}

//***********************************************************************

void FluxKapila::integrateSourceTermsMRF(Cell *cell, const double &dt, const int &numberPhases, const Coord &omega)
{
  sourceConsKap->setToZero(numberPhases);
  //Mass and velocity extraction
  double rho = cell->getMixture()->getDensity();
  Coord u = cell->getMixture()->getVelocity();
  
  //Coriolis acceleration
  sourceConsKap->m_qdm = -2.*rho*Coord::crossProduct(omega,u);
  //Centrifugal acceleration
  sourceConsKap->m_qdm -= rho*Coord::crossProduct(omega, Coord::crossProduct(omega, cell->getPosition()));
  //Centrifugal acceleration work
  sourceConsKap->m_energMixture = Coord::scalarProduct(u, sourceConsKap->m_qdm);

  //Euler integration (order 1)
  m_qdm += dt*sourceConsKap->m_qdm;
  m_energMixture += dt*sourceConsKap->m_energMixture;
}

//***********************************************************************

double FluxKapila::getAlpha(const int &numPhase) const
{
  return m_alpha[numPhase];
}

//***********************************************************************

double FluxKapila::getMasse(const int &numPhase) const
{
  return m_masse[numPhase];
}

//***********************************************************************

double FluxKapila::getEnergy(const int &numPhase) const
{
  return m_energ[numPhase];
}

//***********************************************************************

Coord FluxKapila::getQdm() const
{
  return m_qdm;
}

//***********************************************************************

double FluxKapila::getEnergyMix() const
{
	return m_energMixture;
}

//***********************************************************************

void FluxKapila::setCons(const Flux *cons, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] = cons->getAlpha(k);
    m_masse[k] = cons->getMasse(k);
    m_energ[k] = cons->getEnergy(k);
  }
  m_qdm = cons->getQdm();
  m_energMixture = cons->getEnergyMix();
}

//***********************************************************************