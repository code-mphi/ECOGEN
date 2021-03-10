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
#include "FluxUEq.h"
#include "../Mixture.h"

//***********************************************************************

FluxUEq::FluxUEq(const int& numberPhases)
{
  m_alpha = new double[numberPhases];
  m_masse = new double[numberPhases];
  m_energ = new double[numberPhases];
}

//***********************************************************************

FluxUEq::~FluxUEq()
{
  delete[] m_alpha;
  delete[] m_masse;
  delete[] m_energ;
}

//***********************************************************************

void FluxUEq::printFlux() const
{
  std::cout << m_masse << " " << m_qdm.getX() << " " << m_energ << std::endl;
}

//***********************************************************************

void FluxUEq::addFlux(double coefA, const int& numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += coefA*static_cast<FluxUEq*> (fluxBuff)->m_alpha[k];
    m_masse[k] += coefA*static_cast<FluxUEq*> (fluxBuff)->m_masse[k];
    m_energ[k] += coefA*static_cast<FluxUEq*> (fluxBuff)->m_energ[k];
  }
  m_qdm += coefA*static_cast<FluxUEq*> (fluxBuff)->m_qdm;
  m_energMixture += coefA*static_cast<FluxUEq*> (fluxBuff)->m_energMixture;
}

//***********************************************************************

void FluxUEq::addFlux(Flux* flux, const int& numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += static_cast<FluxUEq*> (flux)->m_alpha[k];
    m_masse[k] += static_cast<FluxUEq*> (flux)->m_masse[k];
    m_energ[k] += static_cast<FluxUEq*> (flux)->m_energ[k];
  }
  m_qdm += static_cast<FluxUEq*> (flux)->m_qdm;
  m_energMixture += static_cast<FluxUEq*> (flux)->m_energMixture;
}

//***********************************************************************

void FluxUEq::subtractFlux(double coefA, const int& numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_alpha[k];
    m_masse[k] -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_masse[k];
    m_energ[k] -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_energ[k];
  }
  m_qdm -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_qdm;
  m_energMixture -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_energMixture;
}

//***********************************************************************

void FluxUEq::multiply(double scalar, const int& numberPhases)
{
    for(int k=0;k<numberPhases;k++)
    {
      m_alpha[k] *= scalar;
      m_masse[k] *= scalar;
      m_energ[k] *= scalar;
    }
    m_qdm *= scalar;
    m_energMixture *= scalar;
}

//***********************************************************************

void FluxUEq::setBufferFlux(Cell& cell, const int& numberPhases)
{
  static_cast<FluxUEq*> (fluxBuff)->buildCons(cell.getPhases(), numberPhases, cell.getMixture());
}

//***********************************************************************

void FluxUEq::buildCons(Phase** phases, const int& numberPhases, Mixture* mixture)
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

void FluxUEq::buildPrim(Phase** phases, Mixture* mixture, const int& numberPhases)
{
  double pressure(0.), energieInterne(0.), rhoMel(0.);
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
      rhoMel += m_masse[k];
      phases[k]->setAlpha(m_alpha[k]);
      phases[k]->setDensity(m_masse[k] / std::max(m_alpha[k], epsilonAlphaNull));
      //Calcul Pressure
      energieInterne = m_energ[k] / std::max(m_masse[k], epsilonAlphaNull);
      pressure = TB->eos[k]->computePressure(phases[k]->getDensity(), energieInterne);
      phases[k]->setPressure(pressure);
      phases[k]->verifyAndCorrectPhase(); //KS//Conservative issue here with UEq because pressures may be corrected and are not reset by total energy afterwards unlike PUEq
  }
  mixture->setVelocity(m_qdm.getX() / rhoMel, m_qdm.getY() / rhoMel, m_qdm.getZ() / rhoMel);
  //Erasing small velocity variations
  if (std::fabs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (std::fabs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (std::fabs(mixture->getW()) < 1.e-8) mixture->setW(0.);
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

void FluxUEq::setToZero(const int& numberPhases)
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

void FluxUEq::setToZeroBufferFlux(const int& numberPhases)
{
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;
}

//***********************************************************************

void FluxUEq::addNonCons(double coefA, const Cell* cell, const int& numberPhases)
{
  Phase* phase;
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] += -coefA*phase->getAlpha()*static_cast<FluxUEq*> (fluxBuff)->m_sM;
    m_energ[k] += coefA*phase->getAlpha()*phase->getPressure()*static_cast<FluxUEq*> (fluxBuff)->m_sM;
  }
}

//***********************************************************************

void FluxUEq::subtractNonCons(double coefA, const Cell* cell, const int& numberPhases)
{
  Phase* phase;
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] -= -coefA*phase->getAlpha()*static_cast<FluxUEq*> (fluxBuff)->m_sM;
    m_energ[k] -= coefA*phase->getAlpha()*phase->getPressure()*static_cast<FluxUEq*> (fluxBuff)->m_sM;
  }
}

//***********************************************************************

void FluxUEq::schemeCorrection(const int& numberPhases) const
{
  double rhoN(0.), rhoNp1(0.), sommeDeltaek(0.);
  for (int k = 0; k < numberPhases; k++) {
    rhoN += static_cast<FluxUEq*> (fluxBuff)->getMasse(k);
    rhoNp1 += m_masse[k];
    TB->Deltaek[k] = m_energ[k] - static_cast<FluxUEq*> (fluxBuff)->getEnergy(k);
    sommeDeltaek += TB->Deltaek[k];
  }

  double DeltaE(0.), DeltarhoU2(0.), Deltae(0.);
  DeltaE = m_energMixture - static_cast<FluxUEq*> (fluxBuff)->getEnergyMix();
  DeltarhoU2 = m_qdm.squaredNorm() / rhoNp1 - static_cast<FluxUEq*> (fluxBuff)->getQdm().squaredNorm() / rhoN;
  Deltae = DeltaE - 0.5*DeltarhoU2; //Note that here the delta from the surface-tension energy is considered as negligible (0.)

  if (std::fabs(sommeDeltaek) > 1.e-8 * std::fabs(m_energMixture)) {
    for (int k = 0; k < numberPhases; k++) {
      m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->getEnergy(k) + TB->Deltaek[k] / sommeDeltaek * Deltae;
    }
  }
}

//***********************************************************************

void FluxUEq::addSymmetricTerms(Phase** phases, Mixture* mixture, const int& numberPhases, const double& r, const double& v)
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

void FluxUEq::prepSourceTermsGravity(const int& numberPhases, const Coord& g)
{
  //Mass and velocity extraction
  double rho(0.); 
  Coord u(0.);
  
  for (int k = 0; k < numberPhases; k++){
    rho += m_masse[k];
  }
  u = m_qdm/rho;

  //Gravity force and work
  m_qdm = rho * g;
  m_energMixture = rho * Coord::scalarProduct(g, u);

  //Null source components
  for (int k = 0; k < numberPhases; k++) {
	  m_alpha[k] = 0.;
	  m_masse[k] = 0.;
	  m_energ[k] = 0.;
  }
}

//***********************************************************************

void FluxUEq::prepSourceTermsHeating(const int& numberPhases, const double& q)
{
  //Version 1 on 2p model, then relax
  m_energMixture = q;
  double sumMasses(0.);
  for (int k = 0; k < numberPhases; k++) {
    sumMasses += m_masse[k];
  }
  for (int k = 0; k < numberPhases; k++) {
    m_energ[k] = m_masse[k] / sumMasses * q;
  }

  //Null source components
  m_qdm = 0.;
  for (int k = 0; k < numberPhases; k++) {
	  m_alpha[k] = 0.;
	  m_masse[k] = 0.;
  }

  //JC//Dev To generalize
  ////Version 2 directly on PUEq (Kapila) model
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
  //} while (std::fabs(f)>1e-12);

  ////Cell update
  //Phase* phase;
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

void FluxUEq::prepSourceTermsMRF(Cell* cell, const int& numberPhases, const Coord& omega)
{
  //Mass and velocity extraction
  double rho(0.); 
  Coord u(0.);
  
  for (int k = 0; k < numberPhases; k++){ 
    rho += m_masse[k];
  }
  u = m_qdm/rho;
  
  //Coriolis acceleration
  m_qdm = -2.*rho*Coord::crossProduct(omega,u);
  //Centrifugal acceleration
  m_qdm -= rho*Coord::crossProduct(omega, Coord::crossProduct(omega, cell->getPosition()));
  //Centrifugal acceleration work
  m_energMixture = Coord::scalarProduct(u, m_qdm);

  //Null source components
  for (int k = 0; k < numberPhases; k++) {
	  m_alpha[k] = 0.;
	  m_masse[k] = 0.;
	  m_energ[k] = 0.;
  }
}

//***********************************************************************

void FluxUEq::setCons(const Flux* cons, const int& numberPhases)
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