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

#include "FluxUEq.h"
#include "../Mixture.h"

//***********************************************************************

FluxUEq::FluxUEq(const int& numbPhases)
{
  m_alpha = new double[numbPhases];
  m_mass = new double[numbPhases];
  m_energ = new double[numbPhases];
}

//***********************************************************************

FluxUEq::~FluxUEq()
{
  delete[] m_alpha;
  delete[] m_mass;
  delete[] m_energ;
}

//***********************************************************************

void FluxUEq::printFlux() const
{
  std::cout << m_mass << " " << m_momentum.getX() << " " << m_energ << std::endl;
}

//***********************************************************************

void FluxUEq::addFlux(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += coefA*static_cast<FluxUEq*> (fluxBuff)->m_alpha[k];
    m_mass[k] += coefA*static_cast<FluxUEq*> (fluxBuff)->m_mass[k];
    m_energ[k] += coefA*static_cast<FluxUEq*> (fluxBuff)->m_energ[k];
  }
  m_momentum += coefA*static_cast<FluxUEq*> (fluxBuff)->m_momentum;
  m_energMixture += coefA*static_cast<FluxUEq*> (fluxBuff)->m_energMixture;
}

//***********************************************************************

void FluxUEq::addFlux(Flux* flux)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += static_cast<FluxUEq*> (flux)->m_alpha[k];
    m_mass[k] += static_cast<FluxUEq*> (flux)->m_mass[k];
    m_energ[k] += static_cast<FluxUEq*> (flux)->m_energ[k];
  }
  m_momentum += static_cast<FluxUEq*> (flux)->m_momentum;
  m_energMixture += static_cast<FluxUEq*> (flux)->m_energMixture;
}

//***********************************************************************

void FluxUEq::subtractFlux(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_alpha[k];
    m_mass[k] -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_mass[k];
    m_energ[k] -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_energ[k];
  }
  m_momentum -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_momentum;
  m_energMixture -= coefA*static_cast<FluxUEq*> (fluxBuff)->m_energMixture;
}

//***********************************************************************

void FluxUEq::addFluxRotatingRegion(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_alpha[k];
    m_mass[k] += coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_mass[k];
    m_energ[k] += coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_energ[k];
  }
  m_momentum += coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum;
  m_energMixture += coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_energMixture;
}

//***********************************************************************

void FluxUEq::subtractFluxRotatingRegion(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] -= coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_alpha[k];
    m_mass[k] -= coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_mass[k];
    m_energ[k] -= coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_energ[k];
  }
  m_momentum -= coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum;
  m_energMixture -= coefA*static_cast<FluxUEq*> (fluxBuffMRF)->m_energMixture;
}

//***********************************************************************

void FluxUEq::multiply(double scalar)
{
    for(int k=0;k<numberPhases;k++)
    {
      m_alpha[k] *= scalar;
      m_mass[k] *= scalar;
      m_energ[k] *= scalar;
    }
    m_momentum *= scalar;
    m_energMixture *= scalar;
}

//***********************************************************************

void FluxUEq::setBufferFlux(Cell& cell)
{
  static_cast<FluxUEq*> (fluxBuff)->buildCons(cell.getPhases(), cell.getMixture());
}

//***********************************************************************

void FluxUEq::buildCons(Phase** phases, Mixture* mixture)
{
  for (int k = 0; k < numberPhases; k++)
	{
    m_alpha[k] = phases[k]->getAlpha();
    m_mass[k] = phases[k]->getAlpha() * phases[k]->getDensity();
    m_energ[k] = phases[k]->getAlpha() * phases[k]->getDensity() * phases[k]->getEos()->computeEnergy(phases[k]->getDensity(), phases[k]->getPressure());
	}
	m_momentum = mixture->getDensity() * mixture->getVelocity();
  m_energMixture = mixture->getDensity() * mixture->getTotalEnergy();
}

//***********************************************************************

void FluxUEq::buildPrim(Phase** phases, Mixture* mixture)
{
  double pressure(0.), rhoMel(0.);
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
      //Calcul Pressure
      TB->ek[k] = m_energ[k] / std::max(m_mass[k], epsilonAlphaNull);
      pressure = TB->eos[k]->computePressure(phases[k]->getDensity(), TB->ek[k]);
      phases[k]->setPressure(pressure);
      phases[k]->verifyAndCorrectPhase(); //Conservative issue arises here with UEq when pressures are corrected but are not reset by total energy afterwards (unlike PUEq)
  }
  mixture->setVelocity(m_momentum.getX() / rhoMel, m_momentum.getY() / rhoMel, m_momentum.getZ() / rhoMel);
  //Erasing small velocity variations
  if (std::fabs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (std::fabs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (std::fabs(mixture->getW()) < 1.e-8) mixture->setW(0.);
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  mixture->computeMixtureVariables(phases);
  //Reconstruction of total energy from total energy equation
  double totalEnergyMixture(0.);
  totalEnergyMixture = m_energMixture / rhoMel;
  mixture->setTotalEnergy(totalEnergyMixture);
}

//***********************************************************************

void FluxUEq::setToZero()
{
  for(int k=0;k<numberPhases;k++){
    m_alpha[k] = 0.;
    m_mass[k] = 0.;
    m_energ[k] = 0.;
  }
  m_momentum = 0.;
  m_energMixture = 0.;
}

//***********************************************************************

void FluxUEq::addNonCons(double coefA, const Cell* cell, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/)
{
  Phase* phase;
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] += -coefA*phase->getAlpha()*static_cast<FluxUEq*> (fluxBuff)->m_sM;
    //m_energ[k] += coefA*phase->getAlpha()*phase->getPressure()*static_cast<FluxUEq*> (fluxBuff)->m_uStar; //Better results when not computed
  }
}

//***********************************************************************

void FluxUEq::subtractNonCons(double coefA, const Cell* cell, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/)
{
  Phase* phase;
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] -= -coefA*phase->getAlpha()*static_cast<FluxUEq*> (fluxBuff)->m_sM;
    //m_energ[k] -= coefA*phase->getAlpha()*phase->getPressure()*static_cast<FluxUEq*> (fluxBuff)->m_uStar; //Better results when not computed
  }
}

//***********************************************************************

void FluxUEq::schemeCorrection(Cell& cell) const
{
  //Un is in fluxBuff and fluxes are in m_cons (this object)
  double rhoN(0.), rhoNp1(0.), sumDeltaAlphaRhoInternalEnergy_k(0.);
  double DeltaRhoU2(0.), DeltaRhoInternalEnergy(0.), errorOnEnergy(0.);
  Coord momentumNp1;

  rhoN = cell.getMixture()->getDensity();
  for (int k = 0; k < numberPhases; k++) {
    rhoNp1 += m_mass[k] + static_cast<FluxUEq*> (fluxBuff)->getMass(k);  //Sum of (alpha_k * rho_k)^(n+1)
    sumDeltaAlphaRhoInternalEnergy_k += m_energ[k];
  }

  momentumNp1 = m_momentum + static_cast<FluxUEq*> (fluxBuff)->getMomentum();
  DeltaRhoU2 = momentumNp1.squaredNorm() / rhoNp1 - static_cast<FluxUEq*> (fluxBuff)->getMomentum().squaredNorm() / rhoN;
  DeltaRhoInternalEnergy = m_energMixture - 0.5*DeltaRhoU2; //Note that here the delta from the surface-tension energy is considered as negligible (0.)

  errorOnEnergy = DeltaRhoInternalEnergy - sumDeltaAlphaRhoInternalEnergy_k;
  for (int k = 0; k < numberPhases; k++) {
    m_energ[k] += cell.getPhase(k)->getAlpha() * errorOnEnergy;
  }
}

//***********************************************************************

void FluxUEq::addFluxSmooth1D(double coefA, const Coord& normal, Cell* cell)
{
  Mixture* mixture(cell->getMixture());

  coefA *= normal.getX(); // Switch sign for inflow boundary
  // Contribution only on x-direction 
  if (std::fabs(normal.getY()) > 1.e-6 || std::fabs(normal.getZ()) > 1.e-6) coefA = 0.;

  m_momentum.setX(m_momentum.getX() - mixture->getPressure()*coefA);
}

//***********************************************************************

void FluxUEq::substractFluxSmooth1D(double coefA, const Coord& normal, Cell* cell)
{
  Mixture* mixture(cell->getMixture());

  coefA *= normal.getX(); // Switch sign for inflow boundary
  // Contribution only on x-direction 
  if (std::fabs(normal.getY()) > 1.e-6 || std::fabs(normal.getZ()) > 1.e-6) coefA = 0.;

  m_momentum.setX(m_momentum.getX() + mixture->getPressure()*coefA);
}

//***********************************************************************

void FluxUEq::addSymmetricTerms(Phase** phases, Mixture* mixture, const double& r, const double& v)
{
  for (int k = 0; k<numberPhases; k++)
  {
    //Note: there is no cylindrical or spherical terms in a transport equation
    m_mass[k] += -phases[k]->getAlpha() * phases[k]->getDensity() * v / r;
    m_energ[k] += -(phases[k]->getAlpha() * phases[k]->getDensity() * phases[k]->getEnergy() + phases[k]->getAlpha() * phases[k]->getPressure()) * v / r;
  }
  m_momentum += -v / r * mixture->getDensity() * mixture->getVelocity();
  m_energMixture += -(mixture->getDensity() * mixture->getTotalEnergy() + mixture->getPressure()) * v / r;
}

//***********************************************************************

void FluxUEq::prepSourceTermsGravity(const Coord& g)
{
  //Mass and velocity extraction
  double rho(0.); 
  Coord u(0.);
  
  for (int k = 0; k < numberPhases; k++){
    rho += m_mass[k];
  }
  u = m_momentum/rho;

  //Gravity force and work
  m_momentum = rho * g;
  m_energMixture = rho * Coord::scalarProduct(g, u);

  //Null source components
  for (int k = 0; k < numberPhases; k++) {
	  m_alpha[k] = 0.;
	  m_mass[k] = 0.;
	  m_energ[k] = 0.;
  }
}

//***********************************************************************

void FluxUEq::prepSourceTermsHeating(const double& q)
{
  //Version 1 on 2p model, then relax
  m_energMixture = q;
  double sumMasses(0.);
  for (int k = 0; k < numberPhases; k++) {
    sumMasses += m_mass[k];
  }
  for (int k = 0; k < numberPhases; k++) {
    m_energ[k] = m_mass[k] / sumMasses * q;
  }

  //Null source components
  m_momentum = 0.;
  for (int k = 0; k < numberPhases; k++) {
	  m_alpha[k] = 0.;
	  m_mass[k] = 0.;
  }

  //FP//Dev To generalize
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
  //    f -= m_mass[k] * vk[k];
  //    df -= m_mass[k] * dvkdp;
  //  }
  //} while (std::fabs(f)>1e-12);

  ////Cell update
  //Phase* phase;
  //double ak;
  //for (int k = 0; k < numberPhases; k++) {
  //  phase = cell->getPhase(k);
  //  ak = m_mass[k] * vk[k];
  //  phase->setAlpha(ak);
  //  phase->setDensity(1./vk[k]);
  //  phase->setPressure(pStar);

  //  //m_alpha[k] = m_mass[k] * vk[k];
  //  //m_energ[k] = m_mass[k] * TB->eos[k]->computeEnergy(1. / vk[k], pStar);
  //}

  //
  ////m_energMixture += dt*q;

  ////cell->getMixture()->setPressure(pStar);
  //cell->fulfillState();
}

//***********************************************************************

void FluxUEq::prepSourceTermsMRF(Cell* cell, const Coord& omega)
{
  //Mass and velocity extraction
  double rho(0.); 
  Coord u(0.);
  
  for (int k = 0; k < numberPhases; k++){ 
    rho += m_mass[k];
  }
  u = m_momentum/rho;
  
  //Coriolis acceleration
  m_momentum = -2.*rho*Coord::crossProduct(omega,u);
  //Centrifugal acceleration
  m_momentum -= rho*Coord::crossProduct(omega, Coord::crossProduct(omega, cell->getPosition()));
  //Centrifugal acceleration work
  m_energMixture = Coord::scalarProduct(u, m_momentum);

  //Null source components
  for (int k = 0; k < numberPhases; k++) {
	  m_alpha[k] = 0.;
	  m_mass[k] = 0.;
	  m_energ[k] = 0.;
  }
}

//***********************************************************************

void FluxUEq::setCons(const Flux* cons)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] = cons->getAlpha(k);
    m_mass[k] = cons->getMass(k);
    m_energ[k] = cons->getEnergy(k);
  }
  m_momentum = cons->getMomentum();
  m_energMixture = cons->getEnergyMix();
}

//***********************************************************************

void FluxUEq::addNonConsMrfFlux(Phase** phase)
{
  for(int k=0; k<numberPhases; k++){
    m_alpha[k] -= phase[k]->getAlpha() * m_uStar; 
    m_energ[k] += phase[k]->getAlpha() * phase[k]->getPressure() * m_uStar;
  }
}

//***********************************************************************