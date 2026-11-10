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

#include "FluxNEq.h"
#include "../../Mixture.h"

//***********************************************************************

FluxNEq::FluxNEq(const int& numbPhases)
{
  m_alpha    = new double[numbPhases];
  m_mass     = new double[numbPhases];
  m_momentum = new Coord[numbPhases];
  m_energ    = new double[numbPhases];
}

//***********************************************************************

FluxNEq::~FluxNEq()
{
  delete[] m_alpha;
  delete[] m_mass;
  delete[] m_momentum;
  delete[] m_energ;
}

//***********************************************************************

void FluxNEq::printFlux() const { std::cout << m_mass << " " << m_momentum << " " << m_energ << std::endl; }

//***********************************************************************

void FluxNEq::addFlux(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    += coefA * static_cast<FluxNEq*>(fluxBuff)->m_alpha[k];
    m_mass[k]     += coefA * static_cast<FluxNEq*>(fluxBuff)->m_mass[k];
    m_momentum[k] += coefA * static_cast<FluxNEq*>(fluxBuff)->m_momentum[k];
    m_energ[k]    += coefA * static_cast<FluxNEq*>(fluxBuff)->m_energ[k];
  }
}

//***********************************************************************

void FluxNEq::addFlux(Flux* flux)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    += static_cast<FluxNEq*>(flux)->m_alpha[k];
    m_mass[k]     += static_cast<FluxNEq*>(flux)->m_mass[k];
    m_momentum[k] += static_cast<FluxNEq*>(flux)->m_momentum[k];
    m_energ[k]    += static_cast<FluxNEq*>(flux)->m_energ[k];
  }
}

//***********************************************************************

void FluxNEq::subtractFlux(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    -= coefA * static_cast<FluxNEq*>(fluxBuff)->m_alpha[k];
    m_mass[k]     -= coefA * static_cast<FluxNEq*>(fluxBuff)->m_mass[k];
    m_momentum[k] -= coefA * static_cast<FluxNEq*>(fluxBuff)->m_momentum[k];
    m_energ[k]    -= coefA * static_cast<FluxNEq*>(fluxBuff)->m_energ[k];
  }
}

//***********************************************************************

void FluxNEq::addFluxRotatingRegion(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    += coefA * static_cast<FluxNEq*>(fluxBuffMRF)->m_alpha[k];
    m_mass[k]     += coefA * static_cast<FluxNEq*>(fluxBuffMRF)->m_mass[k];
    m_momentum[k] += coefA * static_cast<FluxNEq*>(fluxBuffMRF)->m_momentum[k];
    m_energ[k]    += coefA * static_cast<FluxNEq*>(fluxBuffMRF)->m_energ[k];
  }
}

//***********************************************************************

void FluxNEq::subtractFluxRotatingRegion(double coefA)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    -= coefA * static_cast<FluxNEq*>(fluxBuffMRF)->m_alpha[k];
    m_mass[k]     -= coefA * static_cast<FluxNEq*>(fluxBuffMRF)->m_mass[k];
    m_momentum[k] -= coefA * static_cast<FluxNEq*>(fluxBuffMRF)->m_momentum[k];
    m_energ[k]    -= coefA * static_cast<FluxNEq*>(fluxBuffMRF)->m_energ[k];
  }
}

//***********************************************************************

void FluxNEq::multiply(double scalar)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    *= scalar;
    m_mass[k]     *= scalar;
    m_momentum[k] *= scalar;
    m_energ[k]    *= scalar;
  }
}

//***********************************************************************

void FluxNEq::setBufferFlux(Cell& cell) { static_cast<FluxNEq*>(fluxBuff)->buildCons(cell.getPhases(), cell.getMixture()); }

//***********************************************************************

void FluxNEq::buildCons(Phase** phases, Mixture* /*mixture*/)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    = phases[k]->getAlpha();
    m_mass[k]     = phases[k]->getAlpha() * phases[k]->getDensity();
    m_momentum[k] = phases[k]->getAlpha() * phases[k]->getDensity() * phases[k]->getVelocity();
    m_energ[k]    = phases[k]->getAlpha() * phases[k]->getDensity() * phases[k]->getTotalEnergy();
  }
}

//***********************************************************************

void FluxNEq::buildPrim(Phase** phases, Mixture* mixture)
{
  //Verification and correction if needed (alpha and mass for order 2)
  double sumAlpha(0.);
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    for (int k = 0; k < numberPhases; k++) {
      if (m_alpha[k] < 0.) m_alpha[k] = 0.;
      if (m_alpha[k] > 1.) m_alpha[k] = 1.;
      sumAlpha += m_alpha[k];
    }
  }
  else { // alpha = 0 is desactivated (alpha != 0)
    for (int k = 0; k < numberPhases; k++) {
      if (m_alpha[k] <= 1.e-15) m_alpha[k] = 1.e-15;
      if (m_alpha[k] >= 1. - 1.e-15) m_alpha[k] = 1.0 - 1.e-15;
      sumAlpha += m_alpha[k];
    }
  }
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] /= sumAlpha;
  }

  //Phases and mixture variables
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->setAlpha(m_alpha[k]);
    phases[k]->setDensity(m_mass[k] / std::max(m_alpha[k], epsilonAlphaNull));
    phases[k]->verifyAndCorrectDensityMax(m_mass[k]);
    if (m_mass[k] <= epsilonAlphaNull) {
      phases[k]->setVelocity(0.);
    }
    else {
      phases[k]->setVelocity(m_momentum[k] / m_mass[k]);
    }
    //Erasing small velocity variations
    if (std::fabs(phases[k]->getU()) < 1.e-8) phases[k]->setU(0.);
    if (std::fabs(phases[k]->getV()) < 1.e-8) phases[k]->setV(0.);
    if (std::fabs(phases[k]->getW()) < 1.e-8) phases[k]->setW(0.);
    //Calcul Pressure
    TB->ek[k] = m_energ[k] / std::max(m_mass[k], epsilonAlphaNull) - 0.5 * phases[k]->getVelocity().squaredNorm();
    phases[k]->setPressure(TB->eos[k]->computePressure(phases[k]->getDensity(), TB->ek[k]));
    phases[k]->verifyAndCorrectPhase(); //Conservative issue arises here with NEq when pressures are corrected
  }
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  mixture->computeMixtureVariables(phases);
}

//***********************************************************************

void FluxNEq::setToZero()
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    = 0.;
    m_mass[k]     = 0.;
    m_momentum[k] = 0.;
    m_energ[k]    = 0.;
  }
  m_uI = 0.;
  m_pI = 0.;
}

//***********************************************************************

void FluxNEq::addSymmetricTerms(Phase** phases, Mixture* /*mixture */, const double& r, const double& v)
{
  for (int k = 0; k < numberPhases; k++) {
    //Note: there is no cylindrical or spherical terms in a transport equation
    m_mass[k]     += -phases[k]->getAlpha() * phases[k]->getDensity() * v / r;
    m_momentum[k] += -v / r * phases[k]->getAlpha() * phases[k]->getDensity() * phases[k]->getVelocity();
    m_energ[k] += -(phases[k]->getAlpha() * phases[k]->getDensity() * phases[k]->getEnergy() + phases[k]->getAlpha() * phases[k]->getPressure()) * v /
                  r;
  }
}

//***********************************************************************

void FluxNEq::setCons(const Flux* cons)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k]    = cons->getAlpha(k);
    m_mass[k]     = cons->getMass(k);
    m_momentum[k] = cons->getMomentum(k);
    m_energ[k]    = cons->getEnergy(k);
  }
}

//***********************************************************************
