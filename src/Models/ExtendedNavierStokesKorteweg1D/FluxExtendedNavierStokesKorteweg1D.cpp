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

#include "FluxExtendedNavierStokesKorteweg1D.h"

//***********************************************************************

FluxExtendedNavierStokesKorteweg1D::FluxExtendedNavierStokesKorteweg1D() :
  m_mass(0.), m_eqOmega(0.), m_eqDensitySurEta(0.), m_momentum(0.), m_eqVectorP(0.), m_energ(0.)
{}

//***********************************************************************

FluxExtendedNavierStokesKorteweg1D::~FluxExtendedNavierStokesKorteweg1D() {}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::printFlux() const
{
  std::cout << m_mass << " " << m_eqOmega << " " << m_eqDensitySurEta << " " << m_momentum.getX() << " " << m_eqVectorP.getX() << " " << m_energ
            << std::endl;
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::addFlux(double coefA)
{
  m_mass            += coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_mass;
  m_eqOmega         += coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqOmega;
  m_eqDensitySurEta += coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqDensitySurEta;
  m_momentum        += coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum;
  m_eqVectorP       += coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP;
  m_energ           += coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_energ;
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::addFlux(Flux* flux)
{
  m_mass            += static_cast<FluxExtendedNavierStokesKorteweg1D*>(flux)->m_mass;
  m_eqOmega         += static_cast<FluxExtendedNavierStokesKorteweg1D*>(flux)->m_eqOmega;
  m_eqDensitySurEta += static_cast<FluxExtendedNavierStokesKorteweg1D*>(flux)->m_eqDensitySurEta;
  m_momentum        += static_cast<FluxExtendedNavierStokesKorteweg1D*>(flux)->m_momentum;
  m_eqVectorP       += static_cast<FluxExtendedNavierStokesKorteweg1D*>(flux)->m_eqVectorP;
  m_energ           += static_cast<FluxExtendedNavierStokesKorteweg1D*>(flux)->m_energ;
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::subtractFlux(double coefA)
{
  m_mass            -= coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_mass;
  m_eqOmega         -= coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqOmega;
  m_eqDensitySurEta -= coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqDensitySurEta;
  m_momentum        -= coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum;
  m_eqVectorP       -= coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP;
  m_energ           -= coefA * static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_energ;
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::multiply(double scalar)
{
  m_mass            *= scalar;
  m_eqOmega         *= scalar;
  m_eqDensitySurEta *= scalar;
  m_momentum        *= scalar;
  m_eqVectorP       *= scalar;
  m_energ           *= scalar;
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::setBufferFlux(Cell& cell)
{
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->buildCons(cell.getPhases(), cell.getMixture());
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::buildCons(Phase** phases, Mixture* /*mixture*/)
{
  Phase* phase(phases[0]);

  m_mass            = phase->getDensity();
  m_eqOmega         = m_mass * phase->getOmega();
  m_eqDensitySurEta = phase->getDensitySurEta();
  m_momentum        = m_mass * phase->getVelocity();
  m_eqVectorP       = phase->getVectorP();
  m_energ           = m_mass * phase->getTotalEnergy();
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::buildPrim(Phase** phases, Mixture* /*mixture*/)
{
  double densitySurEta(0.), omega(0.), pressure(0.), internalEnergy(0.), totalEnergy(0.), entropy(0.), squareSoundSpeed(0.);
  Phase* phase(phases[0]);
  Eos* eos(phase->getEos());

  phase->setDensity(m_mass);

  omega = m_eqOmega / m_mass;
  phase->setOmega(omega);

  densitySurEta = m_eqDensitySurEta;
  phase->setDensitySurEta(densitySurEta);

  phase->setVectorP(m_eqVectorP.getX(), m_eqVectorP.getY(), m_eqVectorP.getZ());

  phase->setVelocity(m_momentum.getX() / m_mass, m_momentum.getY() / m_mass, m_momentum.getZ() / m_mass);
  //Erasing small velocity variations
  if (std::fabs(phase->getU()) < 1.e-8) phase->setU(0.);
  if (std::fabs(phase->getV()) < 1.e-8) phase->setV(0.);
  if (std::fabs(phase->getW()) < 1.e-8) phase->setW(0.);

  totalEnergy = m_energ / m_mass;
  phase->setTotalEnergy(totalEnergy);

  internalEnergy = totalEnergy - (0.5 * phase->getVelocity().squaredNorm() + 0.5 * betaEK * omega * omega +
                                  0.5 * lambdaEK * (1. - 1. / densitySurEta) * (1. - 1. / densitySurEta) +
                                  0.5 * epsEK * epsEK * eos->kappa(m_mass) / m_mass * phase->getVectorP().squaredNorm());
  phase->setEnergy(internalEnergy);

  pressure = eos->computePressure(m_mass, internalEnergy);
  phase->setPressure(pressure);

  entropy = eos->computeEntropyFromDensityAndEnergy(m_mass, internalEnergy);
  phase->setEntropy(entropy);

  squareSoundSpeed = eos->computeSquareSoundSpeed(m_mass, internalEnergy);
  phase->setSquareSoundSpeed(squareSoundSpeed);
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::setToZero()
{
  m_mass            = 0.;
  m_eqOmega         = 0.;
  m_eqDensitySurEta = 0.;
  m_momentum        = 0.;
  m_eqVectorP       = 0.;
  m_energ           = 0.;
}

//***********************************************************************

void FluxExtendedNavierStokesKorteweg1D::setCons(const Flux* cons)
{
  m_mass            = cons->getMassMix();
  m_eqOmega         = cons->getEqOmega();
  m_eqDensitySurEta = cons->getEqDensitySurEta();
  m_momentum        = cons->getMomentum();
  m_eqVectorP       = cons->getEqVectorP();
  m_energ           = cons->getEnergyMix();
}

//***********************************************************************
