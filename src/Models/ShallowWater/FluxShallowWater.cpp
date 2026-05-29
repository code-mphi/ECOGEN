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

#include "FluxShallowWater.h"

//***********************************************************************

FluxShallowWater::FluxShallowWater() : m_mass(0.), m_momentum(0.) {}

//***********************************************************************

FluxShallowWater::~FluxShallowWater() {}

//***********************************************************************

void FluxShallowWater::addFlux(double coefA)
{
  m_mass     += coefA * static_cast<FluxShallowWater*>(fluxBuff)->m_mass;
  m_momentum += coefA * static_cast<FluxShallowWater*>(fluxBuff)->m_momentum;
}

//***********************************************************************

void FluxShallowWater::addFlux(Flux* flux)
{
  m_mass     += static_cast<FluxShallowWater*>(flux)->m_mass;
  m_momentum += static_cast<FluxShallowWater*>(flux)->m_momentum;
}

//***********************************************************************

void FluxShallowWater::subtractFlux(double coefA)
{
  m_mass     -= coefA * static_cast<FluxShallowWater*>(fluxBuff)->m_mass;
  m_momentum -= coefA * static_cast<FluxShallowWater*>(fluxBuff)->m_momentum;
}

//***********************************************************************

void FluxShallowWater::multiply(double scalar)
{
  m_mass     *= scalar;
  m_momentum *= scalar;
}

//***********************************************************************

void FluxShallowWater::setBufferFlux(Cell& cell) { static_cast<FluxShallowWater*>(fluxBuff)->buildCons(cell.getPhases(), cell.getMixture()); }

//***********************************************************************

void FluxShallowWater::buildCons(Phase** phases, Mixture* /*mixture*/)
{
  Phase* phase{phases[0]};

  // TO FILL

  // Compute conservative variables from primitive ones: m_mass = h, m_momentum = (hu, hv)
  assert(0 && "Implement conservative variables computation");
  //m_mass     =
  //m_momentum =
}

//***********************************************************************

void FluxShallowWater::buildPrim(Phase** phases, Mixture* /*mixture*/)
{
  Phase* phase{phases[0]};
  Eos* eos{phase->getEos()};

  // TO FILL
  assert(0 && "Implement primitive variables computation");

  // Compute primitive variables from conservative ones and assign them to phases: h = m_mass , (u,v) =  m_momentum / h
  // phase->setHeight(?);
  // phase->setVelocity(?);

  //Erasing small velocity variations
  if (std::fabs(phase->getU()) < 1.e-8) phase->setU(0.);
  if (std::fabs(phase->getV()) < 1.e-8) phase->setV(0.);

  // TO FILL
  assert(0 && "Call methods to compute pressure and sound speed (methods of Eos class)");

  // Compute and assign pressure and sound speed
  //double pressure =
  //phase->setPressure(pressure);
  //double soundSpeed =
  //phase->setSoundSpeed(soundSpeed);
}

//***********************************************************************

void FluxShallowWater::setToZero()
{
  m_mass     = 0.;
  m_momentum = 0.;
}

//***********************************************************************
