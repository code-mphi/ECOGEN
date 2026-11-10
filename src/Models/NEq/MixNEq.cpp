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

#include "MixNEq.h"

using namespace tinyxml2;

//***************************************************************************

MixNEq::MixNEq() : m_density(0.), m_pressure(0.), m_velocity(0), m_energy(0.), m_totalEnergy(0.) {}

//***************************************************************************

MixNEq::~MixNEq() {}

//***************************************************************************

void MixNEq::allocateAndCopyMixture(Mixture** mixture) { *mixture = new MixNEq(*this); }

//***************************************************************************

void MixNEq::copyMixture(Mixture& mixture)
{
  m_density     = mixture.getDensity();
  m_pressure    = mixture.getPressure();
  m_velocity    = mixture.getVelocity();
  m_energy      = mixture.getEnergy();
  m_totalEnergy = mixture.getTotalEnergy();
}

//***************************************************************************

double MixNEq::computeDensity(const double* alphak, const double* rhok)
{
  double rho(0.);
  for (int k = 0; k < numberPhases; k++) {
    rho += alphak[k] * rhok[k];
  }
  return rho;
}

//***************************************************************************

double MixNEq::computePressure(const double* alphak, const double* pk)
{
  double p(0.);
  for (int k = 0; k < numberPhases; k++) {
    p += alphak[k] * pk[k];
  }
  return p;
}

//***************************************************************************

double MixNEq::computeVelocity(const double* Yk, const double* uk)
{
  double u(0.);
  for (int k = 0; k < numberPhases; k++) {
    u += Yk[k] * uk[k];
  }
  return u;
}

//***************************************************************************

double MixNEq::computeInternalEnergy(const double* Yk, const double* ek)
{
  double e(0.);
  for (int k = 0; k < numberPhases; k++) {
    e += Yk[k] * ek[k];
  }
  return e;
}

//***************************************************************************

void MixNEq::computeMixtureVariables(Phase** vecPhase)
{
  //mixture density and pressure
  m_density  = 0.;
  m_pressure = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_density  += vecPhase[k]->getAlpha() * vecPhase[k]->getDensity();
    m_pressure += vecPhase[k]->getAlpha() * vecPhase[k]->getPressure();
  }
  //Mass fraction
  for (int k = 0; k < numberPhases; k++) {
    vecPhase[k]->computeMassFraction(m_density);
  }
  //Specific internal energy and speed of sounds
  m_energy   = 0.;
  m_velocity = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_energy   += vecPhase[k]->getY() * vecPhase[k]->getEnergy();
    m_velocity += vecPhase[k]->getY() * vecPhase[k]->getVelocity();
  }
  m_totalEnergy = m_energy + 0.5 * m_velocity.squaredNorm();
}

//***************************************************************************

void MixNEq::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void MixNEq::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//**************************** DATA PRINTING *********************************
//****************************************************************************

double MixNEq::returnScalar(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return m_density;
    break;
  case 2:
    return m_pressure;
    break;
  case 3:
    return m_totalEnergy;
    break;
  case 4:
    return m_velocity.getX();
    break;
  default:
    return 0.;
    break;
  }
}

//***************************************************************************

Coord MixNEq::returnVector(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return m_velocity;
    break;
  default:
    return 0;
    break;
  }
}

//***************************************************************************

std::string MixNEq::returnNameScalar(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return "Density_Mixture";
    break;
  case 2:
    return "Pressure_Mixture";
    break;
  case 3:
    return "Total_energy_Mixture";
    break;
  case 4:
    return "Velocity_x";
    break;
  default:
    return "NoName";
    break;
  }
}

//***************************************************************************

std::string MixNEq::returnNameVector(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return "Velocity_Mixture";
    break;
  default:
    return "NoName";
    break;
  }
}

//****************************************************************************
//**************************** DATA READING **********************************
//****************************************************************************

void MixNEq::setScalar(const int& numVar, const double& value)
{
  switch (numVar) {
  case 1:
    m_density = value;
    break;
  case 2:
    m_pressure = value;
    break;
  case 3:
    m_totalEnergy = value;
    break;
  default:
    Errors::errorMessage("numVar not found in MixNEq::setScalar");
    break;
  }
}

//***************************************************************************

void MixNEq::setVector(const int& numVar, const Coord& value)
{
  switch (numVar) {
  case 1:
    m_velocity = value;
    break;
  default:
    Errors::errorMessage("numVar not found in MixNEq::setVector");
    break;
  }
}

//***************************************************************************
