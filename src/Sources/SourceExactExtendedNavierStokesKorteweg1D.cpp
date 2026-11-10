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

#include "SourceExactExtendedNavierStokesKorteweg1D.h"

//***********************************************************************

SourceExactExtendedNavierStokesKorteweg1D::SourceExactExtendedNavierStokesKorteweg1D(int physicalEntity) : SourceExact(physicalEntity) {}

//***********************************************************************

SourceExactExtendedNavierStokesKorteweg1D::~SourceExactExtendedNavierStokesKorteweg1D() {}

//***********************************************************************

void SourceExactExtendedNavierStokesKorteweg1D::integrationExactSolution(Cell* cell, const double& dt, Prim type)
{
  double rho(0.), rhoSurEta(0.), omega(0.), b(0.), c(0.), delta(0.), mu(0.), diff(0.), newOmega(0.), newDensitySurEta(0.), kappaVar(0.);
  rho       = cell->getPhase(0, type)->getDensity();
  rhoSurEta = cell->getPhase(0, type)->getDensitySurEta();
  omega     = cell->getPhase(0, type)->getOmega();
  kappaVar  = cell->getPhase(0, type)->getEos()->kappa(rho);
  mu        = 2. * rho * std::sqrt(rho * kappaVar);
  b         = tauEK * mu / (betaEK * rho * rho * rho);
  c         = lambdaEK / (betaEK * rho * rho);
  delta     = b * b - 4. * c;
  diff      = rho * (1. / rhoSurEta - 1.);

  if (delta > 0.) {
    double r_1       = 2. * c / (-b - std::sqrt(delta));
    double r_2       = c / r_1;
    double A         = (diff * r_2 - omega) / (r_2 - r_1);
    double B         = (diff * r_1 - omega) / (r_1 - r_2);
    double exp_1     = std::exp(r_1 * dt);
    double exp_2     = std::exp(r_2 * dt);
    newDensitySurEta = rho / (rho + A * exp_1 + B * exp_2);
    newOmega         = r_1 * A * exp_1 + r_2 * B * exp_2;
  }
  else if (delta < 0.) { // business regime !
    double E         = -0.5 * b;
    double F         = 0.5 * std::sqrt(-delta);
    double A_i       = diff;
    double B_i       = (omega - E * A_i) / F;
    double exp_      = std::exp(E * dt);
    double cos_      = std::cos(F * dt);
    double sin_      = std::sin(F * dt);
    newDensitySurEta = rho / (rho + exp_ * (A_i * cos_ + B_i * sin_));
    newOmega         = E * exp_ * (A_i * cos_ + B_i * sin_) + exp_ * F * (-A_i * sin_ + B_i * cos_);
  }
  else {
    double r_0       = -0.5 * b;
    double C         = diff;
    double D         = omega - r_0 * C;
    double exp_      = std::exp(r_0 * dt);
    newDensitySurEta = rho / (rho + (C + D * dt) * exp_);
    newOmega         = (r_0 * (C + D * dt) + D) * exp_;
  }

  // Update of other primitive variables with conservation of the total energy E^n+1 = E^n
  double totalEnergy(0.), SquareNormVecP(0.), SquareNormVelocity(0.), internalEnergy(0.), pressure(0.), entropy(0.), squareSoundSpeed(0.);
  totalEnergy        = cell->getPhase(0, type)->getTotalEnergy();
  SquareNormVecP     = cell->getPhase(0, type)->getVectorP().squaredNorm();
  SquareNormVelocity = cell->getPhase(0, type)->getVelocity().squaredNorm();

  internalEnergy = totalEnergy - (0.5 * SquareNormVelocity + 0.5 * betaEK * newOmega * newOmega +
                                  0.5 * lambdaEK * (1. - 1. / newDensitySurEta) * (1. - 1. / newDensitySurEta) +
                                  0.5 * epsEK * epsEK * kappaVar / rho * SquareNormVecP);

  cell->getPhase(0, type)->setDensitySurEta(newDensitySurEta);
  cell->getPhase(0, type)->setOmega(newOmega);

  cell->getPhase(0, type)->setEnergy(internalEnergy);

  pressure = cell->getPhase(0, type)->getEos()->computePressure(rho, internalEnergy);
  cell->getPhase(0, type)->setPressure(pressure);

  entropy = cell->getPhase(0, type)->getEos()->computeEntropyFromDensityAndEnergy(rho, internalEnergy);
  cell->getPhase(0, type)->setEntropy(entropy);

  squareSoundSpeed = cell->getPhase(0, type)->getEos()->computeSquareSoundSpeed(rho, internalEnergy);
  cell->getPhase(0, type)->setSquareSoundSpeed(squareSoundSpeed);
}

//***********************************************************************
