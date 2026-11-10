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

#include "RiemannSolverHLLC.h"

//***********************************************************************

RiemannSolverHLLC::RiemannSolverHLLC() : RiemannSolver() {}

//***********************************************************************

RiemannSolverHLLC::~RiemannSolverHLLC() {}

//***********************************************************************

void RiemannSolverHLLC::solver(
  const Phase& phaseLeft, const Phase& phaseRight, const double& dxLeft, const double& dxRight, double& dtMax, double& pStar, FluxEuler& m_fluxEuler)
{
  double cL, cR, sL, sR;
  double uL, uR, vL, vR, wL, wR, pL, pR, rhoL, rhoR, EL, ER;

  uL   = phaseLeft.getU();
  vL   = phaseLeft.getV();
  wL   = phaseLeft.getW();
  pL   = phaseLeft.getPressure();
  rhoL = phaseLeft.getDensity();
  cL   = phaseLeft.getSoundSpeed();
  EL   = phaseLeft.getTotalEnergy();

  uR   = phaseRight.getU();
  vR   = phaseRight.getV();
  wR   = phaseRight.getW();
  pR   = phaseRight.getPressure();
  rhoR = phaseRight.getDensity();
  cR   = phaseRight.getSoundSpeed();
  ER   = phaseRight.getTotalEnergy();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if (model->getLowMach()) {
    model->lowMachSoundSpeed(machRefMin, uL, cL, uR, cR);
  }

  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  if (std::fabs(sR) > 1.e-3) dtMax = std::min(dtMax, machRefMin * dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoL * (sL - uL)), mR(rhoR * (sR - uR));
  double sM((pR - pL + mL * uL - mR * uR) / (mL - mR));
  pStar = mL * (sM - uL) + pL;
  if (std::fabs(sM) < 1.e-8) sM = 0.;

  // Sampling
  if (sL > 0.) {
    m_fluxEuler.m_mass = rhoL * uL;
    m_fluxEuler.m_momentum.setX(rhoL * uL * uL + pL);
    m_fluxEuler.m_momentum.setY(rhoL * vL * uL);
    m_fluxEuler.m_momentum.setZ(rhoL * wL * uL);
    m_fluxEuler.m_energ = (rhoL * EL + pL) * uL;
  }
  else if (sR < 0.) {
    m_fluxEuler.m_mass = rhoR * uR;
    m_fluxEuler.m_momentum.setX(rhoR * uR * uR + pR);
    m_fluxEuler.m_momentum.setY(rhoR * vR * uR);
    m_fluxEuler.m_momentum.setZ(rhoR * wR * uR);
    m_fluxEuler.m_energ = (rhoR * ER + pR) * uR;
  }
  else if (sM >= 0.) {
    pStar              = mL * (sM - uL) + pL;
    double rhoStar     = mL / (sL - sM);
    double Estar       = EL + (sM - uL) * (sM + pL / mL);
    m_fluxEuler.m_mass = rhoStar * sM;
    m_fluxEuler.m_momentum.setX(rhoStar * sM * sM + pStar);
    m_fluxEuler.m_momentum.setY(rhoStar * sM * vL);
    m_fluxEuler.m_momentum.setZ(rhoStar * sM * wL);
    m_fluxEuler.m_energ = (rhoStar * Estar + pStar) * sM;
  }
  else {
    pStar              = mR * (sM - uR) + pR;
    double rhoStar     = mR / (sR - sM);
    double Estar       = ER + (sM - uR) * (sM + pR / mR);
    m_fluxEuler.m_mass = rhoStar * sM;
    m_fluxEuler.m_momentum.setX(rhoStar * sM * sM + pStar);
    m_fluxEuler.m_momentum.setY(rhoStar * sM * vR);
    m_fluxEuler.m_momentum.setZ(rhoStar * sM * wR);
    m_fluxEuler.m_energ = (rhoStar * Estar + pStar) * sM;
  }

  //Contact discontinuity velocity
  m_fluxEuler.m_sM = sM;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions ****************
//****************************************************************************

void RiemannSolverHLLC::solverWall(const Phase& phaseLeft, const double& dxLeft, double& dtMax, double& pStar, FluxEuler& m_fluxEuler)
{
  double cL, sL;
  double uL, pL, rhoL;

  uL   = phaseLeft.getU();
  pL   = phaseLeft.getPressure();
  rhoL = phaseLeft.getDensity();
  cL   = phaseLeft.getSoundSpeed();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if (model->getLowMach()) {
    model->lowMachSoundSpeed(machRefMin, uL, cL);
  }

  sL = std::min(uL - cL, -uL - cL);
  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  pStar = rhoL * uL * (uL - sL) + pL;

  m_fluxEuler.m_mass = 0.;
  m_fluxEuler.m_momentum.setX(pStar);
  m_fluxEuler.m_momentum.setY(0.);
  m_fluxEuler.m_momentum.setZ(0.);
  m_fluxEuler.m_energ = 0.;
}
