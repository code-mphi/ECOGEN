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

#include "ModNEqHLL.h"

const std::string ModNEqHLL::NAME = "NONEQ";

//***********************************************************************

ModNEqHLL::ModNEqHLL(const int& numbTransports, const int& numbPhases) : ModNEq(NAME, numbTransports)
{
  fluxBuff    = new FluxNEqHLL(numbPhases);
  fluxBuffMRF = new FluxNEqHLL(numbPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxNEqHLL(numbPhases));
  }
}

//***********************************************************************

ModNEqHLL::~ModNEqHLL() {}

//***********************************************************************

void ModNEqHLL::allocateCons(Flux** cons) { *cons = new FluxNEqHLL(numberPhases); }

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModNEqHLL::solver(
  Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/) const
{
  double pL(0.), uL(0.), vL(0.), wL(0.), rhoL(0.), alphaL(0.), EL(0.);
  double pR(0.), uR(0.), vR(0.), wR(0.), rhoR(0.), alphaR(0.), ER(0.);

  std::vector<double> sL(2);
  std::vector<double> sR(2);

  // Davis estimation for sL and sR
  sL[0] = 1e9;
  sR[0] = -1e9;
  sL[1] = 1e9;
  sR[1] = -1e9;

  for (int k = 0; k < numberPhases; k++) {
    sL[k] = std::min(sL[k], cellLeft.getPhase(k)->getVelocity().getX() - cellLeft.getPhase(k)->getSoundSpeed());
    sL[k] = std::min(sL[k], cellRight.getPhase(k)->getVelocity().getX() - cellRight.getPhase(k)->getSoundSpeed());
    sR[k] = std::max(sR[k], cellLeft.getPhase(k)->getVelocity().getX() + cellLeft.getPhase(k)->getSoundSpeed());
    sR[k] = std::max(sR[k], cellRight.getPhase(k)->getVelocity().getX() + cellRight.getPhase(k)->getSoundSpeed());
  }

  //For same sL and sR for all phases
  for (int k = 0; k < numberPhases; k++) {
    sL[k] = std::min(sL[0], sL[1]);
    sR[k] = std::max(sR[0], sR[1]);
  }

  // Uncomments the 4 line beneath and get Rusanov solver
  // for (int k = 0; k < numberPhases; k++) {
  //   sR[k] = std::max(std::abs(sL[0]), std::abs(sR[0]));
  //   sL[k] = -sR[k];
  // }

  //Sampling the solution
  for (int k = 0; k < numberPhases; k++) {
    alphaL = cellLeft.getPhase(k)->getAlpha();
    pL     = cellLeft.getPhase(k)->getPressure();
    uL     = cellLeft.getPhase(k)->getVelocity().getX();
    vL     = cellLeft.getPhase(k)->getVelocity().getY();
    wL     = cellLeft.getPhase(k)->getVelocity().getZ();
    rhoL   = cellLeft.getPhase(k)->getDensity();
    EL     = cellLeft.getPhase(k)->getTotalEnergy();

    alphaR = cellRight.getPhase(k)->getAlpha();
    pR     = cellRight.getPhase(k)->getPressure();
    uR     = cellRight.getPhase(k)->getVelocity().getX();
    vR     = cellRight.getPhase(k)->getVelocity().getY();
    wR     = cellRight.getPhase(k)->getVelocity().getZ();
    rhoR   = cellRight.getPhase(k)->getDensity();
    ER     = cellRight.getPhase(k)->getTotalEnergy();

    static_cast<FluxNEqHLL*>(fluxBuff)->m_alphaL[k] = alphaL; // HLL non-conservative terms
    static_cast<FluxNEqHLL*>(fluxBuff)->m_alphaR[k] = alphaR; // HLL non-conservative terms
    static_cast<FluxNEqHLL*>(fluxBuff)->m_sL[k]     = sL[k];  // HLL non-conservative terms
    static_cast<FluxNEqHLL*>(fluxBuff)->m_sR[k]     = sR[k];  // HLL non-conservative terms

    if (std::fabs(sL[k]) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL[k]));
    if (std::fabs(sR[k]) > 1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR[k]));

    //Sampling
    if (sL[k] > 0.) {
      static_cast<FluxNEqHLL*>(fluxBuff)->m_alpha[k] = 0.;
      static_cast<FluxNEqHLL*>(fluxBuff)->m_mass[k]  = alphaL * rhoL * uL;
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setX(alphaL * (rhoL * uL * uL + pL));
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setY(alphaL * rhoL * uL * vL);
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setZ(alphaL * rhoL * uL * wL);
      static_cast<FluxNEqHLL*>(fluxBuff)->m_energ[k] = uL * alphaL * (rhoL * EL + pL);
    }
    else if (sR[k] < 0.) {
      static_cast<FluxNEqHLL*>(fluxBuff)->m_alpha[k] = 0.;
      static_cast<FluxNEqHLL*>(fluxBuff)->m_mass[k]  = alphaR * rhoR * uR;
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setX(alphaR * (rhoR * uR * uR + pR));
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setY(alphaR * rhoR * uR * vR);
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setZ(alphaR * rhoR * uR * wR);
      static_cast<FluxNEqHLL*>(fluxBuff)->m_energ[k] = uR * alphaR * (rhoR * ER + pR);
    }
    else {
      static_cast<FluxNEqHLL*>(fluxBuff)->m_alpha[k] = 0.;
      static_cast<FluxNEqHLL*>(fluxBuff)->m_mass[k]  = (alphaL * rhoL * uL * sR[k] - alphaR * rhoR * uR * sL[k] +
                                                       sL[k] * sR[k] * (alphaR * rhoR - alphaL * rhoL)) /
                                                      (sR[k] - sL[k]);
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setX((alphaL * (rhoL * uL * uL + pL) * sR[k] - alphaR * (rhoR * uR * uR + pR) * sL[k] +
                                                              sL[k] * sR[k] * (alphaR * rhoR * uR - alphaL * rhoL * uL)) /
                                                             (sR[k] - sL[k]));
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setY(
        (alphaL * (rhoL * uL * vL) * sR[k] - alphaR * (rhoR * uR * vR) * sL[k] + sL[k] * sR[k] * (alphaR * rhoR * vR - alphaL * rhoL * vL)) /
        (sR[k] - sL[k]));
      static_cast<FluxNEqHLL*>(fluxBuff)->m_momentum[k].setZ(
        (alphaL * (rhoL * uL * wL) * sR[k] - alphaR * (rhoR * uR * wR) * sL[k] + sL[k] * sR[k] * (alphaR * rhoR * wR - alphaL * rhoL * wL)) /
        (sR[k] - sL[k]));
      static_cast<FluxNEqHLL*>(fluxBuff)->m_energ[k] = (alphaL * (rhoL * EL + pL) * uL * sR[k] - alphaR * (rhoR * ER + pR) * uR * sL[k] +
                                                        sL[k] * sR[k] * (alphaR * rhoR * ER - alphaL * rhoL * EL)) /
                                                       (sR[k] - sL[k]);
    }
  }
}

//****************************************************************************
