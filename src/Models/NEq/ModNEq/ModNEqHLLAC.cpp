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

#include "ModNEqHLLAC.h"

const std::string ModNEqHLLAC::NAME = "NONEQ";

//***********************************************************************

ModNEqHLLAC::ModNEqHLLAC(const int& numbTransports, const int& numbPhases) : ModNEq(NAME, numbTransports)
{
  fluxBuff    = new FluxNEqHLLAC(numbPhases);
  fluxBuffMRF = new FluxNEqHLLAC(numbPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxNEqHLLAC(numbPhases));
  }
}

//***********************************************************************

ModNEqHLLAC::~ModNEqHLLAC() {}

//***********************************************************************

void ModNEqHLLAC::allocateCons(Flux** cons) { *cons = new FluxNEqHLLAC(numberPhases); }

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModNEqHLLAC::solver(
  Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/) const
{
  std::vector<double> sL(2), sR(2), uStarL(2), uStarR(2), rhoStarL(2), rhoStarR(2);
  double pStar_0(0.), pStarL_1(0.), pStarR_1(0.);
  double sM_0(0.), uI(0.);
  bool configA(true);
  // Additional variables for estimation of intermediates states
  double HMean(0.), uMean(0.), cMean(0.), gamma(0.);
  double a(0.), b(0.), c(0.), d(0.), e(0.), f(0.), g(0.), h(0.), RHS1(0.), RHS2(0.);
  // Phase variables
  // ---------------
  std::vector<double> alphaL(2), alphaR(2), rhoL(2), rhoR(2), pL(2), pR(2), uL(2), uR(2), HL(2), HR(2), cL(2), cR(2), EL(2), ER(2), coefL(2),
    coefR(2);
  for (int k = 0; k < numberPhases; k++) {
    alphaL[k] = cellLeft.getPhase(k)->getAlpha();
    uL[k]     = cellLeft.getPhase(k)->getVelocity().getX();
    rhoL[k]   = cellLeft.getPhase(k)->getDensity();
    pL[k]     = cellLeft.getPhase(k)->getPressure();
    HL[k]     = cellLeft.getPhase(k)->getEos()->computeTotalEnthalpy(rhoL[k], pL[k], cellLeft.getPhase(k)->getVelocity());
    cL[k]     = cellLeft.getPhase(k)->getSoundSpeed();
    EL[k]     = cellLeft.getPhase(k)->getTotalEnergy();

    alphaR[k] = cellRight.getPhase(k)->getAlpha();
    uR[k]     = cellRight.getPhase(k)->getVelocity().getX();
    rhoR[k]   = cellRight.getPhase(k)->getDensity();
    pR[k]     = cellRight.getPhase(k)->getPressure();
    HR[k]     = cellRight.getPhase(k)->getEos()->computeTotalEnthalpy(rhoR[k], pR[k], cellRight.getPhase(k)->getVelocity());
    cR[k]     = cellRight.getPhase(k)->getSoundSpeed();
    ER[k]     = cellRight.getPhase(k)->getTotalEnergy();

    if (std::fabs(uL[k]) < 1.e-12) uL[k] = 0.;
    if (std::fabs(uR[k]) < 1.e-12) uR[k] = 0.;
  }

  // Speed wave estimations with Roe means
  // -------------------------------------
  for (int k = 0; k < numberPhases; k++) {
    gamma = cellLeft.getPhase(k)->getEos()->getGamma();
    // Enthalpy estimation
    HMean = (std::sqrt(rhoL[k]) * HL[k] + std::sqrt(rhoR[k]) * HR[k]) / (std::sqrt(rhoL[k]) + std::sqrt(rhoR[k]));
    // Velocity estimation
    uMean = (std::sqrt(rhoL[k]) * uL[k] + std::sqrt(rhoR[k]) * uR[k]) / (std::sqrt(rhoL[k]) + std::sqrt(rhoR[k]));
    // Sound speed estimations (only for SG)
    cMean = std::sqrt((gamma - 1.) * (HMean - 0.5 * uMean * uMean));

    sL[k] = std::min(uL[k] - cL[k], uMean - cMean);
    sR[k] = std::max(uR[k] + cR[k], uMean + cMean);
  }

  for (int k = 0; k < numberPhases; k++) {
    if (std::fabs(sL[k]) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL[k]));
    if (std::fabs(sR[k]) > 1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR[k]));
  }

  // Coefficients to compute star states
  // -----------------------------------
  for (int k = 0; k < numberPhases; k++) {
    coefL[k] = rhoL[k] * (sL[k] - uL[k]);
    coefR[k] = rhoR[k] * (sR[k] - uR[k]);
  }

  // Pressure estimations of star states
  // -----------------------------------
  RHS1 = pL[1] / coefL[1] - uL[1] - pR[1] / coefR[1] + uR[1];
  RHS2 = alphaL[0] * (pL[0] / coefL[0] - uL[0] - pL[1] / coefL[1] + uL[1]) - alphaR[0] * (pR[0] / coefR[0] - uR[0] - pR[1] / coefR[1] + uR[1]);

  a = 1. / coefL[1];
  b = -1. / coefR[1];
  c = alphaL[0] / coefL[0] - alphaR[0] / coefR[0];
  d = -alphaL[0] / coefL[1];
  e = alphaR[0] / coefR[1];
  f = alphaL[0] - alphaR[0];
  g = alphaL[1];
  h = -alphaR[1];

  pStarR_1 = (-f * RHS2 / (h * c) + f * d * RHS1 / (a * c * h) - g * RHS1 / (h * a)) /
             (1. + f * d * b / (a * c * h) - f * e / (c * h) - g * b / (h * a));
  pStar_0  = RHS2 / c - d / (a * c) * (RHS1 - b * pStarR_1) - e * pStarR_1 / c;
  pStarL_1 = (RHS1 - b * pStarR_1) / a;

  // Estimations of velocities of star states
  // ----------------------------------------
  uStarL[0] = uL[0] + (pStar_0 - pL[0]) / coefL[0];
  uStarR[0] = uR[0] + (pStar_0 - pR[0]) / coefR[0];

  uStarL[1] = uL[1] + (pStarL_1 - pL[1]) / coefL[1];
  uStarR[1] = uR[1] + (pStarR_1 - pR[1]) / coefR[1];

  // Estimations of densities of star states
  // ---------------------------------------
  for (int k = 0; k < numberPhases; k++) {
    rhoStarL[k] = rhoL[k] * (sL[k] - uL[k]) / (sL[k] - uStarL[k]);
    rhoStarR[k] = rhoR[k] * (sR[k] - uR[k]) / (sR[k] - uStarR[k]);
  }

  // Configuration of the Riemann problem
  // ------------------------------------
  uI = uStarL[1];               // = sM_0
  if (uStarL[0] >= uStarL[1]) { // configuration A
    sM_0 = uStarR[0];
  }
  else { // configuration B
    sM_0    = uStarL[0];
    configA = false;
  }

  // Correction of contact discontinuities speed
  // -------------------------------------------
  if (std::fabs(uI) < 1.e-12) uI = 0.;
  if (std::fabs(sM_0) < 1.e-12) sM_0 = 0.;

  // Sampling the solution
  // --------------------
  for (int k = 0; k < numberPhases; k++) {
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alpha[k]  = 0.;
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaL[k] = alphaL[k];
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_alphaR[k] = alphaR[k];
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[k].setY(0.);
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[k].setZ(0.);
  }

  // Phase 1 in the paper
  if (sL[0] >= 0.) { //Supersonic left
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[0] = alphaL[0] * rhoL[0] * uL[0];
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].setX(alphaL[0] * (rhoL[0] * uL[0] * uL[0] + pL[0]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[0] = alphaL[0] * uL[0] * (rhoL[0] * EL[0] + pL[0]);
  }
  else if (sR[0] <= 0.) { //Supersonic right
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[0] = alphaR[0] * rhoR[0] * uR[0];
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].setX(alphaR[0] * (rhoR[0] * uR[0] * uR[0] + pR[0]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[0] = alphaR[0] * uR[0] * (rhoR[0] * ER[0] + pR[0]);
  }
  else if ((sL[0] <= 0. && uI >= 0. && configA) || (sL[0] <= 0. && sM_0 >= 0. && !configA)) { //Star left
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[0] = alphaL[0] * rhoL[0] * uL[0] +
                                                      sL[0] * (alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0) - alphaL[0] * rhoL[0]);
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].setX(
      alphaL[0] * (rhoL[0] * uL[0] * uL[0] + pL[0]) +
      sL[0] * (alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0) * sM_0 - alphaL[0] * rhoL[0] * uL[0]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[0] = alphaL[0] * uL[0] * (rhoL[0] * EL[0] + pL[0]) +
                                                       sL[0] * (alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0) *
                                                                  (EL[0] + ((sM_0 - uL[0]) * (sM_0 + pL[0] / (rhoL[0] * (sL[0] - uL[0]))))) -
                                                                alphaL[0] * rhoL[0] * EL[0]);
  }
  else if ((sM_0 < 0. && sR[0] >= 0. && configA) || (uI < 0. && sR[0] >= 0. && !configA)) { //Star right
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[0] = alphaR[0] * rhoR[0] * uR[0] +
                                                      sR[0] * (alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0) - alphaR[0] * rhoR[0]);
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].setX(
      alphaR[0] * (rhoR[0] * uR[0] * uR[0] + pR[0]) +
      sR[0] * (alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0) * sM_0 - alphaR[0] * rhoR[0] * uR[0]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[0] = alphaR[0] * uR[0] * (rhoR[0] * ER[0] + pR[0]) +
                                                       sR[0] * (alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0) *
                                                                  (ER[0] + ((sM_0 - uR[0]) * (sM_0 + pR[0] / (rhoR[0] * (sR[0] - uR[0]))))) -
                                                                alphaR[0] * rhoR[0] * ER[0]);
  }
  else if (configA) { //Star middle config A
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[0] = alphaR[0] * rhoR[0] * uR[0] +
                                                      sR[0] * (alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0) - alphaR[0] * rhoR[0]);
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].setX(
      alphaR[0] * (rhoR[0] * uR[0] * uR[0] + pR[0]) +
      sR[0] * (alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0) * sM_0 - alphaR[0] * rhoR[0] * uR[0]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[0] = alphaR[0] * uR[0] * (rhoR[0] * ER[0] + pR[0]) +
                                                       sR[0] * (alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0) *
                                                                  (ER[0] + ((sM_0 - uR[0]) * (sM_0 + pR[0] / (rhoR[0] * (sR[0] - uR[0]))))) -
                                                                alphaR[0] * rhoR[0] * ER[0]);

    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[0] += sM_0 * (alphaR[0] * rhoStarL[0] - alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].setX(
      static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].getX() +
      sM_0 * (alphaR[0] * rhoStarL[0] * uStarR[0] - alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0) * sM_0));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[0] +=
      sM_0 *
      (alphaR[0] * (rhoStarL[0] * cellLeft.getPhase(0)->getEos()->computeEnergy(rhoStarL[0], pStar_0) + 0.5 * rhoStarL[0] * uStarR[0] * uStarR[0]) -
       alphaR[0] * rhoR[0] * (sR[0] - uR[0]) / (sR[0] - sM_0) * (ER[0] + ((sM_0 - uR[0]) * (sM_0 + pR[0] / (rhoR[0] * (sR[0] - uR[0]))))));
  }
  else { //Star middle config B
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[0] = alphaL[0] * rhoL[0] * uL[0] +
                                                      sL[0] * (alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0) - alphaL[0] * rhoL[0]);
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].setX(
      alphaL[0] * (rhoL[0] * uL[0] * uL[0] + pL[0]) +
      sL[0] * (alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0) * sM_0 - alphaL[0] * rhoL[0] * uL[0]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[0] = alphaL[0] * uL[0] * (rhoL[0] * EL[0] + pL[0]) +
                                                       sL[0] * (alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0) *
                                                                  (EL[0] + ((sM_0 - uL[0]) * (sM_0 + pL[0] / (rhoL[0] * (sL[0] - uL[0]))))) -
                                                                alphaL[0] * rhoL[0] * EL[0]);

    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[0] += sM_0 * (alphaL[0] * rhoStarR[0] - alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].setX(
      static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[0].getX() +
      sM_0 * (alphaL[0] * rhoStarR[0] * uStarL[0] - alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0) * sM_0));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[0] +=
      sM_0 *
      (alphaL[0] * (rhoStarR[0] * cellLeft.getPhase(0)->getEos()->computeEnergy(rhoStarR[0], pStar_0) + 0.5 * rhoStarR[0] * uStarL[0] * uStarL[0]) -
       alphaL[0] * rhoL[0] * (sL[0] - uL[0]) / (sL[0] - sM_0) * (EL[0] + ((sM_0 - uL[0]) * (sM_0 + pL[0] / (rhoL[0] * (sL[0] - uL[0]))))));
  }

  // Phase 2 in the paper
  if (sL[1] >= 0.) { //Supersonic left
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[1] = alphaL[1] * rhoL[1] * uL[1];
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[1].setX(alphaL[1] * (rhoL[1] * uL[1] * uL[1] + pL[1]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[1] = alphaL[1] * uL[1] * (rhoL[1] * EL[1] + pL[1]);
  }
  else if (sR[1] <= 0.) { //Supersonic right
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[1] = alphaR[1] * rhoR[1] * uR[1];
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[1].setX(alphaR[1] * (rhoR[1] * uR[1] * uR[1] + pR[1]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[1] = alphaR[1] * uR[1] * (rhoR[1] * ER[1] + pR[1]);
  }
  else if (sL[1] <= 0. && uI >= 0.) { //Star left
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[1] = alphaL[1] * rhoL[1] * uL[1] +
                                                      sL[1] * (alphaL[1] * rhoL[1] * (sL[1] - uL[1]) / (sL[1] - uI) - alphaL[1] * rhoL[1]);
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[1].setX(
      alphaL[1] * (rhoL[1] * uL[1] * uL[1] + pL[1]) +
      sL[1] * (alphaL[1] * rhoL[1] * (sL[1] - uL[1]) / (sL[1] - uI) * uI - alphaL[1] * rhoL[1] * uL[1]));
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[1] = alphaL[1] * uL[1] * (rhoL[1] * EL[1] + pL[1]) +
                                                       sL[1] * (alphaL[1] * rhoL[1] * (sL[1] - uL[1]) / (sL[1] - uI) *
                                                                  (EL[1] + (uI - uL[1]) * (uI + pL[1] / (rhoL[1] * (sL[1] - uL[1])))) -
                                                                alphaL[1] * rhoL[1] * EL[1]);
  }
  else { //Star right
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_mass[1] = alphaR[1] * rhoR[1] * uR[1] +
                                                      sR[1] * (alphaR[1] * rhoR[1] * (sR[1] - uR[1]) / (sR[1] - uI) - alphaR[1] * rhoR[1]);
    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_momentum[1].setX(
      alphaR[1] * (rhoR[1] * uR[1] * uR[1] + pR[1]) +
      sR[1] * (alphaR[1] * rhoR[1] * (sR[1] - uR[1]) / (sR[1] - uI) * uI - alphaR[1] * rhoR[1] * uR[1]));

    static_cast<FluxNEqHLLAC*>(fluxBuff)->m_energ[1] = alphaR[1] * uR[1] * (rhoR[1] * ER[1] + pR[1]) +
                                                       sR[1] * (alphaR[1] * rhoR[1] * (sR[1] - uR[1]) / (sR[1] - uI) *
                                                                  (ER[1] + (uI - uR[1]) * (uI + pR[1] / (rhoR[1] * (sR[1] - uR[1])))) -
                                                                alphaR[1] * rhoR[1] * ER[1]);
  }

  // Contact discontinuity velocity
  // ------------------------------
  static_cast<FluxNEqHLLAC*>(fluxBuff)->m_uI     = uI;
  static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarL = pStarL_1;
  static_cast<FluxNEqHLLAC*>(fluxBuff)->m_pStarR = pStarR_1;
}

//****************************************************************************
