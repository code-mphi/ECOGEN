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

#include "RiemannSolverExact.h"

//***********************************************************************

RiemannSolverExact::RiemannSolverExact() : RiemannSolver() {}

//***********************************************************************

RiemannSolverExact::~RiemannSolverExact() {}

//***********************************************************************

void RiemannSolverExact::solver(
  const Phase& phaseLeft, const Phase& phaseRight, const double& dxLeft, const double& dxRight, double& dtMax, double& pStar, FluxEuler& m_fluxEuler)
{
  // WARNING: Only written for IG and SG EOS (not generalized)

  double cL, cR, gammaL, gammaR, pInfL, pInfR, uStar, uStarL(1.), uStarR(2.);
  double uL, uR, vL, vR, wL, wR, pL, pR, rhoL, rhoR, waveFrontR(0.), waveFrontL(0.), rhoStarR, rhoStarL, cStarL, cStarR;
  double eRefFlux(0.), eRefL, eRefR;
  double magnitudeR(0.), magnitudeL(0.);
  double coef;

  // Initialisation of the left state
  // -----------------------------------------------------------
  uL     = phaseLeft.getU();
  vL     = phaseLeft.getV();
  wL     = phaseLeft.getW();
  pL     = phaseLeft.getPressure();
  rhoL   = phaseLeft.getDensity();
  cL     = phaseLeft.getSoundSpeed();
  eRefL  = phaseLeft.getEos()->getERef();
  gammaL = phaseLeft.getEos()->getGamma();
  pInfL  = phaseLeft.getEos()->getPInf();

  // Initialisation of the right state
  // -----------------------------------------------------------
  uR     = phaseRight.getU();
  vR     = phaseRight.getV();
  wR     = phaseRight.getW();
  pR     = phaseRight.getPressure();
  rhoR   = phaseRight.getDensity();
  cR     = phaseRight.getSoundSpeed();
  eRefR  = phaseRight.getEos()->getERef();
  gammaR = phaseRight.getEos()->getGamma();
  pInfR  = phaseRight.getEos()->getPInf();

  // Initial calculus of the pressure of the intermediate states
  // -----------------------------------------------------------
  double zL = rhoL * cL;
  double zR = rhoR * cR;
  pStar     = (zL * pR + zR * pL + zL * zR * (uL - uR)) / (zL + zR);

  double pmin = std::max(-pInfL, -pInfR) + 1.e-3;
  if (pStar <= pmin) {
    pStar = pmin;
  }

  // Determination of the star pressure by using Newton algorithm
  // ------------------------------------------------------------
  double deltaUStar = uStarR - uStarL, absDeltaUStar = std::fabs(deltaUStar), dDeltaUStar;
  double epsilon(1.e-6);

  int maxIter = 100;
  int iter    = 1;
  double waveSpeedL, dWaveSpeedL, waveSpeedR, dWaveSpeedR, duStarL, duStarR;

  while (absDeltaUStar > epsilon && iter < maxIter) {
    computeWaveSpeed(pL, rhoL, cL, gammaL, pInfL, pStar, waveSpeedL, dWaveSpeedL);
    uStarL  = uL - waveSpeedL;
    duStarL = -dWaveSpeedL;

    computeWaveSpeed(pR, rhoR, cR, gammaR, pInfR, pStar, waveSpeedR, dWaveSpeedR);
    uStarR  = uR + waveSpeedR;
    duStarR = dWaveSpeedR;

    deltaUStar    = uStarR - uStarL;
    absDeltaUStar = std::fabs(deltaUStar);

    if (absDeltaUStar >= epsilon) {
      dDeltaUStar = duStarR - duStarL;
      pStar       = pStar - deltaUStar / dDeltaUStar;
      if (pStar <= pmin) {
        pStar = pmin;
      }
    }
    iter++;
  }

  if (iter == maxIter) {
    warnings.push_back(Errors("Convergence issue of the exact Riemann solver in ModNEqDEM::EulerExact.", __FILE__, __LINE__));
  }

  // Determination of the star velocity
  // ----------------------------------
  uStar = 0.5 * (uStarL + uStarR);
  if (std::fabs(uStar) <= epsilon) {
    uStar = 0.;
  }

  // Pre-calculations and determination of dtMax
  // ------------------------------------------

  double rapR = (pStar + pInfR) / (pR + pInfR);
  double rapL = (pStar + pInfL) / (pL + pInfL);

  // Left shock wave
  // ---------------
  if (rapL > 1.) {
    computeCoefficientRatio(rapL, gammaL, coef);
    magnitudeL = std::sqrt(gammaL * rhoL * (pL + pInfL)) * coef;
    waveFrontL = uL - magnitudeL / rhoL;
    rhoStarL   = magnitudeL / (uStar - waveFrontL);
    cStarL     = std::sqrt(gammaL * (pStar + pInfL) / rhoStarL);
    dtMax      = std::fmin(dxLeft / std::fabs(waveFrontL), dtMax);
  }
  // Left rarefaction wave
  // ---------------------
  else {
    rhoStarL = rhoL * std::pow(((pStar + pInfL) / (pL + pInfL)), 1. / gammaL);
    cStarL   = std::sqrt(gammaL * (pStar + pInfL) / rhoStarL);
    dtMax    = std::fmin(dxLeft / std::fabs(uL - cL), dtMax);
  }
  // Right shock wave
  // ----------------
  if (rapR > 1.) {
    computeCoefficientRatio(rapR, gammaR, coef);
    magnitudeR = std::sqrt(gammaR * rhoR * (pR + pInfR)) * coef;
    waveFrontR = uR + magnitudeR / rhoR;
    rhoStarR   = -magnitudeR / (uStar - waveFrontR);
    cStarR     = std::sqrt(gammaR * (pStar + pInfR) / rhoStarR);
    dtMax      = std::fmin(dxRight / std::fabs(waveFrontR), dtMax);
  }
  // Right rarefaction wave
  // ---------------------
  else {
    rhoStarR = rhoR * std::pow(((pStar + pInfR) / (pR + pInfR)), 1. / gammaR);
    cStarR   = std::sqrt(gammaR * (pStar + pInfR) / rhoStarR);
    dtMax    = std::fmin(dxRight / std::fabs(uR + cR), dtMax);
  }

  // Sampling
  // --------
  double pFlux(0.), cFlux(0.), uFlux(0.), rhoFlux(0.), gammaFlux(1.01), pinfFlux(0.), vFlux(0.), wFlux(0.);

  // Left shock wave
  // ---------------
  if (uStar > 0.) {
    if (pL <= pStar) {
      if (0. < waveFrontL) {
        pFlux   = pL;
        cFlux   = cL;
        uFlux   = uL;
        rhoFlux = rhoL;
      }
      else {
        pFlux   = pStar;
        cFlux   = cStarL;
        uFlux   = uStar;
        rhoFlux = rhoStarL;
      }
      gammaFlux = gammaL;
      pinfFlux  = pInfL;
      vFlux     = vL;
      wFlux     = wL;
      eRefFlux  = eRefL;
    }

    // Left rarefaction wave
    // ---------------------
    else {
      if (uL - cL >= 0.) {
        pFlux   = pL;
        cFlux   = cL;
        uFlux   = uL;
        rhoFlux = rhoL;
      }
      else {
        if (uStar - cStarL > 0.) {
          uFlux   = ((gammaL - 1.) * uL + 2. * cL) / (gammaL + 1.);
          cFlux   = uFlux;
          rhoFlux = std::pow((cFlux * cFlux * std::pow(rhoL, gammaL)) / ((pInfL + pL) * gammaL), 1. / (gammaL - 1.));
          pFlux   = cFlux * cFlux * rhoFlux / gammaL - pInfL;
        }
        else {
          pFlux   = pStar;
          cFlux   = cStarL;
          uFlux   = uStar;
          rhoFlux = rhoStarL;
        }
      }
      vFlux     = vL;
      wFlux     = wL;
      gammaFlux = gammaL;
      pinfFlux  = pInfL;
      eRefFlux  = eRefL;
    }
  }
  else {
    // Right shock wave
    // ----------------
    if (pStar >= pR) {

      if (waveFrontR >= 0.) {
        pFlux   = pStar;
        cFlux   = cStarR;
        uFlux   = uStar;
        rhoFlux = rhoStarR;
      }
      else {
        pFlux   = pR;
        cFlux   = cStarR;
        uFlux   = uR;
        rhoFlux = rhoR;
      }
      vFlux     = vR;
      wFlux     = wR;
      gammaFlux = gammaR;
      pinfFlux  = pInfR;
      eRefFlux  = eRefR;
    }
    // Right rarefaction wave
    // ----------------------
    else {
      if (uStar + cStarR >= 0.) {
        pFlux   = pStar;
        cFlux   = cStarR;
        uFlux   = uStar;
        rhoFlux = rhoStarR;
      }
      else {
        if (uR + cR > 0.) {
          uFlux   = ((gammaR - 1.) * uR + 2. * (-cR)) / (gammaR + 1.);
          cFlux   = -uFlux;
          rhoFlux = std::pow(cFlux * cFlux * std::pow(rhoR, gammaR) / ((pInfR + pR) * gammaR), (1. / (gammaR - 1.)));
          pFlux   = cFlux * cFlux * rhoFlux / gammaR - pInfR;
        }
        else {
          pFlux   = pR;
          cFlux   = cR;
          uFlux   = uR;
          rhoFlux = rhoR;
        }
      }
      vFlux     = vR;
      wFlux     = wR;
      gammaFlux = gammaR;
      pinfFlux  = pInfR;
      eRefFlux  = eRefR;
    }
  }

  // Flux
  // ----
  m_fluxEuler.m_mass = rhoFlux * uFlux;
  m_fluxEuler.m_momentum.setX(rhoFlux * uFlux * uFlux + pFlux);
  m_fluxEuler.m_momentum.setY(rhoFlux * uFlux * vFlux);
  m_fluxEuler.m_momentum.setZ(rhoFlux * uFlux * wFlux);
  m_fluxEuler.m_energ = uFlux * ((pFlux + gammaFlux * pinfFlux) / (gammaFlux - 1.) + rhoFlux * eRefFlux + rhoFlux * uFlux * uFlux / 2. + pFlux);

  // Interfacial condition
  // ----------------------
  m_fluxEuler.m_sM = uStar;
}

//****************************************************************************

void RiemannSolverExact::computeWaveSpeed(
  double p0, double rho0, double soundSpeed0, double gamma, double pInf, double targetPressure, double& waveSpeed, double& dWaveSpeed)
{
  double gammaPlus     = (gamma + 1.) / (2. * gamma);
  double gammaMinus    = (gamma - 1.) / (2. * gamma);
  double pressureRatio = (targetPressure + pInf) / (p0 + pInf);

  // Shock wave
  if (targetPressure > p0) {
    double numerator             = targetPressure - p0;
    double denominator           = rho0 * soundSpeed0 * sqrt(gammaPlus * pressureRatio + gammaMinus);
    waveSpeed                    = numerator / denominator;
    double denominatorDerivative = rho0 * soundSpeed0 * (gammaPlus / (p0 + pInf)) / (2. * sqrt(gammaPlus * pressureRatio + gammaMinus));
    dWaveSpeed                   = (denominator - numerator * denominatorDerivative) / (denominator * denominator);
  }
  // Rarefaction wave
  else {
    waveSpeed  = 2. * soundSpeed0 / (gamma - 1.) * (std::pow(pressureRatio, gammaMinus) - 1.);
    dWaveSpeed = 2. * soundSpeed0 / (gamma - 1.) * gammaMinus / (p0 + pInf) * std::pow(pressureRatio, -gammaPlus);
  }
}

//****************************************************************************

void RiemannSolverExact::computeCoefficientRatio(double pressureRatio, double gamma, double& coefficientRatio)
{
  double threshold = 1. - 1.e-8;
  double factor    = (gamma - 1.) / (2. * gamma);

  if (pressureRatio < threshold) {
    coefficientRatio = factor * (1. - pressureRatio) / (1. - std::pow(pressureRatio, factor));
  }
  else {
    coefficientRatio = sqrt((gamma + 1.) / 2. / gamma * pressureRatio + factor);
  }
}

//****************************************************************************
