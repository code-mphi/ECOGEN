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

#include "ModNEqHLLC4Waves.h"

const std::string ModNEqHLLC4Waves::NAME = "NONEQ";

//***********************************************************************

ModNEqHLLC4Waves::ModNEqHLLC4Waves(const int& numbTransports, const int& numbPhases) : ModNEq(NAME, numbTransports)
{
  fluxBuff    = new FluxNEqHLLC4Waves(numbPhases);
  fluxBuffMRF = new FluxNEqHLLC4Waves(numbPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxNEqHLLC4Waves(numbPhases));
  }
}

//***********************************************************************

ModNEqHLLC4Waves::~ModNEqHLLC4Waves() {}

//***********************************************************************

void ModNEqHLLC4Waves::allocateCons(Flux** cons) { *cons = new FluxNEqHLLC4Waves(numberPhases); }

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModNEqHLLC4Waves::solver(
  Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/) const
{
  double alphaStar(0.), pStar(0.), rhoStar(0.), EStar(0.), uStar(0.), vStar(0.), wStar(0.);
  double pL(0.), uL(0.), rhoL(0.), cL(0.), alphaL(0.), EL(0.);
  double pR(0.), uR(0.), rhoR(0.), cR(0.), alphaR(0.), ER(0.);
  double sL(0.), sR(0.), sM(0.), sLmax(1.e5), sRmax(-1.e5);

  //uI and pI sampling depending on the left and right volume fraction proportions
  alphaL = cellLeft.getPhase(0)->getAlpha();
  alphaR = cellRight.getPhase(0)->getAlpha();
  int phaseL(0), phaseR(1); //default contact 1-2
  if (alphaL < alphaR) {    //contact 2-1
    phaseL = 1;
    phaseR = 0;
  }
  pL   = cellLeft.getPhase(phaseL)->getPressure();
  uL   = cellLeft.getPhase(phaseL)->getVelocity().getX();
  rhoL = cellLeft.getPhase(phaseL)->getDensity();
  cL   = cellLeft.getPhase(phaseL)->getSoundSpeed();
  pR   = cellRight.getPhase(phaseR)->getPressure();
  uR   = cellRight.getPhase(phaseR)->getVelocity().getX();
  rhoR = cellRight.getPhase(phaseR)->getDensity();
  cR   = cellRight.getPhase(phaseR)->getSoundSpeed();
  sL   = std::min(uL - cL, cellRight.getPhase(phaseL)->getVelocity().getX() - cellRight.getPhase(phaseL)->getSoundSpeed());
  sR   = std::max(cellLeft.getPhase(phaseR)->getVelocity().getX() + cellLeft.getPhase(phaseR)->getSoundSpeed(), uR + cR);

  double uI(0.), pI(0.);
  uI = (rhoR * uR * uR + pR - rhoL * uL * uL - pL + sL * rhoL * uL - sR * rhoR * uR) / (rhoR * uR - rhoL * uL + sL * rhoL - sR * rhoR);
  pI = rhoR * (uR - sR) * (uR - uI) + pR;

  //Determine sR and sL max
  for (int k = 0; k < numberPhases; k++) {
    uL    = cellLeft.getPhase(k)->getVelocity().getX();
    cL    = cellLeft.getPhase(k)->getSoundSpeed();
    uR    = cellRight.getPhase(k)->getVelocity().getX();
    cR    = cellRight.getPhase(k)->getSoundSpeed();
    sL    = std::min(uL - cL, uR - cR);
    sR    = std::max(uL + cL, uR + cR);
    sLmax = std::min(sL, sLmax);
    sRmax = std::max(sR, sRmax);
  }

  //Loop of "independant" Riemann problems for phases
  for (int k = 0; k < numberPhases; k++) {
    alphaL = cellLeft.getPhase(k)->getAlpha();
    alphaR = cellRight.getPhase(k)->getAlpha();
    pL     = cellLeft.getPhase(k)->getPressure();
    uL     = cellLeft.getPhase(k)->getVelocity().getX();
    rhoL   = cellLeft.getPhase(k)->getDensity();
    cL     = cellLeft.getPhase(k)->getSoundSpeed();
    EL     = cellLeft.getPhase(k)->getTotalEnergy();
    pR     = cellRight.getPhase(k)->getPressure();
    uR     = cellRight.getPhase(k)->getVelocity().getX();
    rhoR   = cellRight.getPhase(k)->getDensity();
    cR     = cellRight.getPhase(k)->getSoundSpeed();
    ER     = cellRight.getPhase(k)->getTotalEnergy();
    sL     = std::min(uL - cL, uR - cR);
    sR     = std::max(uL + cL, uR + cR);
    if (uI < 0.9 * sL) sL = sLmax;
    if (uI > 0.9 * sR) sR = sRmax;
    sM = (alphaR * (rhoR * uR * uR + pR) - alphaL * (rhoL * uL * uL + pL) + sL * alphaL * rhoL * uL - sR * alphaR * rhoR * uR +
          (alphaL - alphaR) * pI) /
         (alphaR * rhoR * uR - alphaL * rhoL * uL + sL * alphaL * rhoL - sR * alphaR * rhoR);
    if (std::fabs(sM) < 1.e-8) sM = 0.;

    if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
    if (std::fabs(sR) > 1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR));

    //Sampling
    if (sL > 0.) { //supersonic, left state solution
      double vitY = cellLeft.getPhase(k)->getVelocity().getY();
      double vitZ = cellLeft.getPhase(k)->getVelocity().getZ();

      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_alpha[k] = alphaL * uI;
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_mass[k]  = alphaL * rhoL * uL;
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setX(alphaL * (rhoL * uL * uL + pL - pI));
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setY(alphaL * rhoL * vitY * uL);
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setZ(alphaL * rhoL * vitZ * uL);
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_energ[k] = alphaL * ((rhoL * EL + pL) * uL - pI * uI);
    }
    else if (sR < 0.) { //supersonic, right state solution
      double vitY = cellRight.getPhase(k)->getVelocity().getY();
      double vitZ = cellRight.getPhase(k)->getVelocity().getZ();

      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_alpha[k] = alphaR * uI;
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_mass[k]  = alphaR * rhoR * uR;
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setX(alphaR * (rhoR * uR * uR + pR - pI));
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setY(alphaR * rhoR * vitY * uR);
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setZ(alphaR * rhoR * vitZ * uR);
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_energ[k] = alphaR * ((rhoR * ER + pR) * uR - pI * uI);
    }
    else {
      uStar = sM;
      if (sM < uI) {
        if (sM > 0.) { //State left star solution
          alphaStar = alphaL;
          rhoStar   = rhoL * (uL - sL) / (sM - sL);
          pStar     = pL + rhoL * (uL - sL) * (uL - sM);
          EStar     = EL + (pL * uL - pStar * sM) / rhoL / (uL - sL);
          vStar     = cellLeft.getPhase(k)->getVelocity().getY();
          wStar     = cellLeft.getPhase(k)->getVelocity().getZ();
        }
        else if (uI > 0.) { //State star star solution
          alphaStar            = alphaL;
          double alphaRhoStarR = alphaR * rhoR * (uR - sR) / (sM - sR);
          rhoStar              = alphaRhoStarR / alphaStar;
          pStar                = pL + rhoL * (uL - sL) * (uL - sM);
          double pStarR        = pR + rhoR * (uR - sR) * (uR - sM);
          EStar                = (ER + (pR * uR - pStarR * sM) / rhoR / (uR - sR)) - (alphaL - alphaR) / alphaRhoStarR * pI;
          vStar                = cellRight.getPhase(k)->getVelocity().getY();
          wStar                = cellRight.getPhase(k)->getVelocity().getZ();
        }
        else { //State right star solution
          alphaStar = alphaR;
          rhoStar   = rhoR * (uR - sR) / (sM - sR);
          pStar     = pR + rhoR * (uR - sR) * (uR - sM);
          EStar     = ER + (pR * uR - pStar * sM) / rhoR / (uR - sR);
          vStar     = cellRight.getPhase(k)->getVelocity().getY();
          wStar     = cellRight.getPhase(k)->getVelocity().getZ();
        }
      } // end sM < uI
      else {
        if (uI > 0.) { //State left star solution
          alphaStar = alphaL;
          rhoStar   = rhoL * (uL - sL) / (sM - sL);
          pStar     = pL + rhoL * (uL - sL) * (uL - sM);
          EStar     = EL + (pL * uL - pStar * sM) / rhoL / (uL - sL);
          vStar     = cellLeft.getPhase(k)->getVelocity().getY();
          wStar     = cellLeft.getPhase(k)->getVelocity().getZ();
        }
        else if (sM > 0.) { //State star star solution
          alphaStar            = alphaR;
          double alphaRhoStarL = alphaL * rhoL * (uL - sL) / (sM - sL);
          rhoStar              = alphaRhoStarL / alphaStar;
          pStar                = pR + rhoR * (uR - sR) * (uR - sM);
          double pStarL        = pL + rhoL * (uL - sL) * (uL - sM);
          EStar                = (EL + (pL * uL - pStarL * sM) / rhoL / (uL - sL)) - (alphaR - alphaL) / alphaRhoStarL * pI;
          vStar                = cellLeft.getPhase(k)->getVelocity().getY();
          wStar                = cellLeft.getPhase(k)->getVelocity().getZ();
        }
        else { //State right star solution
          alphaStar = alphaR;
          rhoStar   = rhoR * (uR - sR) / (sM - sR);
          pStar     = pR + rhoR * (uR - sR) * (uR - sM);
          EStar     = ER + (pR * uR - pStar * sM) / rhoR / (uR - sR);
          vStar     = cellRight.getPhase(k)->getVelocity().getY();
          wStar     = cellRight.getPhase(k)->getVelocity().getZ();
        }
      }
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_alpha[k] = alphaStar * uI;
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_mass[k]  = alphaStar * rhoStar * uStar;
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setX(alphaStar * (rhoStar * uStar * uStar + pStar - pI));
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setY(alphaStar * rhoStar * vStar * uStar);
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_momentum[k].setZ(alphaStar * rhoStar * wStar * uStar);
      static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_energ[k] = alphaStar * ((rhoStar * EStar + pStar) * uStar - pI * uI);
    }
  }

  //Contact discontinuity velocity
  static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_uI = uI;
  static_cast<FluxNEqHLLC4Waves*>(fluxBuff)->m_pI = pI;
}

//****************************************************************************
