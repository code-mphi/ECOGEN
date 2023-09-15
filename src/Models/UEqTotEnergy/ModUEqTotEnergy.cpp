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

#include <cmath>
#include <algorithm>
#include "ModUEqTotEnergy.h"
#include "PhaseUEqTotEnergy.h"

const std::string ModUEqTotEnergy::NAME = "VELOCITYEQTOTENERGY";

//***********************************************************************

ModUEqTotEnergy::ModUEqTotEnergy(const int& numbTransports, const int& numbPhases) : Model(NAME, numbTransports)
{
  fluxBuff = new FluxUEqTotEnergy(numbPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxUEqTotEnergy(numbPhases));
  }
}

//***********************************************************************

ModUEqTotEnergy::ModUEqTotEnergy(const std::string& name, const int& numbTransports) : Model(name, numbTransports){}

//***********************************************************************

ModUEqTotEnergy::~ModUEqTotEnergy()
{
  delete fluxBuff;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//***********************************************************************

void ModUEqTotEnergy::allocateCons(Flux** cons)
{
  *cons = new FluxUEqTotEnergy(numberPhases);
}

//***********************************************************************

void ModUEqTotEnergy::allocatePhase(Phase** phase)
{
  *phase = new PhaseUEqTotEnergy;
}

//***********************************************************************

void ModUEqTotEnergy::allocateMixture(Mixture** mixture)
{
  *mixture = new MixUEqTotEnergy;
}

//***********************************************************************

void ModUEqTotEnergy::fulfillState(Phase** phases, Mixture* mixture)
{
  //Complete phases state
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  //Complete mixture variables using phases variable
  mixture->computeMixtureVariables(phases);
}

//***********************************************************************

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModUEqTotEnergy::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double> &boundData) const
{
  Phase* vecPhase;
  double sL, sR;
  double pStar(0.), rhoStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getFrozenSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uR = cellRight.getMixture()->getVelocity().getX(), cR = cellRight.getMixture()->getFrozenSoundSpeed(), pR = cellRight.getMixture()->getPressure(), rhoR = cellRight.getMixture()->getDensity();

  //Davies
  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);

  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR)), mkL, mkR;
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (std::fabs(sM)<1.e-8) sM = 0.;

  //Solution sampling
  if (sL >= 0.){
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellLeft.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double density = vecPhase->getDensity();
      double energie = vecPhase->getTotalEnergy();
      double pressure = vecPhase->getPressure();
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_mass[k] = alpha*density*uL;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_totEnerg[k] = alpha*density*energie*uL + alpha*pressure*uL;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alphap[k] = alpha*pressure;
    }
    double vitY = cellLeft.getMixture()->getVelocity().getY(); double vitZ = cellLeft.getMixture()->getVelocity().getZ();
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setX(rhoL*uL*uL + pL);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setY(rhoL*vitY*uL);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setZ(rhoL*vitZ*uL);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_uStar = uL;
    
    // Boundary data for output
    boundData[VarBoundary::p] = pL;
    boundData[VarBoundary::rho] = rhoL;
    boundData[VarBoundary::velU] = uL;
    boundData[VarBoundary::velV] = vitY;
    boundData[VarBoundary::velW] = vitZ;
  }
  else if (sR <= 0.){
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellRight.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double density = vecPhase->getDensity();
      double energie = vecPhase->getTotalEnergy();
      double pressure = vecPhase->getPressure();
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_mass[k] = alpha*density*uR;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_totEnerg[k] = alpha*density*energie*uR + alpha*pressure*uR;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alphap[k] = alpha*pressure;
    }
    double vitY = cellRight.getMixture()->getVelocity().getY(); double vitZ = cellRight.getMixture()->getVelocity().getZ();
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setX(rhoR*uR*uR + pR);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setY(rhoR*vitY*uR);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setZ(rhoR*vitZ*uR);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_uStar = uR;

    // Boundary data for output
    boundData[VarBoundary::p] = pR;
    boundData[VarBoundary::rho] = rhoR;
    boundData[VarBoundary::velU] = uR;
    boundData[VarBoundary::velV] = vitY;
    boundData[VarBoundary::velW] = vitZ;
  }
  else if (sM >= 0.){
    //Compute left solution state
    double vitY = cellLeft.getMixture()->getVelocity().getY(); double vitZ = cellLeft.getMixture()->getVelocity().getZ();
    rhoStar = mL / (sL - sM);
    pStar = mL*(sM - uL) + pL;
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellLeft.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double density = vecPhase->getDensity();
      double pressure = vecPhase->getPressure();
      mkL = density*(sL - uL);
      TB->rhokStar[k] = mkL / (sL - sM);
      TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
      TB->pkStar[k] = TB->eos[k]->computePressureIsentropic(pressure, density, TB->rhokStar[k]);
      //TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
      TB->EkStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]) + 0.5 * (sM*sM + vitY*vitY + vitZ*vitZ);
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_mass[k] = alpha* TB->rhokStar[k] * sM;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_totEnerg[k] = alpha* TB->rhokStar[k] * TB->EkStar[k] * sM + alpha*TB->pkStar[k]*sM;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alphap[k] = alpha* TB->pkStar[k];
    }
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setY(rhoStar*vitY*sM);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setZ(rhoStar*vitZ*sM);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_uStar = sM;

    // Boundary data for output
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vitY;
    boundData[VarBoundary::velW] = vitZ;                                    
  }
  else {
    //Compute right solution state
    double vitY = cellRight.getMixture()->getVelocity().getY(); double vitZ = cellRight.getMixture()->getVelocity().getZ();
    rhoStar = mR / (sR - sM);
    pStar = mR*(sM - uR) + pR;
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellRight.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double density = vecPhase->getDensity();
      double pressure = vecPhase->getPressure();
      mkR = density*(sR - uR);
      TB->rhokStar[k] = mkR / (sR - sM);
      TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
      TB->pkStar[k] = TB->eos[k]->computePressureIsentropic(pressure, density, TB->rhokStar[k]);
      //TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
      TB->EkStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]) + 0.5 * (sM*sM + vitY*vitY + vitZ*vitZ);
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_mass[k] = alpha* TB->rhokStar[k] * sM;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_totEnerg[k] = alpha* TB->rhokStar[k] * TB->EkStar[k] * sM + alpha*TB->pkStar[k]*sM;
      static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_alphap[k] = alpha* TB->pkStar[k];
    }
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setY(rhoStar*vitY*sM);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.setZ(rhoStar*vitZ*sM);
    static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_uStar = sM;

    // Boundary data for output
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vitY;
    boundData[VarBoundary::velW] = vitZ;
  }

  //Contact discontinuity velocity
  static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModUEqTotEnergy::selectScalar(Phase** phases, Mixture* mixture, Transport* transports, Variable nameVariable, int num) const
{
  switch (nameVariable) {
    case Variable::pressure:
      if (num < 0) {
        return mixture->getPressure();
      }
      else {
        return phases[num]->getPressure();
      }
      break;
    case Variable::density:
      if (num < 0) {
        return mixture->getDensity();
      }
      else {
        return phases[num]->getDensity();
      }
      break;
    case Variable::alpha:
      return phases[num]->getAlpha();
      break;
    case Variable::velocityU:
      return mixture->getVelocity().getX();
      break;
    case Variable::velocityV:
      return mixture->getVelocity().getY();
      break;
    case Variable::velocityW:
      return mixture->getVelocity().getZ();
      break;
    case Variable::velocityMag:
      return mixture->getVelocity().norm();
      break;
    case Variable::transport:
      return transports[num].getValue();
      break;
    case Variable::temperature:
      return phases[num]->getTemperature();
      break;
    default:
      Errors::errorMessage("nameVariable unknown in selectScalar"); return 0;
      break;
  }
}

//****************************************************************************

const double& ModUEqTotEnergy::getSM()
{
  return static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_sM;
}

//****************************************************************************
//***************************** others methods *******************************
//****************************************************************************

void ModUEqTotEnergy::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  static_cast<FluxUEqTotEnergy*> (fluxBuff)->m_momentum.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************