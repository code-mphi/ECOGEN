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

#include "ModEulerHomogeneous.h"
#include "PhaseEulerHomogeneous.h"
#include "GradPhaseEulerHomogeneous.h"
#include "GradMixEulerHomogeneous.h"

const std::string ModEulerHomogeneous::NAME = "EULERHOMOGENEOUS";

//****************************************************************************

ModEulerHomogeneous::ModEulerHomogeneous(const int& numbTransports, const int liquid, const int vapor) :
  Model(NAME, numbTransports), m_liq(liquid), m_vap(vapor)
{
  fluxBuff = new FluxEulerHomogeneous(); 
}

//****************************************************************************

ModEulerHomogeneous::~ModEulerHomogeneous()
{
  delete fluxBuff;
}

//****************************************************************************

void ModEulerHomogeneous::allocateCons(Flux** cons)
{
  *cons = new FluxEulerHomogeneous();
}

//***********************************************************************

void ModEulerHomogeneous::allocatePhase(Phase** phase)
{
  *phase = new PhaseEulerHomogeneous;
}

//***********************************************************************

void ModEulerHomogeneous::allocateMixture(Mixture** mixture)
{
  *mixture = new MixEulerHomogeneous;
}

//***********************************************************************

void ModEulerHomogeneous::allocatePhaseGradient(GradPhase** phase)
{
  *phase = new GradPhaseEulerHomogeneous;
}

//***********************************************************************

void ModEulerHomogeneous::allocateMixtureGradient(GradMixture** mixture)
{
  *mixture = new GradMixEulerHomogeneous;
}

//***********************************************************************

void ModEulerHomogeneous::fulfillState(Phase** phases, Mixture* mixture)
{
  //Temperature calculus
  double Tsat = mixture->computeTsat(phases[m_liq]->getEos(), phases[m_vap]->getEos(), mixture->getPressure());

  //Complete phases and mixture states from : alphak, pressure, temperature
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->setPressure(mixture->getPressure());
    phases[k]->setDensity(phases[k]->getEos()->computeDensitySaturation(mixture->getPressure(), Tsat, Tsat));
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  mixture->setTemperature(Tsat);
  mixture->computeMixtureVariables(phases);
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModEulerHomogeneous::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double> &boundData) const
{
  double sL, sR;
  
  //FP//TODO//look for sound speed
  double uL = cellLeft.getMixture()->getVelocity().getX(), vL = cellLeft.getMixture()->getVelocity().getY(), wL = cellLeft.getMixture()->getVelocity().getZ(), cL = cellLeft.getMixture()->getMixSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uR = cellRight.getMixture()->getVelocity().getX(), vR = cellRight.getMixture()->getVelocity().getY(), wR = cellRight.getMixture()->getVelocity().getZ(), cR = cellRight.getMixture()->getMixSoundSpeed(), pR = cellRight.getMixture()->getPressure(), rhoR = cellRight.getMixture()->getDensity();
  double EL = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm(), ER = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();

  //Davis
  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (std::fabs(sM)<1.e-8) sM = 0.;

  if (sL > 0.){
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_mass = rhoL*uL;
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setX(rhoL*uL*uL + pL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setY(rhoL*vL*uL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setZ(rhoL*wL*uL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = (rhoL*EL + pL)*uL;

    // Boundary data for output
    boundData[VarBoundary::p] = pL;
    boundData[VarBoundary::rho] = rhoL;
    boundData[VarBoundary::velU] = uL;
    boundData[VarBoundary::velV] = vL;
    boundData[VarBoundary::velW] = wL;
  }
  else if (sR < 0.){
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_mass = rhoR*uR;
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setX(rhoR*uR*uR + pR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setY(rhoR*vR*uR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setZ(rhoR*wR*uR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = (rhoR*ER + pR)*uR;

    // Boundary data for output
    boundData[VarBoundary::p] = pR;
    boundData[VarBoundary::rho] = rhoR;
    boundData[VarBoundary::velU] = uR;
    boundData[VarBoundary::velV] = vR;
    boundData[VarBoundary::velW] = wR;
  }

  ////1) Option HLL
  //else if (std::fabs(sR - sL)>1.e-3)
  //{
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_mass = (rhoR*uR*sL - rhoL*uL*sR + sL*sR*(rhoL - rhoR)) / (sL - sR);
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setX(((rhoR*uR*uR + pR)*sL - (rhoL*uL*uL + pL)*sR + sL*sR*(rhoL*uL - rhoR*uR)) / (sL - sR));
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setY((rhoR*uR*vR*sL - rhoL*uL*vL*sR + sL*sR*(rhoL*vL - rhoR*vR)) / (sL - sR));
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setZ((rhoR*uR*wR*sL - rhoL*uL*wL*sR + sL*sR*(rhoL*wL - rhoR*wR)) / (sL - sR));
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = ((rhoR*ER + pR)*uR*sL - (rhoL*EL + pL)*uL*sR + sL*sR*(rhoL*EL - rhoR*ER)) / (sL - sR);
  //}

  //2) Option HLLC
  else if (sM >= 0.) {
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_mass = rhoStar*sM;
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM+pStar);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setY(rhoStar*sM*vL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setZ(rhoStar*sM*wL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;

    // Boundary data for output
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vL;
    boundData[VarBoundary::velW] = wL;
  }
  else {
    double pStar = mR*(sM - uR) + pR;
    double rhoStar = mR / (sR - sM);
    double Estar = ER + (sM - uR)*(sM + pR / mR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_mass = rhoStar*sM;
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setY(rhoStar*sM*vR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.setZ(rhoStar*sM*wR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;

    // Boundary data for output
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vR;
    boundData[VarBoundary::velW] = wR;
  }

  //Contact discontinuity velocity
  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions ****************
//****************************************************************************

//Not implemented yet

//****************************************************************************

const double& ModEulerHomogeneous::getSM()
{
  return static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_sM;
}

//****************************************************************************

void ModEulerHomogeneous::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_momentum.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************

int ModEulerHomogeneous::getLiq()
{
  return m_liq;
}

//****************************************************************************

int ModEulerHomogeneous::getVap()
{
  return m_vap;
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModEulerHomogeneous::selectScalar(Phase** phases, Mixture* mixture, Transport* transports, Variable nameVariable, int num) const
{
  switch (nameVariable) {
    case Variable::pressure:
      return mixture->getPressure();
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
    case Variable::density:
      return mixture->getDensity();
      break;
    case Variable::transport:
      return transports[num].getValue();
      break;
    case Variable::temperature:
      return mixture->getTemperature();
      break;
    default:
      Errors::errorMessage("nameVariable unknown in selectScalar"); return 0;
      break;
  }
}

//***********************************************************************