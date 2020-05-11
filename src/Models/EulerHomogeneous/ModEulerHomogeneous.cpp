//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

//! \file      ModEulerHomogeneousous.cpp
//! \author    F. Petitpas, K. Schmidmayer, J. Caze
//! \version   1.1
//! \date      November 19 2019

#include <cmath>
#include <algorithm>
#include "ModEulerHomogeneous.h"
#include "PhaseEulerHomogeneous.h"

const std::string ModEulerHomogeneous::NAME = "EULERHOMOGENEOUS";

//****************************************************************************

ModEulerHomogeneous::ModEulerHomogeneous(const int &numberTransports, const int liquid, const int vapor) :
  m_liq(liquid), m_vap(vapor), Model(NAME, numberTransports)
{
  fluxBuff = new FluxEulerHomogeneous(); 
}

//****************************************************************************

ModEulerHomogeneous::~ModEulerHomogeneous()
{
  delete fluxBuff;
}

//****************************************************************************

void ModEulerHomogeneous::allocateCons(Flux **cons, const int &numberPhases)
{
  *cons = new FluxEulerHomogeneous(this);
}

//***********************************************************************

void ModEulerHomogeneous::allocatePhase(Phase **phase)
{
  *phase = new PhaseEulerHomogeneous;
}

//***********************************************************************

void ModEulerHomogeneous::allocateMixture(Mixture **mixture)
{
  *mixture = new MixEulerHomogeneous;
}

//***********************************************************************

void ModEulerHomogeneous::fulfillState(Phase **phases, Mixture *mixture, const int &numberPhases, Prim type)
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
  mixture->computeMixtureVariables(phases, numberPhases);
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModEulerHomogeneous::solveRiemannIntern(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const
{
  double sL, sR;
  
  //FP//TODO//look for sound speed
  double uL = cellLeft.getMixture()->getVelocity().getX(), vL = cellLeft.getMixture()->getVelocity().getY(), wL = cellLeft.getMixture()->getVelocity().getZ(), cL = cellLeft.getMixture()->getMixSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uR = cellRight.getMixture()->getVelocity().getX(), vR = cellRight.getMixture()->getVelocity().getY(), wR = cellRight.getMixture()->getVelocity().getZ(), cR = cellRight.getMixture()->getMixSoundSpeed(), pR = cellRight.getMixture()->getPressure(), rhoR = cellRight.getMixture()->getDensity();
  double EL = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm(), ER = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();

  //Davies
  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (std::fabs(sM)<1.e-8) sM = 0.;

  if (sL > 0.){
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_masse = rhoL*uL;
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setX(rhoL*uL*uL + pL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setY(rhoL*vL*uL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setZ(rhoL*wL*uL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = (rhoL*EL + pL)*uL;
  }
  else if (sR < 0.){
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_masse = rhoR*uR;
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setX(rhoR*uR*uR + pR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setY(rhoR*vR*uR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setZ(rhoR*wR*uR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = (rhoR*ER + pR)*uR;
  }

  ////1) Option HLL
  //else if (std::fabs(sR - sL)>1.e-3)
  //{
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_masse = (rhoR*uR*sL - rhoL*uL*sR + sL*sR*(rhoL - rhoR)) / (sL - sR);
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setX(((rhoR*uR*uR + pR)*sL - (rhoL*uL*uL + pL)*sR + sL*sR*(rhoL*uL - rhoR*uR)) / (sL - sR));
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setY((rhoR*uR*vR*sL - rhoL*uL*vL*sR + sL*sR*(rhoL*vL - rhoR*vR)) / (sL - sR));
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setZ((rhoR*uR*wR*sL - rhoL*uL*wL*sR + sL*sR*(rhoL*wL - rhoR*wR)) / (sL - sR));
  //  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = ((rhoR*ER + pR)*uR*sL - (rhoL*EL + pL)*uL*sR + sL*sR*(rhoL*EL - rhoR*ER)) / (sL - sR);
  //}

  //2) Option HLLC
  else if (sM >= 0.) {
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_masse = rhoStar*sM;
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM+pStar);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setY(rhoStar*sM*vL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setZ(rhoStar*sM*wL);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;
  }
  else {
    double pStar = mR*(sM - uR) + pR;
    double rhoStar = mR / (sR - sM);
    double Estar = ER + (sM - uR)*(sM + pR / mR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_masse = rhoStar*sM;
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setY(rhoStar*sM*vR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setZ(rhoStar*sM*wR);
    static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;
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
  Coord fluxProjete;
  fluxProjete.setX(normal.getX()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getX() + tangent.getX()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getY() + binormal.getX()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getZ());
  fluxProjete.setY(normal.getY()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getX() + tangent.getY()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getY() + binormal.getY()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getZ());
  fluxProjete.setZ(normal.getZ()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getX() + tangent.getZ()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getY() + binormal.getZ()*static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.getZ());
  static_cast<FluxEulerHomogeneous*> (fluxBuff)->m_qdm.setXYZ(fluxProjete.getX(), fluxProjete.getY(), fluxProjete.getZ());
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