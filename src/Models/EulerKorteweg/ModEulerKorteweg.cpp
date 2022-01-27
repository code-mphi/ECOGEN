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

#include "ModEulerKorteweg.h"
#include "PhaseEulerKorteweg.h"

const std::string ModEulerKorteweg::NAME = "EULERKORTEWEG";

//****************************************************************************

ModEulerKorteweg::ModEulerKorteweg(const int& numbTransports) :
  Model(NAME, numbTransports)
{
  fluxBuff = new FluxEulerKorteweg();
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxEulerKorteweg());
  }
}

//***********************************************************************

ModEulerKorteweg::ModEulerKorteweg(const int& numbTransports, const double& alpha, const double &beta, const double &temperature, const double &kappa) :
  Model(NAME, numbTransports)
{
  fluxBuff = new FluxEulerKorteweg();
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxEulerKorteweg());
  }
  alphaEK = alpha;
  betaEK = beta;
  temperatureEK = temperature;
  kappaEK = kappa;
}

//***********************************************************************

ModEulerKorteweg::ModEulerKorteweg(const std::string& name, const int& numbTransports, const double& alpha, const double &beta, const double &temperature, const double &kappa) :
  Model(name, numbTransports)
{
  alphaEK = alpha;
  betaEK = beta;
  temperatureEK = temperature;
  kappaEK = kappa;
}

//****************************************************************************

ModEulerKorteweg::~ModEulerKorteweg()
{
  delete fluxBuff;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//****************************************************************************

void ModEulerKorteweg::allocateCons(Flux** cons)
{
  *cons = new FluxEulerKorteweg;
}

//***********************************************************************

void ModEulerKorteweg::allocatePhase(Phase** phase)
{
  *phase = new PhaseEulerKorteweg;
}

//***********************************************************************

void ModEulerKorteweg::allocateMixture(Mixture** mixture)
{
  *mixture = new MixEulerKorteweg;
}

//***********************************************************************

void ModEulerKorteweg::initializeAugmentedVariables(Cell* cell)
{
  Phase* phase(cell->getPhase(0));

  //Eta
  phase->setEta(phase->getDensity());
  //Omega
  Coord gradRho(cell->computeGradient(density, 0));
  double omega(0.); //This variable may need to be manually initialized
  omega = phase->getVelocity().norm() * gradRho.norm();
  phase->setOmega(omega);
  //VectorP
  phase->setVectorP(gradRho);
  //Pressure (not for NLS)
  if (phase->getEos() != nullptr) phase->setPressure(phase->getEos()->computePressure(phase->getDensity(), temperatureEK));
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModEulerKorteweg::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/) const
{
  //Get variables
  double rhoL, rhoR, omegaL, omegaR, etaL, etaR;
  double uL, uR, vL, vR;
  double vecPxL, vecPxR, vecPyL, vecPyR;
  double kappaL, kappaR;

  Phase* phaseLeft(0), *phaseRight(0);
  phaseLeft = cellLeft.getPhase(0);
  phaseRight = cellRight.getPhase(0);

  rhoL = phaseLeft->getDensity();
  omegaL = phaseLeft->getOmega();
  etaL = phaseLeft->getEta();
  uL = phaseLeft->getU(); vL = phaseLeft->getV(); // wL = phaseLeft->getW();
  vecPxL = phaseLeft->getVectorPX(); vecPyL = phaseLeft->getVectorPY(); // vecPzL = phaseLeft->getVectorPZ();
  kappaL = this->kappa(rhoL);

  rhoR = phaseRight->getDensity();
  omegaR = phaseRight->getOmega();
  etaR = phaseRight->getEta();
  uR = phaseRight->getU(); vR = phaseRight->getV(); // wR = phaseRight->getW();
  vecPxR = phaseRight->getVectorPX(); vecPyR = phaseRight->getVectorPY(); // vecPzR = phaseRight->getVectorPZ();
  kappaR = this->kappa(rhoR);

  //Left and right tensor P in momentum equation
  double piL = rhoL*rhoL*this->epsilonPrime(cellLeft, rhoL) + 0.5*(rhoL*this->kappaPrime(rhoL) - kappaL)*phaseLeft->getVectorP().squaredNorm() + etaL/alphaEK*(1. - etaL/rhoL);
  double tensorPxxL= piL + kappaL * vecPxL * vecPxL;
  double tensorPyxL= kappaL * vecPxL * vecPyL;
  // double tensorPzxL= kappaL * vecPxL * vecPzL;

  double piR = rhoR*rhoR*this->epsilonPrime(cellRight, rhoR) + 0.5*(rhoR*this->kappaPrime(rhoR) - kappaR)*phaseRight->getVectorP().squaredNorm() + etaR/alphaEK*(1. - etaR/rhoR);
  double tensorPxxR= piR + kappaR * vecPxR * vecPxR;
  double tensorPyxR= kappaR * vecPxR * vecPyR;
  // double tensorPzxR= kappaR * vecPxR * vecPzR;

  //Compute maximal wave speed using Davis approximation and eigenvalues of hyperbolic equations
  double maxWaveSpeed = this->computeMaxWaveSpeed(cellLeft, cellRight, rhoL, rhoR, uL, uR, etaL, etaR, vecPxL, vecPxR, vecPyL, vecPyR);
  if (maxWaveSpeed > 1.e-3) {
    dtMax = std::min(dtMax, dxLeft / maxWaveSpeed);
    dtMax = std::min(dtMax, dxRight / maxWaveSpeed);
  }

  //Compute Rusanov approximate solver
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_mass = 0.5 * (rhoL*uL + rhoR*uR - maxWaveSpeed*(rhoR - rhoL));
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqOmega = 0.5 * (rhoL*omegaL*uL - kappaL*vecPxL/betaEK + rhoR*omegaR*uR - kappaR*vecPxR/betaEK - maxWaveSpeed*(rhoR*omegaR - rhoL*omegaL));
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqEta = 0.5 * (rhoL*etaL*uL + rhoR*etaR*uR - maxWaveSpeed*(rhoR*etaR - rhoL*etaL));
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.setX(0.5 * (rhoL*uL*uL + tensorPxxL + rhoR*uR*uR + tensorPxxR - maxWaveSpeed*(rhoR*uR - rhoL*uL)));
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.setY(0.5 * (rhoL*uL*vL + tensorPyxL + rhoR*uR*vR + tensorPyxR - maxWaveSpeed*(rhoR*vR - rhoL*vL)));
  // static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.setZ(0.5 * (rhoL*uL*wL + tensorPzxL + rhoR*uR*wR + tensorPzxR - maxWaveSpeed*(rhoR*wR - rhoL*wL)));
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.setX(0.5 * (vecPxL*uL + vecPyL*vL - omegaL + vecPxR*uR + vecPyR*vR  - omegaR - maxWaveSpeed*(vecPxR - vecPxL)));
  // static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.setX(0.5 * (vecPxL*uL + vecPyL*vL + vecPzL*wL - omegaL + vecPxR*uR + vecPyR*vR + vecPzR*wR - omegaR - maxWaveSpeed*(vecPxR - vecPxL)));
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.setY(- 0.5 * maxWaveSpeed*(vecPyR - vecPyL));
  // static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.setZ(- 0.5 * maxWaveSpeed*(vecPzR - vecPzL));
}

void ModEulerKorteweg::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  Coord fluxProjete;
  fluxProjete.setX(normal.getX()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getX() + tangent.getX()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getY() + binormal.getX()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getZ());
  fluxProjete.setY(normal.getY()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getX() + tangent.getY()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getY() + binormal.getY()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getZ());
  fluxProjete.setZ(normal.getZ()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getX() + tangent.getZ()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getY() + binormal.getZ()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.getZ());
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_momentum.setXYZ(fluxProjete.getX(), fluxProjete.getY(), fluxProjete.getZ());

  fluxProjete.setX(normal.getX()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getX() + tangent.getX()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getY() + binormal.getX()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getZ());
  fluxProjete.setY(normal.getY()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getX() + tangent.getY()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getY() + binormal.getY()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getZ());
  fluxProjete.setZ(normal.getZ()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getX() + tangent.getZ()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getY() + binormal.getZ()*static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.getZ());
  static_cast<FluxEulerKorteweg*> (fluxBuff)->m_eqVectorP.setXYZ(fluxProjete.getX(), fluxProjete.getY(), fluxProjete.getZ());
}

//****************************************************************************
//********************* Methods specific to Euler-Korteweg *******************
//****************************************************************************

double ModEulerKorteweg::kappa(const double& /*density*/) const
{
  return kappaEK;
}

//****************************************************************************

double ModEulerKorteweg::kappaPrime(const double& /*density*/) const
{
  return 0.;
}

//****************************************************************************

double ModEulerKorteweg::kappaSecond(const double& /*density*/) const
{
  return 0.;
}

//****************************************************************************

double ModEulerKorteweg::epsilonPrime(Cell& cell, const double& density) const
{
  return cell.getPhase(0)->getEos()->dedrho(density, temperatureEK);
}

//****************************************************************************

double ModEulerKorteweg::epsilonSecond(Cell& cell, const double& density) const
{
  return cell.getPhase(0)->getEos()->dedrhoSecond(density, temperatureEK);
}

//****************************************************************************

double ModEulerKorteweg::computeMaxWaveSpeed(Cell& cellLeft, Cell& cellRight, const double& rhoL, const double& rhoR,
  const double& uL, const double& uR, const double& etaL, const double& etaR, const double& vecPxL, const double& vecPxR,
  const double& vecPyL, const double& vecPyR) const
{
  //Left
  double aL = 2.*this->epsilonPrime(cellLeft, rhoL) + 0.5*(vecPxL*vecPxL + vecPyL*vecPyL)*this->kappaSecond(rhoL) + rhoL*this->epsilonSecond(cellLeft, rhoL) + etaL*etaL/(alphaEK*rhoL*rhoL*rhoL);
  double bL = rhoL*aL + (vecPxL*vecPxL + vecPyL*vecPyL + 1./betaEK)*this->kappa(rhoL)/rhoL + 2.*vecPxL*vecPxL*this->kappaPrime(rhoL);
  double cL = (vecPyL*vecPyL + 1./betaEK) * (this->kappa(rhoL) * aL - 2.*vecPxL*vecPxL*std::pow(this->kappaPrime(rhoL), 2.));
  double soundSpeedMaxL = std::sqrt((bL + std::sqrt(bL*bL - 4.*cL)) / 2.);
  
  //Right
  double aR = 2.*this->epsilonPrime(cellRight, rhoR) + 0.5*(vecPxR*vecPxR + vecPyR*vecPyR)*this->kappaSecond(rhoR) + rhoR*this->epsilonSecond(cellRight, rhoR) + etaR*etaR/(alphaEK*rhoR*rhoR*rhoR);
  double bR = rhoR*aR + (vecPxR*vecPxR + vecPyR*vecPyR + 1./betaEK)*this->kappa(rhoR)/rhoR + 2.*vecPxR*vecPxR*this->kappaPrime(rhoR);
  double cR = (vecPyR*vecPyR + 1./betaEK) * (this->kappa(rhoR) * aR - 2.*vecPxR*vecPxR*std::pow(this->kappaPrime(rhoR), 2.));
  double soundSpeedMaxR = std::sqrt((bR + std::sqrt(bR*bR - 4.*cR)) / 2.);
  
  //Max wave speed
  double S = std::max(std::abs(uL), std::abs(uR)) + std::max(soundSpeedMaxL, soundSpeedMaxR);
  return S;
}

//****************************************************************************