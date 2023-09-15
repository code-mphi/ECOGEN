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

#include "ModEuler.h"

const std::string ModEuler::NAME = "EULER";

//****************************************************************************

ModEuler::ModEuler(const int& numbTransports) :
  Model(NAME, numbTransports)
{
  fluxBuff = new FluxEuler();
  fluxBuffMRF = new FluxEuler();
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxEuler());
  }
}

//****************************************************************************

ModEuler::~ModEuler()
{
  delete fluxBuff;
  delete fluxBuffMRF;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//****************************************************************************

void ModEuler::allocateCons(Flux** cons)
{
  *cons = new FluxEuler;
}

//***********************************************************************

void ModEuler::allocatePhase(Phase** phase)
{
  *phase = new PhaseEuler;
}

//***********************************************************************

void ModEuler::allocateMixture(Mixture** mixture)
{
  *mixture = new MixEuler;
}

//***********************************************************************

void ModEuler::allocatePhaseGradient(GradPhase** phase)
{
  *phase = new GradPhaseEuler;
}

//***********************************************************************

void ModEuler::allocateMixtureGradient(GradMixture** mixture)
{
  *mixture = new GradMixEuler;
}

//***********************************************************************

void ModEuler::fulfillState(Phase** phases, Mixture* /*mixture*/)
{
  phases[0]->extendedCalculusPhase(phases[0]->getVelocity());
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModEuler::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double> &boundData) const
{
  double cL, cR, sL, sR;
  double uL, uR, vL, vR, wL, wR, pL, pR, rhoL, rhoR, EL, ER;

  Phase* phaseLeft(0), *phaseRight(0);
  phaseLeft = cellLeft.getPhase(0);
  phaseRight = cellRight.getPhase(0);

  uL = phaseLeft->getU(); vL = phaseLeft->getV(); wL = phaseLeft->getW();
  pL = phaseLeft->getPressure();
  rhoL = phaseLeft->getDensity();
  cL = phaseLeft->getSoundSpeed();
  EL = phaseLeft->getTotalEnergy();

  uR = phaseRight->getU(); vR = phaseRight->getV(); wR = phaseRight->getW();
  pR = phaseRight->getPressure();
  rhoR = phaseRight->getDensity();
  cR = phaseRight->getSoundSpeed();
  ER = phaseRight->getTotalEnergy();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL, uR, cR);
  }

  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (std::fabs(sM)<1.e-8) sM = 0.;

  if (sL > 0.){
    static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoL*uL;
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoL*uL*uL + pL);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoL*vL*uL);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoL*wL*uL);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoL*EL + pL)*uL;

    // Boundary data for output
    boundData[VarBoundary::p] = pL;
    boundData[VarBoundary::rho] = rhoL;
    boundData[VarBoundary::velU] = uL;
    boundData[VarBoundary::velV] = vL;
    boundData[VarBoundary::velW] = wL;
  }
  else if (sR < 0.){
    static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoR*uR;
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoR*uR*uR + pR);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoR*vR*uR);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoR*wR*uR);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoR*ER + pR)*uR;

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
  //  static_cast<FluxEuler*> (fluxBuff)->m_mass = (rhoR*uR*sL - rhoL*uL*sR + sL*sR*(rhoL - rhoR)) / (sL - sR);
  //  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(((rhoR*uR*uR + pR)*sL - (rhoL*uL*uL + pL)*sR + sL*sR*(rhoL*uL - rhoR*uR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY((rhoR*uR*vR*sL - rhoL*uL*vL*sR + sL*sR*(rhoL*vL - rhoR*vR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ((rhoR*uR*wR*sL - rhoL*uL*wL*sR + sL*sR*(rhoL*wL - rhoR*wR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_energ = ((rhoR*ER + pR)*uR*sL - (rhoL*EL + pL)*uL*sR + sL*sR*(rhoL*EL - rhoR*ER)) / (sL - sR);
  //}

  //2) Option HLLC
  else if (sM >= 0.) {
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar*sM;
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM+pStar);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar*sM*vL);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar*sM*wL);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;

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
    static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar*sM;
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar*sM*vR);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar*sM*wR);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;

    // Boundary data for output
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vR;
    boundData[VarBoundary::velW] = wR;
  }

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
//*** Half Riemann solver for MRF interface between static/rotating region ***
//****************************************************************************

void ModEuler::solveRiemannInternMRF(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, const Coord& omega, const Coord& normal, const Coord& tangent, const Coord& binormal, const Coord& position) const
{
  double cL, cR, sL, sR;
  double uL, uR, vL, vR, wL, wR, pL, pR, rhoL, rhoR, EL, ER;
  Coord velocityStar;

  Phase* phaseLeft(0), *phaseRight(0);
  phaseLeft = cellLeft.getPhase(0);
  phaseRight = cellRight.getPhase(0);

  uL = phaseLeft->getU(); vL = phaseLeft->getV(); wL = phaseLeft->getW();
  pL = phaseLeft->getPressure();
  rhoL = phaseLeft->getDensity();
  cL = phaseLeft->getSoundSpeed();
  EL = phaseLeft->getTotalEnergy();

  uR = phaseRight->getU(); vR = phaseRight->getV(); wR = phaseRight->getW();
  pR = phaseRight->getPressure();
  rhoR = phaseRight->getDensity();
  cR = phaseRight->getSoundSpeed();
  ER = phaseRight->getTotalEnergy();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL, uR, cR);
  }

  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (std::fabs(sM)<1.e-8) sM = 0.;

  if (sL > 0.){
    // Flux for static region
    static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoL*uL;
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoL*uL*uL + pL);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoL*vL*uL);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoL*wL*uL);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoL*EL + pL)*uL;
    
    // Flux for rotating region
    velocityStar.setXYZ(uL, vL, wL);
    velocityStar.buildRelativeVelForRiemannMRF(omega, normal, tangent, binormal, position);
    static_cast<FluxEuler*> (fluxBuffMRF)->m_mass = rhoL*velocityStar.getX();
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setX(rhoL*velocityStar.getX()*velocityStar.getX() + pL);
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setY(rhoL*velocityStar.getY()*velocityStar.getX());
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setZ(rhoL*velocityStar.getZ()*velocityStar.getX());
    EL = phaseLeft->getEnergy() + 0.5 * velocityStar.squaredNorm();
    static_cast<FluxEuler*> (fluxBuffMRF)->m_energ = (rhoL*EL + pL)*velocityStar.getX();
  }
  else if (sR < 0.){
    // Flux for static region
    static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoR*uR;
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoR*uR*uR + pR);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoR*vR*uR);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoR*wR*uR);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoR*ER + pR)*uR;
    
    // Flux for rotating region
    velocityStar.setXYZ(uR, vR, wR);
    velocityStar.buildRelativeVelForRiemannMRF(omega, normal, tangent, binormal, position);
    static_cast<FluxEuler*> (fluxBuffMRF)->m_mass = rhoR*velocityStar.getX();
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setX(rhoR*velocityStar.getX()*velocityStar.getX() + pR);
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setY(rhoR*velocityStar.getY()*velocityStar.getX());
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setZ(rhoR*velocityStar.getZ()*velocityStar.getX());
    ER = phaseRight->getEnergy() + 0.5 * velocityStar.squaredNorm();
    static_cast<FluxEuler*> (fluxBuffMRF)->m_energ = (rhoR*ER + pR)*velocityStar.getX();
  }

  ////1) Option HLL
  //else if (std::fabs(sR - sL)>1.e-3)
  //{
  //  static_cast<FluxEuler*> (fluxBuff)->m_mass = (rhoR*uR*sL - rhoL*uL*sR + sL*sR*(rhoL - rhoR)) / (sL - sR);
  //  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(((rhoR*uR*uR + pR)*sL - (rhoL*uL*uL + pL)*sR + sL*sR*(rhoL*uL - rhoR*uR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY((rhoR*uR*vR*sL - rhoL*uL*vL*sR + sL*sR*(rhoL*vL - rhoR*vR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ((rhoR*uR*wR*sL - rhoL*uL*wL*sR + sL*sR*(rhoL*wL - rhoR*wR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_energ = ((rhoR*ER + pR)*uR*sL - (rhoL*EL + pL)*uL*sR + sL*sR*(rhoL*EL - rhoR*ER)) / (sL - sR);
  //}

  //2) Option HLLC
  else if (sM >= 0.) {
    // Flux for static region
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar*sM;
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM+pStar);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar*sM*vL);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar*sM*wL);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;
    
    // Flux for rotating region
    velocityStar.setXYZ(sM, vL, wL);
    velocityStar.buildRelativeVelForRiemannMRF(omega, normal, tangent, binormal, position);
    static_cast<FluxEuler*> (fluxBuffMRF)->m_mass = rhoStar*velocityStar.getX();
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setX(rhoStar*velocityStar.getX()*velocityStar.getX() + pStar);
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setY(rhoStar*velocityStar.getY()*velocityStar.getX());
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setZ(rhoStar*velocityStar.getZ()*velocityStar.getX());
    Estar = Estar - 0.5 * (sM*sM+vL*vL+wL*wL);
    Estar += 0.5*velocityStar.squaredNorm(); 
    static_cast<FluxEuler*> (fluxBuffMRF)->m_energ = (rhoStar*Estar + pStar)*velocityStar.getX();
  }
  else {
    // Flux for static region
    double pStar = mR*(sM - uR) + pR;
    double rhoStar = mR / (sR - sM);
    double Estar = ER + (sM - uR)*(sM + pR / mR);
    static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar*sM;
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar*sM*vR);
    static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar*sM*wR);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;

    // Flux for rotating region
    velocityStar.setXYZ(sM, vR, wR);
    velocityStar.buildRelativeVelForRiemannMRF(omega, normal, tangent, binormal, position);
    static_cast<FluxEuler*> (fluxBuffMRF)->m_mass = rhoStar*velocityStar.getX();
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setX(rhoStar*velocityStar.getX()*velocityStar.getX() + pStar);
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setY(rhoStar*velocityStar.getY()*velocityStar.getX());
    static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setZ(rhoStar*velocityStar.getZ()*velocityStar.getX());
    Estar = Estar - 0.5 * (sM*sM+vR*vR+wR*wR);
    Estar += 0.5*velocityStar.squaredNorm(); 
    static_cast<FluxEuler*> (fluxBuffMRF)->m_energ = (rhoStar*Estar + pStar)*velocityStar.getX();
  }

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions** **************
//****************************************************************************

void ModEuler::solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double> &boundData) const
{
  double cL, sL;
  double uL, pL, rhoL;
  double pStar(0.);

  Phase* phaseLeft(0);
  phaseLeft = cellLeft.getPhase(0);

  uL = phaseLeft->getU();
  pL = phaseLeft->getPressure();
  rhoL = phaseLeft->getDensity();
  cL = phaseLeft->getSoundSpeed();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  sL = std::min(uL - cL, -uL - cL);
  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  pStar = rhoL*uL*(uL - sL) + pL;

  static_cast<FluxEuler*> (fluxBuff)->m_mass = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = 0.;

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = 0.;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = 0.;
  boundData[VarBoundary::velU] = 0.;
  boundData[VarBoundary::velV] = 0.;
  boundData[VarBoundary::velW] = 0.;
}

//****************************************************************************

void ModEuler::solveRiemannInletTank(Cell& cellLeft, const double& dxLeft, double& dtMax, const double* /*ak0*/, const double* rhok0, const double& p0, const double& /*T0*/, std::vector<double> &boundData) const
{
  Eos* eos;

  double cL, sL, zL;
  double uL, pL, rhoL, vL, wL;
  double uStar(0.), rhoStar(0.), pStar(0.), eStar(0.), vStar(0.), wStar(0.);

  Phase* phaseLeft(0);
  phaseLeft = cellLeft.getPhase(0);

  eos = phaseLeft->getEos();
  uL = phaseLeft->getU();
  vL = phaseLeft->getV();
  wL = phaseLeft->getW();
  pL = phaseLeft->getPressure();
  rhoL = phaseLeft->getDensity();
  cL = phaseLeft->getSoundSpeed();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL*cL;

  //Left wave velocity estimation using pStar = p0
  //----------------------------------------------
  pStar = p0;
  double v(0.), vmv0, mL, u;
  v = 1./eos->computeDensityIsentropic(pL, rhoL, pStar); 
  //v = 1. / eos->computeDensityHugoniot(pL, rhoL, pStar); //Other possibility
  vmv0 = v - 1. / rhoL;
  if (std::fabs(vmv0) > 1e-10) { mL = sqrt((pL - p0) / vmv0); }
  else { mL = zL; }
  sL = uL - mL / rhoL;
  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  u = uL + mL*vmv0;

  //Pathologic cases
  //----------------
  if (sL >= 0.) { //supersonic outflow => left state solution
    uStar = uL;
    pStar = pL;
    rhoStar = rhoL;
    vStar = vL;
    wStar = wL;
  }
  else if (u >= -1e-3) { //subsonic outflow => star left state solution
    uStar = u;
    pStar = p0;
    rhoStar = 1. / v;
    vStar = vL;
    wStar = wL;
  }
  //Tank case
  //---------
  else { //tank inflow => star right state solution
    //Total enthalpy in tank state
    double H0(0.);
    double v0 = 0.;
    H0 = eos->computeTotalEnthalpy(rhok0[0], p0, v0);

    //ITERATIVE PROCESS FOR PRESSURE DETERMINATION 
    //--------------------------------------------
    int iteration(0);
    double p(0.5*p0);
    double f(0.), df(1.);
    double dv, h, dh, drho;
    double uStarR(0.), duStarR(0.), uStarL(0.), duStarL(0.);
    do {
      p -= f / df; iteration++;
      if (iteration > 50) Errors::errorMessage("solveRiemannInletTank not converged in modEuler");
      //Physical pressure ?
      eos->verifyAndModifyPressure(p);
      if (p > p0) { p = p0 - 1e-6; }

      //Tank rekations in the right (H=cte et s=cste)
      rhoStar = eos->computeDensityIsentropic(p0, rhok0[0], p);
      h = eos->computeEnthalpyIsentropic(p0, rhok0[0], p, &dh);
      uStarR = -sqrt(2.*(H0 - h));
      duStarR = -dh / uStarR; ;

      //Isentropic relations on the left //FP//DEV// chocs a mettre
      double dmL;
      v = 1.0 / eos->computeDensityIsentropic(pL, rhoL, p, &drho);
      //v = 1.0 / eos->computeDensityHugoniot(pL, rhoL, p, &drho); //Other possibility
      dv = - v*v*drho;
      vmv0 = v - 1. / rhoL;
      if (std::fabs(vmv0) > 1e-10) {
        mL = sqrt((pL - p) / vmv0);
        dmL = 0.5*(-vmv0 + (p - pL)*dv) / (vmv0*vmv0) / mL;
      }
      else { //if limit density overpassed under shock => acoustic relations instead
        mL = zL;
        dmL = 0.;
      }
      sL = uL - mL / rhoL;
      // For low-Mach (for general purpose machRefMin set to 1)
      if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
      uStarL = uL + mL*vmv0;
      duStarL = dmL*vmv0 + mL*dv;

      //solved function
      f = uStarR - uStarL;
      df = duStarR - duStarL;

    } while (std::fabs(f)>1e-3);

    pStar = p;
    uStar = 0.5*(uStarL + uStarR);
    vStar = 0.;
    wStar = 0.;
  }

  eStar = eos->computeEnergy(rhoStar, pStar);

  static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar*uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar*uStar*vStar);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar*uStar*wStar);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vStar*vStar + wStar*wStar)) + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = vStar;
  boundData[VarBoundary::velW] = wStar;
}

//****************************************************************************

void ModEuler::solveRiemannInletInjStagState(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, const double* /*ak0*/, const double* rhok0, const double* pk0, std::vector<double> &boundData) const
{
  Eos* eos;
  double H0, u0;

  double cL, sL, zL;
  double uL, pL, rhoL, vL, wL;
  double uStar(0.), rhoStar(0.), pStar(0.), eStar(0.);

  Phase* phaseLeft(0);
  phaseLeft = cellLeft.getPhase(0);

  eos = phaseLeft->getEos();
  uL = phaseLeft->getU();
  vL = phaseLeft->getV();
  wL = phaseLeft->getW();
  pL = phaseLeft->getPressure();
  rhoL = phaseLeft->getDensity();
  cL = phaseLeft->getSoundSpeed();
  
  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL*cL;
  sL = uL - cL;

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  u0 = m0 / rhok0[0];
  H0 = eos->computeTotalEnthalpy(rhok0[0], pk0[0], u0);
  double c0 = eos->computeSoundSpeed(rhok0[0], pk0[0]);

  //Null Mass flow
  //--------------
  if (fabs(u0) < 1.e-6) {
    rhoStar = rhoL;
    pStar = pL;
    uStar = 0.;
  }
  //Supersonic inflow
  //-----------------
  else if (u0 < -c0) {
    rhoStar = rhok0[0];
    pStar = pk0[0];
    uStar = u0;
  }
  //Subsonic inflow
  //---------------
  else {
    int iteration(0);
    pStar = pL;
    double f(0.), df(1.);
    double u, du, v, dv, h;
    do {
      pStar -= f / df; iteration++;
      if (iteration > 50) Errors::errorMessage("solveRiemannInflow not converged in modEuler");
      //physical pressure ?
      eos->verifyAndModifyPressure(pStar);
      //Acoustic relation in the left (can be modified by shock relations)
      u = uL + (pL - pStar) / zL;
      if (u >= -1e-6) u = -1e-6;
      du = -1. / zL;
      f = u / u0; df = du / u0;
      //Compute from m0, H0 on the right
      h = H0 - 0.5 * u * u;
      v = eos->vfpfh(pStar, h);
      dv = eos->dvdpch(pStar, h) - eos->dvdhcp(pStar) * u * du;
      f -= v * rhok0[0];
      df -= dv * rhok0[0];

    } while (std::fabs(f) > 1e-8 && iteration <= 50);
    uStar = u;
    rhoStar = m0 / uStar;
  }
  eStar = eos->computeEnergy(rhoStar, pStar);

  ////IG or SG version only   =>   exact solver
  ////-----------------------------------------
  //double* dataEos;
  //double a, b, c, delta, u1, u2, gammaTemp;
  //eos->sendInfo(dataEos);
  //gammaTemp = (dataEos[0] - 1)* m0 / dataEos[0];
  //a = 0.5*gammaTemp - zL;
  //b = pL + zL*uL;
  //c = -gammaTemp*H0;
  //delta = b*b - 4 * a*c;
  //u1 = (-b - sqrt(delta)) / (2 * a);
  //u2 = (-b + sqrt(delta)) / (2 * a);
  //uStar = std::min(u1, u2);
  //pStar = pL + zL*uL - zL*uStar;
  //rhoStar = m0 / uStar;
  //eStar = eos->computeEnergy(rhoStar, pStar);

  static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar*uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar*uStar*vL);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar*uStar*wL);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vL*vL + wL*wL)) + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = vL;
  boundData[VarBoundary::velW] = wL;
}

//****************************************************************************

void ModEuler::solveRiemannInletInjTemp(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, const double* Tk0, const double* /*ak0*/, std::vector<double> &boundData) const
{
  Eos* eos;
  Phase* leftPhase(0);
  double uL, vL, wL, pL, rhoL;
  double cL, sL, zL;

  leftPhase = cellLeft.getPhase(0);
  eos = leftPhase->getEos();
  uL = leftPhase->getU();
  vL = leftPhase->getV();
  wL = leftPhase->getW();
  pL = leftPhase->getPressure();
  rhoL = leftPhase->getDensity();
  cL = leftPhase->getSoundSpeed();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL * cL;

  int it(0);
  double p(pL), u(0.), du(0.), f(0.), df(1.);
  double mL(zL), dmL(0.), vSmvL(0.);
  double rhoStarL(0.), drhoStarL(0.), vStarL(0.), dvStarL(0.);
  double eStar(0.), rhoStar(0.), pStar(0.), uStar(0.);

  // ITERATIVE PROCESS FOR PRESSURE DETERMINATION
  // --------------------------------------------
  do {
    p -= f / df;
    it++;
    if (it > 50) Errors::errorMessage("solveRiemannInletInjTemp not converged in modEuler");
    eos->verifyAndModifyPressure(p);

    // Left intermediate state
    // rhoStarL = eos->computeDensityIsentropic(pL, rhoL, p, &drhoStarL);
    rhoStarL = eos->computeDensityHugoniot(pL, rhoL, p, &drhoStarL);
    vStarL = 1. / rhoStarL;
    dvStarL = - drhoStarL / (rhoStarL * rhoStarL);
    vSmvL = vStarL - 1. / rhoL;
    if (std::fabs(vSmvL) > 1.e-10) { // Rankine-Hugoniot
      mL = std::sqrt((pL - p) / vSmvL);
      dmL = 0.5 * (- vSmvL + (p - pL) * dvStarL) / (vSmvL * vSmvL) / mL;
    }
    else { // Acoustic
      mL = zL;
      dmL = 0.; 
    }
    sL = uL - mL / rhoL;
    
    u = uL + mL * vSmvL;
    du = dmL * vSmvL + mL * dvStarL;

    f = m0 * vStarL - u;
    df = m0 * dvStarL - du;

  } while (std::fabs(f) > 1.e-5 && it <= 50);

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  pStar = p;
  rhoStar = eos->computeDensity(pStar, Tk0[0]);
  eStar = eos->computeEnergy(rhoStar, pStar);
  uStar = u;

  static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar * uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar * uStar * uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar * uStar * vL);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar * uStar * wL);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar * (eStar + 0.5 * (uStar * uStar + vL * vL + wL * wL)) + pStar) * uStar;

  // Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = vL;
  boundData[VarBoundary::velW] = wL;
}

//****************************************************************************

void ModEuler::solveRiemannOutletPressure(Cell& cellLeft, const double& dxLeft, double& dtMax, const double p0, std::vector<double> &boundData) const
{
  double cL, sL, zL;
  double uL, pL, rhoL, vL, wL;
  double uStar(0.), rhoStar(0.), pStar(p0), eStar(0.);

  Phase* phaseLeft(0);
  phaseLeft = cellLeft.getPhase(0);
  
  uL = phaseLeft->getU();
  vL = phaseLeft->getV();
  wL = phaseLeft->getW();
  pL = phaseLeft->getPressure();
  rhoL = phaseLeft->getDensity();
  cL = phaseLeft->getSoundSpeed();
  
  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL*cL;

  // Perturbed state
  double vStar(0.), vSmvL(0.), mL(zL), sM(0.);
  //rhoStar = TB->eos[0]->computeDensityIsentropic(pL, rhoL, pStar); // Isentropic relations if needed
  rhoStar = TB->eos[0]->computeDensityHugoniot(pL, rhoL, pStar);
  vStar = 1. / std::max(rhoStar, 1e-10);
  vSmvL = vStar - 1. / rhoL;

  // Jump accross left wave is described by RH or acoustic relations (if specific volume jump is small)
  if (std::abs(vSmvL) > 1e-10) { mL = sqrt((pL - pStar) / vSmvL); }
  else { mL = zL; }
  sL = uL - mL / rhoL;
  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  sM = uL + mL * vSmvL;

  // Pathologic case I: sL>0
  if (sL >= 0.) { //Supersonic outflow => Left state solution
    uStar = uL;
    pStar = pL;
    rhoStar = rhoL;
  }
  // Pathologic case II: inflow conditions, we temporarly keep the specific mass
  else if (sM < 0) {
    uStar = sM;
    rhoStar = rhoL;
  }
  // Imposed pressure outflow
  else {
    uStar = sM;
  }

  // Flux completion
  eStar = TB->eos[0]->computeEnergy(rhoStar, pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar*uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar*uStar*vL);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar*uStar*wL);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vL*vL + wL*wL)) + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = vL;
  boundData[VarBoundary::velW] = wL;
}

//****************************************************************************

void ModEuler::solveRiemannOutletMassflow(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, std::vector<double>& boundData) const
{
  Eos* eos;
  Phase* leftPhase(0);
  double uL, vL, wL, pL, rhoL;
  double cL, sL, zL;

  leftPhase = cellLeft.getPhase(0);
  eos = leftPhase->getEos();
  uL = leftPhase->getU();
  vL = leftPhase->getV();
  wL = leftPhase->getW();
  pL = leftPhase->getPressure();
  rhoL = leftPhase->getDensity();
  cL = leftPhase->getSoundSpeed();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL * cL;

  int it(0);
  double p(0.8*pL), u(0.), du(0.), f(0.), df(1.);
  double mL(zL), dmL(0.), vSmvL(0.);
  double rhoStarL(0.), drhoStarL(0.), vStarL(0.), dvStarL(0.);
  double eStar(0.), rhoStar(0.), pStar(0.), uStar(0.);

  // ITERATIVE PROCESS FOR PRESSURE DETERMINATION
  // --------------------------------------------
  do {
    p -= f / df;
    it++;
    if (it > 50) { 
      warnings.push_back(Errors("solveRiemannOutletMassflow not converged in ModEuler", __FILE__, __LINE__));
    }
    // Check physical pressure
    eos->verifyAndModifyPressure(p);

    // Left intermediate state
    // rhoStarL = eos->computeDensityIsentropic(pL, rhoL, p, &drhoStarL);
    rhoStarL = eos->computeDensityHugoniot(pL, rhoL, p, &drhoStarL);
    vStarL = 1. / rhoStarL;
    dvStarL = - drhoStarL / (rhoStarL * rhoStarL);
    vSmvL = vStarL - 1. / rhoL;
    if (std::fabs(vSmvL) > 1.e-10) { // Rankine-Hugoniot
      mL = std::sqrt((pL - p) / vSmvL);
      dmL = 0.5 * (- vSmvL + (p - pL) * dvStarL) / (vSmvL * vSmvL) / mL;
    }
    else { // Acoustic
      mL = zL;
      dmL = 0.; 
    }

    sL = uL - mL / rhoL;

    u = uL + mL * vSmvL;
    du = dmL * vSmvL + mL * dvStarL;

    f = m0 - rhoStarL * u;
    df = - drhoStarL * u - rhoStarL * du;
  } while (std::fabs(f) > 1.e-3 && it <= 50);

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  // Solution state
  if (sL >= 0.) { // Pathologic case: supersonic outflow
    uStar = uL;
    rhoStar = rhoL;
    pStar = pL;
  }
  else if (u < 0.) { // Pathologic case: Back flow (we keep the mass temporarly)
    uStar = u;
    rhoStar = rhoL;
    pStar = p;
  }
  else { // Subsonic outflow
    pStar = p;
    rhoStar = rhoStarL;
    eStar = eos->computeEnergy(rhoStar, pStar);
    uStar = u;
  }

  // Flux completion
  static_cast<FluxEuler*> (fluxBuff)->m_mass = rhoStar * uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(rhoStar * uStar * uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(rhoStar * uStar * vL);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(rhoStar * uStar * wL);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar * (eStar + 0.5 * (uStar * uStar + vL * vL + wL * wL)) + pStar) * uStar;

  // Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = vL;
  boundData[VarBoundary::velW] = wL;
}

//****************************************************************************

void ModEuler::solveRiemannNullFlux() const
{
  static_cast<FluxEuler*> (fluxBuff)->m_mass = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setX(0.);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = 0.;
}

//****************************************************************************
//********************** Transport Riemann solvers ***************************
//****************************************************************************

void ModEuler::solveRiemannTransportIntern(Cell& cellLeft, Cell& cellRight)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemann(cellLeft.getTransport(k).getValue(), cellRight.getTransport(k).getValue(), static_cast<FluxEuler*> (fluxBuff)->m_sM);
  }
}

//****************************************************************************

void ModEuler::solveRiemannTransportWall()
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannWall();
  }
}

//****************************************************************************

void ModEuler::solveRiemannTransportPiston(Cell& cellLeft, double uPiston)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannPiston(cellLeft.getTransport(k).getValue(), uPiston);
  }
}


//****************************************************************************

void ModEuler::solveRiemannTransportInletTank(Cell& cellLeft, double* valueTransports)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannInletTank(cellLeft.getTransport(k).getValue(), static_cast<FluxEuler*> (fluxBuff)->m_sM, valueTransports[k]);
  }
}

//****************************************************************************

void ModEuler::solveRiemannTransportInletInjStagState(Cell& cellLeft, double* valueTransports)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannInletInjStagState(cellLeft.getTransport(k).getValue(), static_cast<FluxEuler*> (fluxBuff)->m_sM, valueTransports[k]);
  }
}

//****************************************************************************

void ModEuler::solveRiemannTransportOutletPressure(Cell& cellLeft, double* valueTransports)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannOutletPressure(cellLeft.getTransport(k).getValue(), static_cast<FluxEuler*> (fluxBuff)->m_sM, valueTransports[k]);
  }
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModEuler::selectScalar(Phase** phases, Mixture* /*mixture*/, Transport* transports, Variable nameVariable, int num) const
{
  switch (nameVariable) {
    case Variable::pressure:
      return phases[0]->getPressure();
      break;
    case Variable::density:
      return phases[0]->getDensity();
      break;
    case Variable::velocityU:
      return phases[0]->getU();
      break;
    case Variable::velocityV:
      return phases[0]->getV();
      break;
    case Variable::velocityW:
      return phases[0]->getW();
      break;
    case Variable::velocityMag:
      return phases[0]->getVelocity().norm();
      break;
    case Variable::transport:
      return transports[num].getValue();
      break;
    case Variable::temperature:
      return phases[0]->getTemperature();
      break;
    default:
      Errors::errorMessage("nameVariable unknown in selectScalar"); return 0;
      break;
  }
}

//****************************************************************************

const double& ModEuler::getSM()
{
  return static_cast<FluxEuler*> (fluxBuff)->m_sM;
}

//****************************************************************************
//***************************** Others methods *******************************
//****************************************************************************

void ModEuler::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  static_cast<FluxEuler*> (fluxBuff)->m_momentum.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************

void ModEuler::lowMachSoundSpeed(double& machRef, const double& uL, double& cL, const double& uR, double& cR) const
{
  // Low-Mach preconditioning
  // ------------------------
  // See eq. (24) of LeMartelot, S., Nkonga, B., & Saurel, R. (2013). Liquid and liquid-gas flows at all speeds. 
  // Journal of Computational Physics, 255, 53-82.

  // --- Mref ---
  // double machRefMin = 0.1; 
  // cL = (1. - machRefMin*machRefMin)*uL + sqrt((machRefMin*machRefMin - 1.)*(machRefMin*machRefMin - 1.) * uL*uL + 4.* machRefMin*machRefMin*cL*cL);
  // cL *= 0.5;
  // cR = (machRefMin*machRefMin - 1.)*uR + sqrt((machRefMin*machRefMin - 1.)*(machRefMin*machRefMin - 1.) * uR*uR + 4.* machRefMin*machRefMin*cR*cR);
  // cR *= 0.5;
  // machRef = machRefMin; // For CFL criteria
  
  // --- Mlocal ---
  double machLimitComp(0.3);
  double machRefL(0.);
  if (std::fabs(uL)/cL >= machLimitComp) machRefL = 1.;
  else if (std::fabs(uL)/cL > m_machRefMin) machRefL = std::fabs(uL)/cL;
  else machRefL = m_machRefMin;
  // Right
  double machRefR(0.);
  if (std::fabs(uR)/cR >= machLimitComp) machRefR = 1.;
  else if (std::fabs(uR)/cR > m_machRefMin) machRefR = std::fabs(uR)/cR;
  else machRefR = m_machRefMin;
  machRef = std::max(machRefL, machRefR);
  //Caution: Keep right before left for boundary condition Riemann solver
  cR = (machRef*machRef - 1.)*uR + sqrt((machRef*machRef - 1.)*(machRef*machRef - 1.) * uR*uR + 4.* machRef*machRef*cR*cR);
  cR *= 0.5;
  cL = (1. - machRef*machRef)*uL + sqrt((machRef*machRef - 1.)*(machRef*machRef - 1.) * uL*uL + 4.* machRef*machRef*cL*cL);
  cL *= 0.5;
}

//****************************************************************************

void ModEuler::addNonConsMrfFlux(Phase** phases)
{
  static_cast<FluxEuler*> (fluxBuffMRF)->addNonConsMrfFlux(phases);
}

//****************************************************************************

void ModEuler::reverseProjectionMrfFlux(const Coord normal, const Coord tangent, const Coord binormal) const
{
  Coord projectedFlux;
  projectedFlux.setX(normal.getX()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getX() + tangent.getX()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getY() + binormal.getX()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getZ());
  projectedFlux.setY(normal.getY()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getX() + tangent.getY()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getY() + binormal.getY()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getZ());
  projectedFlux.setZ(normal.getZ()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getX() + tangent.getZ()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getY() + binormal.getZ()*static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.getZ());
  static_cast<FluxEuler*> (fluxBuffMRF)->m_momentum.setXYZ(projectedFlux.getX(), projectedFlux.getY(), projectedFlux.getZ());
}

//****************************************************************************