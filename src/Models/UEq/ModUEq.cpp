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

#include "ModUEq.h"

const std::string ModUEq::NAME = "VELOCITYEQ";

//***********************************************************************

ModUEq::ModUEq(const int& numbTransports, const int& numbPhases) : Model(NAME, numbTransports)
{
  fluxBuff = new FluxUEq(numbPhases);
  fluxBuffMRF = new FluxUEq(numbPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxUEq(numbPhases));
  }
}

//***********************************************************************

ModUEq::ModUEq(const std::string& name, const int& numbTransports) : Model(name, numbTransports){}

//***********************************************************************

ModUEq::~ModUEq()
{
  delete fluxBuff;
  delete fluxBuffMRF;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//***********************************************************************

void ModUEq::allocateCons(Flux** cons)
{
  *cons = new FluxUEq(numberPhases);
}

//***********************************************************************

void ModUEq::allocatePhase(Phase** phase)
{
  *phase = new PhaseUEq;
}

//***********************************************************************

void ModUEq::allocateMixture(Mixture** mixture)
{
  *mixture = new MixUEq;
}

//***********************************************************************

void ModUEq::allocatePhaseGradient(GradPhase** phase)
{
  *phase = new GradPhaseUEq;
}

//***********************************************************************

void ModUEq::allocateMixtureGradient(GradMixture** mixture)
{
  *mixture = new GradMixUEq;
}

//***********************************************************************

void ModUEq::fulfillState(Phase** phases, Mixture* mixture)
{
  //Complete phases state
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  //Complete mixture variables using phases variable
  mixture->computeMixtureVariables(phases);
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModUEq::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double> &boundData) const
{
  Phase* vecPhase;
  double sL, sR;
  double pStar(0.), rhoStar(0.), EStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getFrozenSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uR = cellRight.getMixture()->getVelocity().getX(), cR = cellRight.getMixture()->getFrozenSoundSpeed(), pR = cellRight.getMixture()->getPressure(), rhoR = cellRight.getMixture()->getDensity();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL, uR, cR);
  }

  //Davies
  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxRight / std::fabs(sR));

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
      double energy = vecPhase->getEnergy();
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = alpha*density*uL;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha*density*energy*uL;
    }
    double vitY = cellLeft.getMixture()->getVelocity().getY(); double vitZ = cellLeft.getMixture()->getVelocity().getZ();
    double totalEnergy = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoL*uL*uL + pL);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoL*vitY*uL);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoL*vitZ*uL);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoL*totalEnergy + pL)*uL;
    static_cast<FluxUEq*> (fluxBuff)->m_uStar = uL;

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
      double energy = vecPhase->getEnergy();
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = alpha*density*uR;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha*density*energy*uR;
    }
    double vitY = cellRight.getMixture()->getVelocity().getY(); double vitZ = cellRight.getMixture()->getVelocity().getZ();
    double totalEnergy = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoR*uR*uR + pR);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoR*vitY*uR);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoR*vitZ*uR);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoR*totalEnergy + pR)*uR;
    static_cast<FluxUEq*> (fluxBuff)->m_uStar = uR;

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
    double totalEnergy = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
    rhoStar = mL / (sL - sM);
    EStar = totalEnergy + (sM - uL)*(sM + pL / mL);
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
      // TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
      if(!m_lowMach) TB->ekStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]);
      else { TB->ekStar[k] = vecPhase->getEnergy(); }
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = alpha* TB->rhokStar[k] * sM;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha* TB->rhokStar[k] * TB->ekStar[k] * sM;
    }
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoStar*vitY*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoStar*vitZ*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*sM;
    static_cast<FluxUEq*> (fluxBuff)->m_uStar = sM;

    // Boundary data for output
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vitY;
    boundData[VarBoundary::velW] = vitZ;
  }
  else{
    //Compute right solution state
    double vitY = cellRight.getMixture()->getVelocity().getY(); double vitZ = cellRight.getMixture()->getVelocity().getZ();
    double totalEnergy = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();
    rhoStar = mR / (sR - sM);
    EStar = totalEnergy + (sM - uR)*(sM + pR / mR);
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
      // TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
      if(!m_lowMach) TB->ekStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]);
      else { TB->ekStar[k] = vecPhase->getEnergy(); }
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = alpha* TB->rhokStar[k] * sM;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha* TB->rhokStar[k] * TB->ekStar[k] * sM;
    }
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoStar*vitY*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoStar*vitZ*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*sM;
    static_cast<FluxUEq*> (fluxBuff)->m_uStar = sM;

    // Boundary data for output
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vitY;
    boundData[VarBoundary::velW] = vitZ;
  }

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
//*** Half Riemann solver for MRF interface between static/rotating region ***
//****************************************************************************

void ModUEq::solveRiemannInternMRF(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, const Coord& omega, const Coord& normal, const Coord& tangent, const Coord& binormal, const Coord& position) const
{
  Phase* vecPhase;
  double sL, sR;
  double pStar(0.), rhoStar(0.), EStar(0.);
  Coord velocityStar;

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getFrozenSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uR = cellRight.getMixture()->getVelocity().getX(), cR = cellRight.getMixture()->getFrozenSoundSpeed(), pR = cellRight.getMixture()->getPressure(), rhoR = cellRight.getMixture()->getDensity();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL, uR, cR);
  }

  //Davies
  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR)), mkL, mkR;
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (std::fabs(sM)<1.e-8) sM = 0.;

  //Solution sampling
  if (sL >= 0.){
    double vitY = cellLeft.getMixture()->getVelocity().getY(); double vitZ = cellLeft.getMixture()->getVelocity().getZ();
    velocityStar.setXYZ(uL, vitY, vitZ);
    velocityStar.buildRelativeVelForRiemannMRF(omega, normal, tangent, binormal, position);
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellLeft.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double density = vecPhase->getDensity();
      double energy = vecPhase->getEnergy();
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = alpha*density*uL;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha*density*energy*uL;
      static_cast<FluxUEq*> (fluxBuffMRF)->m_alpha[k] = alpha*velocityStar.getX(); //FP//Should be related to sM (choice on the other components?)
      static_cast<FluxUEq*> (fluxBuffMRF)->m_mass[k] = alpha*density*velocityStar.getX();
      static_cast<FluxUEq*> (fluxBuffMRF)->m_energ[k] = alpha*density*energy*velocityStar.getX();
    }
    double totalEnergy = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoL*uL*uL + pL);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoL*vitY*uL);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoL*vitZ*uL);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoL*totalEnergy + pL)*uL;
    static_cast<FluxUEq*> (fluxBuff)->m_uStar = uL;
    totalEnergy = cellLeft.getMixture()->getEnergy() + 0.5*velocityStar.squaredNorm();
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setX(rhoL*velocityStar.getX()*velocityStar.getX() + pL);
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setY(rhoL*velocityStar.getY()*velocityStar.getX());
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setZ(rhoL*velocityStar.getZ()*velocityStar.getX());
    static_cast<FluxUEq*> (fluxBuffMRF)->m_energMixture = (rhoL*totalEnergy + pL)*velocityStar.getX();
    static_cast<FluxUEq*> (fluxBuffMRF)->m_uStar = velocityStar.getX();
  }
  else if (sR <= 0.){
    double vitY = cellRight.getMixture()->getVelocity().getY(); double vitZ = cellRight.getMixture()->getVelocity().getZ();
    velocityStar.setXYZ(uR, vitY, vitZ);
    velocityStar.buildRelativeVelForRiemannMRF(omega, normal, tangent, binormal, position);
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellRight.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double density = vecPhase->getDensity();
      double energy = vecPhase->getEnergy();
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = alpha*density*uR;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha*density*energy*uR;
      static_cast<FluxUEq*> (fluxBuffMRF)->m_alpha[k] = alpha*velocityStar.getX(); //FP//Should be related to sM (choice on the other components?)
      static_cast<FluxUEq*> (fluxBuffMRF)->m_mass[k] = alpha*density*velocityStar.getX();
      static_cast<FluxUEq*> (fluxBuffMRF)->m_energ[k] = alpha*density*energy*velocityStar.getX();
    }
    double totalEnergy = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoR*uR*uR + pR);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoR*vitY*uR);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoR*vitZ*uR);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoR*totalEnergy + pR)*uR;
    static_cast<FluxUEq*> (fluxBuff)->m_uStar = uR;
    totalEnergy = cellRight.getMixture()->getEnergy() + 0.5*velocityStar.squaredNorm();
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setX(rhoR*velocityStar.getX()*velocityStar.getX() + pR);
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setY(rhoR*velocityStar.getY()*velocityStar.getX());
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setZ(rhoR*velocityStar.getZ()*velocityStar.getX());
    static_cast<FluxUEq*> (fluxBuffMRF)->m_energMixture = (rhoR*totalEnergy + pR)*velocityStar.getX();
    static_cast<FluxUEq*> (fluxBuffMRF)->m_uStar = velocityStar.getX();    
  }
  else if (sM >= 0.){
    //Compute left solution state
    double vitY = cellLeft.getMixture()->getVelocity().getY(); double vitZ = cellLeft.getMixture()->getVelocity().getZ();
    velocityStar.setXYZ(sM, vitY, vitZ);
    velocityStar.buildRelativeVelForRiemannMRF(omega, normal, tangent, binormal, position);
    double totalEnergy = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
    rhoStar = mL / (sL - sM);
    EStar = totalEnergy + (sM - uL)*(sM + pL / mL);
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
      // TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
      if(!m_lowMach) TB->ekStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]);
      else { TB->ekStar[k] = vecPhase->getEnergy(); }
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = alpha* TB->rhokStar[k] * sM;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha* TB->rhokStar[k] * TB->ekStar[k] * sM;
      static_cast<FluxUEq*> (fluxBuffMRF)->m_alpha[k] = alpha*velocityStar.getX();
      static_cast<FluxUEq*> (fluxBuffMRF)->m_mass[k] = alpha*TB->rhokStar[k]*velocityStar.getX();
      static_cast<FluxUEq*> (fluxBuffMRF)->m_energ[k] = alpha*TB->rhokStar[k]*TB->ekStar[k]*velocityStar.getX();
    }
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoStar*vitY*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoStar*vitZ*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*sM;
    static_cast<FluxUEq*> (fluxBuff)->m_uStar = sM;
    totalEnergy = EStar - 0.5*(sM*sM+vitY*vitY+vitZ*vitZ);
    totalEnergy += 0.5*velocityStar.squaredNorm();
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setX(rhoStar*velocityStar.getX()*velocityStar.getX() + pStar);
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setY(rhoStar*velocityStar.getY()*velocityStar.getX());
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setZ(rhoStar*velocityStar.getZ()*velocityStar.getX());
    static_cast<FluxUEq*> (fluxBuffMRF)->m_energMixture = (rhoStar*totalEnergy + pStar)*velocityStar.getX();
    static_cast<FluxUEq*> (fluxBuffMRF)->m_uStar = velocityStar.getX();
  }
  else{
    //Compute right solution state
    double vitY = cellRight.getMixture()->getVelocity().getY(); double vitZ = cellRight.getMixture()->getVelocity().getZ();
    velocityStar.setXYZ(sM, vitY, vitZ);
    velocityStar.buildRelativeVelForRiemannMRF(omega, normal, tangent, binormal, position);
    double totalEnergy = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();
    rhoStar = mR / (sR - sM);
    EStar = totalEnergy + (sM - uR)*(sM + pR / mR);
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
      // TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
      if(!m_lowMach) TB->ekStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]);
      else { TB->ekStar[k] = vecPhase->getEnergy(); }
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = alpha* TB->rhokStar[k] * sM;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha* TB->rhokStar[k] * TB->ekStar[k] * sM;
      static_cast<FluxUEq*> (fluxBuffMRF)->m_alpha[k] = alpha*velocityStar.getX();
      static_cast<FluxUEq*> (fluxBuffMRF)->m_mass[k] = alpha*TB->rhokStar[k]*velocityStar.getX();
      static_cast<FluxUEq*> (fluxBuffMRF)->m_energ[k] = alpha*TB->rhokStar[k]*TB->ekStar[k]*velocityStar.getX();
    }
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoStar*vitY*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoStar*vitZ*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*sM;
    static_cast<FluxUEq*> (fluxBuff)->m_uStar = sM;
    totalEnergy = EStar - 0.5*(sM*sM+vitY*vitY+vitZ*vitZ);
    totalEnergy += 0.5*velocityStar.squaredNorm();
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setX(rhoStar*velocityStar.getX()*velocityStar.getX() + pStar);
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setY(rhoStar*velocityStar.getY()*velocityStar.getX());
    static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setZ(rhoStar*velocityStar.getZ()*velocityStar.getX());
    static_cast<FluxUEq*> (fluxBuffMRF)->m_energMixture = (rhoStar*totalEnergy + pStar)*velocityStar.getX();
    static_cast<FluxUEq*> (fluxBuffMRF)->m_uStar = velocityStar.getX();
  }

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions ****************
//****************************************************************************

void ModUEq::solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double> &boundData) const
{
  double sL;
  double pStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getFrozenSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  sL = std::min(uL - cL, -uL - cL);
  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  pStar = rhoL*(uL - sL)*uL + pL;

  for (int k = 0; k < numberPhases; k++)
  {
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = 0.;
  static_cast<FluxUEq*> (fluxBuff)->m_uStar = 0.;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = 0.;
  boundData[VarBoundary::velU] = 0.;
  boundData[VarBoundary::velV] = 0.;
  boundData[VarBoundary::velW] = 0.;
}

//****************************************************************************

void ModUEq::solveRiemannInletTank(Cell& cellLeft, const double& dxLeft, double& dtMax, const double* ak0, const double* rhok0, const double& p0, const double& /*T0*/, std::vector<double> &boundData) const
{
  double tabp[50], tabf[50];
  double sL, zL, sM, vmv0, mL;
  double pStar(0.), uStar(0.), rhoStar(0.), vStar(0.), uyStar(0.), uzStar(0.);
  Phase* vecPhase;

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getWoodSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uyL = cellLeft.getMixture()->getVelocity().getY(), uzL = cellLeft.getMixture()->getVelocity().getZ();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL*cL;

  //1) Left wave velocity estimation using pStar = p0
  //-------------------------------------------------
  pStar = p0; vStar = 0.;
  for (int k = 0; k < numberPhases; k++) {
    vecPhase = cellLeft.getPhase(k);
    //TB->rhokStar[k] = TB->eos[k]->computeDensityIsentropic(vecPhase->getPressure(), vecPhase->getDensity(), pStar); //other possiblity
    TB->rhokStar[k] = TB->eos[k]->computeDensityHugoniot(vecPhase->getPressure(), vecPhase->getDensity(), pStar);
    TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
    vStar += vecPhase->getAlpha()*vecPhase->getDensity() / rhoL / std::max(TB->rhokStar[k], epsilonAlphaNull);
  }
  vmv0 = vStar - 1. / rhoL;
  if (std::fabs(vmv0) > 1e-10) { mL = sqrt((pL - pStar) / vmv0); }
  else { mL = zL; }
  sL = uL - mL / rhoL;
  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  sM = uL + mL*vmv0;

  //2) Check for pathologic cases
  //-----------------------------
  if (sL >= 0.) { //supersonic outflow => left state solution
    uStar = uL;
    pStar = pL;
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellLeft.getPhase(k);
      TB->rhokStar[k] = vecPhase->getDensity();
      TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
      TB->YkStar[k] = vecPhase->getAlpha()*vecPhase->getDensity() / rhoL;
    }
    rhoStar = rhoL;
    uyStar = uyL;
    uzStar = uzL;
  }
  else if (sM >= -1e-3) { //subsonic outflow => star left state solution
    uStar = sM;
    pStar = p0;  //approximation
    for (int k = 0; k < numberPhases; k++) {
      // TB->rhokStar[k] unchanged : see 1)
      TB->YkStar[k] = cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() / rhoL;
    }
    rhoStar = 1. / vStar;
    uyStar = uyL;
    uzStar = uzL;
  }

  //3) Tank
  //-------
  else { //tank inflow => star right state solution
    //Total enthalpy in tank state
    double H0(0.);
    Coord u0(0.);
    double rho0 = cellLeft.getMixture()->computeDensity(ak0, rhok0);
    for (int k = 0;k < numberPhases;k++) {
      TB->Yk0[k] = ak0[k] * rhok0[k] / rho0;
      H0 += TB->Yk0[k] * TB->eos[k]->computeTotalEnthalpy(rhok0[k], p0, u0.norm());  //default zero velocity in tank
    }
    //ITERATIVE PROCESS FOR PRESSURE DETERMINATION 
    //--------------------------------------------
    int iteration(0);
    double p(0.5*p0);
    double f(0.), df(1.);
    double hk, dhk, rhok, drhok, dmL, YkL;
    double uStarR(0.), duStarR(0.), uStarL(0.), duStarL(0.);
    double vStarL(0.), dvStarL(0.);
    do {
      p -= f / df; iteration++;
      //Physical pressure ?
      for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(p); }
      if (p > p0) { p = p0 - 1e-6; }
      tabp[iteration - 1] = p; tabf[iteration - 1] = f;
      if (iteration > 50) {
        for (int i = 0; i < 50; i++) { std::cout << tabp[i] << " " << tabf[i] << std::endl; }
        Errors::errorMessage("solveRiemannInletTank not converged in ModUEq");
      }
      //R) Tank rekations in the right (H=cte et sk=cste)
      uStarR = H0; duStarR = 0.;
      for (int k = 0; k < numberPhases; k++) {
        TB->rhokStar[k] = TB->eos[k]->computeDensityIsentropic(p0, rhok0[k], p);
        TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
        hk = TB->eos[k]->computeEnthalpyIsentropic(p0, rhok0[k], p, &dhk);
        uStarR -= TB->Yk0[k] * hk;
        duStarR -= TB->Yk0[k] * dhk;
      }
      uStarR = -sqrt(2.*uStarR);
      duStarR = duStarR / uStarR;
      //L) Left relations sk=cste (could be R-H if needed)
      vStarL = 0.; dvStarL = 0.;
      for (int k = 0; k < numberPhases; k++) {
        vecPhase = cellLeft.getPhase(k);
        rhok = TB->eos[k]->computeDensityIsentropic(vecPhase->getPressure(), vecPhase->getDensity(), p, &drhok); //other possiblity
        //rhok = TB->eos[k]->computeDensityHugoniot(vecPhase->getPressure(), vecPhase->getDensity(), p, &drhok);
        TB->eos[k]->verifyAndCorrectDensityMax(rhok);
        YkL = vecPhase->getAlpha()*vecPhase->getDensity() / rhoL;
        vStarL += YkL / std::max(rhok, epsilonAlphaNull);
        dvStarL -= YkL / std::max((rhok * rhok), epsilonAlphaNull) * drhok;
      }
      vmv0 = vStarL - 1. / rhoL;
      if (std::fabs(vmv0) > 1e-10) {
        mL = sqrt((pL - p) / vmv0);
        dmL = 0.5*(-vmv0 + (p - pL)*dvStarL) / (vmv0*vmv0) / mL;
      }
      else {
        mL = zL;
        dmL = 0.;
      }
      sL = uL - mL / rhoL;
      // For low-Mach (for general purpose machRefMin set to 1)
      if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
      uStarL = uL + mL*vmv0;
      duStarL = dmL*vmv0 + mL*dvStarL;
      //solved function
      f = uStarR - uStarL;
      df = duStarR - duStarL;
    } while (std::fabs(f)>1e-2); //End iterative loop
    pStar = p;
    uStar = 0.5*(uStarL + uStarR);
    rhoStar = 0.;
    for (int k = 0; k < numberPhases; k++) { 
      TB->YkStar[k] = TB->Yk0[k];
      rhoStar += TB->YkStar[k] / std::max(TB->rhokStar[k], epsilonAlphaNull);
    }
    rhoStar = 1. / rhoStar;
    uyStar = 0.;
    uzStar = 0.;
    sM = uStar;

  } //End tank case

  //4) Flux completion
  //------------------
  double akStar;
  double EStar(0.5*(uStar*uStar + uyStar*uyStar + uzStar*uzStar)), ek;
  for (int k = 0; k < numberPhases; k++) {
    if(!m_lowMach) ek = TB->eos[k]->computeEnergy(TB->rhokStar[k], pStar); 
    else { ek = cellLeft.getPhase(k)->getEnergy(); }
    EStar += TB->YkStar[k] * ek;
    akStar = TB->YkStar[k] * rhoStar / std::max(TB->rhokStar[k], epsilonAlphaNull);
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = akStar * sM;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = akStar * uStar * TB->rhokStar[k];
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->m_mass[k] * ek;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(rhoStar*uStar*uyStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(rhoStar*uStar*uzStar);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = sM;
  static_cast<FluxUEq*> (fluxBuff)->m_uStar = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = uyStar;
  boundData[VarBoundary::velW] = uzStar;  
}

//****************************************************************************

void ModUEq::solveRiemannInletInjStagState(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, const double* ak0, const double* rhok0, const double* pk0, std::vector<double> &boundData) const
{
  double sL, zL, sM;
  double pStar(0.), uStar(0.), rhoStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getWoodSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  // double vL = cellLeft.getMixture()->getVelocity().getY(), wL = cellLeft.getMixture()->getVelocity().getZ();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  //Compute total enthalpy of injected fluid and speed of sound
  double rho0 = cellLeft.getMixture()->computeDensity(ak0, rhok0);
  double u0 = m0 / rho0;
  double c0(0.), p0(0.);
  for (int k = 0;k < numberPhases; k++) {
    p0 += ak0[k] * pk0[k];
    TB->Hk0[k] = TB->eos[k]->computeTotalEnthalpy(rhok0[k], pk0[k], u0);
    TB->Yk0[k] = ak0[k] * rhok0[k] / rho0;
    double ck = cellLeft.getPhase(k)->getEos()->computeSoundSpeed(rhok0[k], pk0[k]);
    c0 += ak0[k]/ std::max((rhok0[k] * ck * ck), epsilonAlphaNull);
  }
  c0 = 1./ sqrt(rho0 * c0);

  //Estimates for acoustic wave sL
  sL = uL - cL;
  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  zL = rhoL*cL;

  //Null Mass flow
  //--------------
  if (fabs(u0) < 1.e-6) {
    uStar = 0.;
    sM = 0.;
    pStar = pL;
    rhoStar = rhoL;
    for (int k = 0; k < numberPhases; k++) {
      TB->vkStar[k] = 1. / cellLeft.getPhase(k)->getDensity();
    }
  }
  //Supersonic inflow
  //-----------------
  else if (u0 < -c0) {
    uStar = u0;
    sM = u0;
    pStar = p0;
    rhoStar = rho0;
    for (int k = 0; k < numberPhases; k++) {
      TB->vkStar[k] = 1. / rhok0[k];
    }
  }
  else {
    //Subsonic inflow
    //---------------
    int iteration(0);
    pStar = pL;
    double f(0.), df(1.);
    double u, du, hk;
    double vStar(0.), dvStar(0.), dvk(0.);

    do {
      pStar -= f / df; iteration++;
      if (iteration > 50) { exit(0); Errors::errorMessage("solveRiemannInletInjStagState not converged in ModUEq"); }
      //Physical pressure ?
      for (int k = 0; k < numberPhases; k++) {
        TB->eos[k]->verifyAndModifyPressure(pStar);
      }
      if (pStar - pL - zL*uL < 0.) pStar = (1. + 1.e-6)*(pL + zL * uL);
      //Left acoustic relations
      u = uL + (pL - pStar) / zL;
      du = -1. / zL;
      
      f = u; df = du;

      //Compute from m0, Hk0, Yk0 on the right
      vStar = 0.; dvStar = 0.;
      for (int k = 0; k < numberPhases; k++) {
        hk = TB->Hk0[k] - 0.5 * u * u;
        TB->vkStar[k] = TB->eos[k]->vfpfh(pStar, hk);
        dvk = TB->eos[k]->dvdpch(pStar, hk) - TB->eos[k]->dvdhcp(pStar) * u * du;
        vStar += TB->Yk0[k] * TB->vkStar[k]; // Option 1
        dvStar += TB->Yk0[k] * dvk;
        // f -= TB->Yk0[k] * TB->vkStar[k]; // Option 2
        // df -= TB->Yk0[k] * dvk;
      }
      f -= m0 * vStar; // Option 1
      df -= m0 * dvStar;
    } while (std::fabs(f) > 1e-8 && iteration <= 50);
    uStar = u;
    sM = u;
    rhoStar = m0 / uStar;
  }

  //Flux completion
  double akStar;
  double Estar(0.5*(uStar * uStar)), ek, rhok;
  for (int k = 0; k<numberPhases; k++) {
    rhok = 1. / TB->vkStar[k];
    if(!m_lowMach) ek = TB->eos[k]->computeEnergy(rhok, pStar);
    else { ek = cellLeft.getPhase(k)->getEnergy(); }
    Estar += TB->Yk0[k] * ek;
    akStar = TB->Yk0[k] * TB->vkStar[k] * rhoStar;
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = akStar * sM;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = akStar * uStar * rhok;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->m_mass[k] * ek;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(uStar * uStar * rhoStar + pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (Estar * rhoStar + pStar)* uStar;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = sM;
  static_cast<FluxUEq*> (fluxBuff)->m_uStar = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = 0.;
  boundData[VarBoundary::velW] = 0.;
}

//****************************************************************************

void ModUEq::solveRiemannInletInjTemp(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, const double* Tk0, const double* ak0, std::vector<double> &boundData) const
{
  double sL, zL, sM;

  double uL = cellLeft.getMixture()->getVelocity().getX();
  double pL = cellLeft.getMixture()->getPressure();
  double rhoL = cellLeft.getMixture()->getDensity();
  double cL = cellLeft.getMixture()->getWoodSoundSpeed(); 

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL*cL;

  int it(0);
  double pStar(pL), uStar(0.), rhoStar(0.), f(0.), df(1.);
  double u(0.), du(0.);
  double rhok(0.), drhok(0.), vStarL(0.), dvStarL(0.), vSmvL(0.), mL(zL), dmL(0.), YkL(0.);

  Phase* vecPhase;

  do {
    pStar -= f / df;
    it++;
    if (it > 50) { Errors::errorMessage("solveRiemannInletInjTemp not converged in ModUEq"); }
    // Check physical pressure
    for (int k = 0; k < numberPhases; k++) {
      TB->eos[k]->verifyAndModifyPressure(pStar);
    }
    // Left acoustic relations
    vStarL = 0.; dvStarL = 0.;
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellLeft.getPhase(k);
      rhok = TB->eos[k]->computeDensityIsentropic(vecPhase->getPressure(), vecPhase->getDensity(), pStar, &drhok); //other possiblity
      // rhok = TB->eos[k]->computeDensityHugoniot(vecPhase->getPressure(), vecPhase->getDensity(), pStar, &drhok);
      YkL = vecPhase->getAlpha()*vecPhase->getDensity() / rhoL;
      vStarL += YkL / std::max(rhok, epsilonAlphaNull);
      dvStarL -= YkL / std::max((rhok *rhok), epsilonAlphaNull) * drhok;
    }

    vSmvL = vStarL - 1. / rhoL;
    if (std::fabs(vSmvL) > 1.e-10) { 
      mL = sqrt((pL - pStar) / vSmvL); 
      dmL = 0.5*(-vSmvL + (p - pL)*dvStarL) / (vSmvL*vSmvL) / mL;
    }
    else { 
      mL = zL; 
      dmL = 0.;
    }
    mL = zL; // Force acoustic relations
    dmL = 0.;
    
    sL = uL - mL / rhoL;
    
    u = uL + mL * vSmvL;
    du = dmL * vSmvL + mL * dvStarL;

    f = m0 * vStarL - u;
    df = m0 * dvStarL - du;
    
  } while (std::fabs(f) > 1e-8 && it <= 50);
  uStar = u;
  sM = u;
  rhoStar = m0 / uStar;

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  //Flux completion
  double Estar(0.5*(uStar * uStar)), ekStar;
  for (int k = 0; k<numberPhases; k++) {
    TB->rhokStar[k] = cellLeft.getPhase(k)->getEos()->computeDensity(pStar, Tk0[k]);
    if(!m_lowMach) ekStar = TB->eos[k]->computeEnergy(TB->rhokStar[k], pStar); 
    else { ekStar = cellLeft.getPhase(k)->getEnergy(); }
    TB->YkStar[k] = ak0[k] * TB->rhokStar[k] / rhoStar;
    Estar +=  TB->YkStar[k] * ekStar;
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = ak0[k] * sM;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = ak0[k] * uStar * TB->rhokStar[k];
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->m_mass[k] * ekStar;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(uStar * uStar * rhoStar + pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (Estar * rhoStar + pStar)* uStar;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = sM;
  static_cast<FluxUEq*> (fluxBuff)->m_uStar = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = 0.;
  boundData[VarBoundary::velW] = 0.;
}

//****************************************************************************

void ModUEq::solveRiemannOutletPressure(Cell& cellLeft, const double& dxLeft, double& dtMax, const double p0, std::vector<double> &boundData) const
{
  double sL, sM, zL;
  double pStar(p0), EStar(0.), vStar(0.), uStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getWoodSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uyL = cellLeft.getMixture()->getVelocity().getY(), uzL = cellLeft.getMixture()->getVelocity().getZ();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL*cL;

  //Left wave : isentropic wave assumption
  //--------------------------------------
  double vSmvL, mL(zL);
  for (int k = 0; k < numberPhases; k++) {
    //TB->rhokStar[k] = TB->eos[k]->computeDensityIsentropic(pL, vecPhase->getDensity(), pStar); //other possiblity
    TB->rhokStar[k] = TB->eos[k]->computeDensityHugoniot(pL, cellLeft.getPhase(k)->getDensity(), pStar);
    TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
    vStar += cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() /rhoL / std::max(TB->rhokStar[k], epsilonAlphaNull);
  }
  vSmvL = vStar - 1. / rhoL;
  if (std::fabs(vSmvL) > 1e-10) { mL = sqrt((pL - pStar) / vSmvL); }
  sL = uL - mL / rhoL;
  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  sM = uL + mL*vSmvL;

  //Pathologic case sL>0
  if (sL >= 0.) { //Supersonic outflow => Left state solution
    uStar = uL;
    pStar = pL;
    for (int k = 0; k < numberPhases; k++) { TB->rhokStar[k] = cellLeft.getPhase(k)->getDensity(); }
    vStar = 1. / rhoL;
  }
  else if (sM < 0) { //Inflow conditions : the outflow assumption is not adapted
    uStar = sM;
    for (int k = 0; k < numberPhases; k++) { TB->rhokStar[k] = cellLeft.getPhase(k)->getDensity(); }
    vStar = 1. / rhoL;
  }
  else { //imposed pressure outflow OK
    uStar = sM;
  }

  //Flux completion
  double ekStar, YkL, akStar;
  EStar = 0.5*(uStar*uStar + uyL * uyL + uzL * uzL);
  for (int k = 0; k < numberPhases; k++) {
    YkL = cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() / rhoL;
    akStar = YkL / std::max(TB->rhokStar[k], epsilonAlphaNull) / vStar;

    if(!m_lowMach) ekStar = TB->eos[k]->computeEnergy(TB->rhokStar[k], pStar);
    else { ekStar = cellLeft.getPhase(k)->getEnergy(); }

    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = akStar * sM;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = akStar * uStar * TB->rhokStar[k];
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->m_mass[k] * ekStar;

    EStar += YkL * ekStar;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(uStar*uStar / vStar + pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(uStar*uyL / vStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(uStar*uzL / vStar);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (EStar / vStar + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = sM;
  static_cast<FluxUEq*> (fluxBuff)->m_uStar = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = 1./vStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = uyL;
  boundData[VarBoundary::velW] = uzL;
}

//****************************************************************************

void ModUEq::solveRiemannOutletMassflow(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, std::vector<double>& boundData) const
{
  Phase* vecPhase;
  double sL, zL;
  double uL = cellLeft.getMixture()->getVelocity().getX();
  double uyL = cellLeft.getMixture()->getVelocity().getY(); 
  double uzL = cellLeft.getMixture()->getVelocity().getZ();
  double cL = cellLeft.getMixture()->getWoodSoundSpeed();
  double pL = cellLeft.getMixture()->getPressure();
  double rhoL = cellLeft.getMixture()->getDensity();

  // Low-Mach preconditioning
  double machRefMin(1.); // Default value without low-Mach preco.
  if(m_lowMach){
    lowMachSoundSpeed(machRefMin, uL, cL);
  }

  zL = rhoL*cL;

  // ITERATIVE PROCESS FOR PRESSURE DETERMINATION
  // --------------------------------------------
  int it(0);
  double p(pL), u(0.), du(0.), f(0.), df(1.);
  double vStar(0.), dvStar(0.), vSmvL(0.), mL(zL), dmL(0.);
  double drhokStar(0.), YkL(0.);

  do {
    p -= f / df;
    it++;
    if (it > 50) { 
      warnings.push_back(Errors("solveRiemannOutletMassflow not converged in ModUEq", __FILE__, __LINE__));
    }
    // Check physical pressure
    for (int k = 0; k < numberPhases; k++) {
      TB->eos[k]->verifyAndModifyPressure(p);
    }
    vStar = 0.; dvStar = 0.;
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellLeft.getPhase(k);
      //TB->rhokStar[k] = TB->eos[k]->computeDensityIsentropic(vecPhase->getPressure(), vecPhase->getDensity(), pStar, &drhokStar); //other possiblity
      TB->rhokStar[k] = TB->eos[k]->computeDensityHugoniot(vecPhase->getPressure(), vecPhase->getDensity(), p, &drhokStar);
      TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
      YkL = vecPhase->getAlpha()*vecPhase->getDensity() / rhoL;
      vStar += YkL / std::max(TB->rhokStar[k], epsilonAlphaNull);
      dvStar -= YkL / std::max(TB->rhokStar[k] * TB->rhokStar[k], epsilonAlphaNull) * drhokStar;
    }
    vSmvL = vStar - 1. / rhoL;
    if (std::fabs(vSmvL) > 1.e-10) { // Rankine-Hugoniot
      mL = sqrt((pL - p) / vSmvL); 
      dmL = 0.5 * (- vSmvL + (p - pL) * dvStar) / (vSmvL * vSmvL) / mL;
    }
    else { // Acoustic
      mL = zL; 
      dmL = 0.;
    }
    
    sL = uL - mL / rhoL;

    u = uL + mL * vSmvL;
    du = dmL * vSmvL + mL * dvStar;

    f = m0 - u / vStar;
    df = - (du * vStar - u * dvStar) / std::max(vStar * vStar, epsilonAlphaNull);
  } while (std::fabs(f) > 1.e-3 && it <= 50);

  // For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));

  // Solution state
  double uStar, pStar;
  if (sL >= 0.) { // Pathologic case: supersonic outflow
    uStar = uL;
    vStar = 1. / rhoL;
    pStar = pL;
  }
  else if (u < 0.) { // Pathologic case: Back flow (we keep the mass temporarly)
    uStar = u;
    vStar = 1. / rhoL;
    pStar = p;
  }
  else { // Subsonic outflow
    pStar = p;
    uStar = u;
  }

  // Flux completion
  double ekStar, EStar, akStar;
  EStar = 0.5*(uStar*uStar + uyL * uyL + uzL * uzL);
  for (int k = 0; k < numberPhases; k++) {
    YkL = cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() / rhoL;
    akStar = YkL / std::max(TB->rhokStar[k], epsilonAlphaNull) / vStar;

    if(!m_lowMach) ekStar = TB->eos[k]->computeEnergy(TB->rhokStar[k], pStar);
    else { ekStar = cellLeft.getPhase(k)->getEnergy(); }

    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = akStar * u;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = akStar * uStar * TB->rhokStar[k];
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->m_mass[k] * ekStar;

    EStar += YkL * ekStar;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(uStar*uStar / vStar + pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(uStar*uyL / vStar);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(uStar*uzL / vStar);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (EStar / vStar + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = u;
  static_cast<FluxUEq*> (fluxBuff)->m_uStar = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = 1./vStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = uyL;
  boundData[VarBoundary::velW] = uzL;
}

//****************************************************************************

void ModUEq::solveRiemannNullFlux() const
{
  for (int k = 0; k < numberPhases; k++) {
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;
  static_cast<FluxUEq*> (fluxBuff)->m_sM = 0.;
  static_cast<FluxUEq*> (fluxBuff)->m_uStar = 0.;
}

//****************************************************************************
//********************** Transport Riemann solvers ***************************
//****************************************************************************

void ModUEq::solveRiemannTransportIntern(Cell& cellLeft, Cell& cellRight)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemann(cellLeft.getTransport(k).getValue(), cellRight.getTransport(k).getValue(), static_cast<FluxUEq*> (fluxBuff)->m_sM);
  }
}

//****************************************************************************

void ModUEq::solveRiemannTransportWall()
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannWall();
  }
}

//****************************************************************************

void ModUEq::solveRiemannTransportPiston(Cell& cellLeft, double uPiston)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannPiston(cellLeft.getTransport(k).getValue(), uPiston);
  }
}


//****************************************************************************

void ModUEq::solveRiemannTransportInletTank(Cell& cellLeft, double* valueTransports)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannInletTank(cellLeft.getTransport(k).getValue(), static_cast<FluxUEq*> (fluxBuff)->m_sM, valueTransports[k]);
  }
}

//****************************************************************************

void ModUEq::solveRiemannTransportInletInjStagState(Cell& cellLeft, double* valueTransports)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannInletInjStagState(cellLeft.getTransport(k).getValue(), static_cast<FluxUEq*> (fluxBuff)->m_sM, valueTransports[k]);
  }
}

//****************************************************************************

void ModUEq::solveRiemannTransportOutletPressure(Cell& cellLeft, double* valueTransports)
{
  for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannOutletPressure(cellLeft.getTransport(k).getValue(), static_cast<FluxUEq*> (fluxBuff)->m_sM, valueTransports[k]);
  }
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModUEq::selectScalar(Phase** phases, Mixture* mixture, Transport* transports, Variable nameVariable, int num) const
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

//***********************************************************************

const double& ModUEq::getSM()
{
  return static_cast<FluxUEq*> (fluxBuff)->m_sM;
}

//****************************************************************************
//***************************** Others methods *******************************
//****************************************************************************

void ModUEq::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************

void ModUEq::lowMachSoundSpeed(double& machRef, const double& uL, double& cL, const double& uR, double& cR) const
{
  // Low-Mach preconditioning
  // ------------------------
  // See eq. (24) of LeMartelot, S., Nkonga, B., & Saurel, R. (2013). Liquid and liquid-gas flows at all speeds. 
  // Journal of Computational Physics, 255, 53-82.

  // --- Mref ---
  // double machRefMin = 0.1; 
  // cR = (machRefMin*machRefMin - 1.)*uR + sqrt((machRefMin*machRefMin - 1.)*(machRefMin*machRefMin - 1.) * uR*uR + 4.* machRefMin*machRefMin*cR*cR);
  // cR *= 0.5;
  // cL = (1. - machRefMin*machRefMin)*uL + sqrt((machRefMin*machRefMin - 1.)*(machRefMin*machRefMin - 1.) * uL*uL + 4.* machRefMin*machRefMin*cL*cL);
  // cL *= 0.5;
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

void ModUEq::addNonConsMrfFlux(Phase** phases)
{
  static_cast<FluxUEq*> (fluxBuffMRF)->addNonConsMrfFlux(phases);
}

//****************************************************************************

void ModUEq::reverseProjectionMrfFlux(const Coord normal, const Coord tangent, const Coord binormal) const
{
  Coord projectedFlux;
  projectedFlux.setX(normal.getX()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getX() + tangent.getX()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getY() + binormal.getX()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getZ());
  projectedFlux.setY(normal.getY()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getX() + tangent.getY()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getY() + binormal.getY()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getZ());
  projectedFlux.setZ(normal.getZ()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getX() + tangent.getZ()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getY() + binormal.getZ()*static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.getZ());
  static_cast<FluxUEq*> (fluxBuffMRF)->m_momentum.setXYZ(projectedFlux.getX(), projectedFlux.getY(), projectedFlux.getZ());
}

//****************************************************************************