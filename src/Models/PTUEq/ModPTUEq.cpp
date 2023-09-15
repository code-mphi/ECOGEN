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

#include "ModPTUEq.h"
#include "PhasePTUEq.h"
#include "GradPhasePTUEq.h"
#include "GradMixPTUEq.h"

const std::string ModPTUEq::NAME = "TEMPERATUREPRESSUREVELOCITYEQ";

//***********************************************************************

ModPTUEq::ModPTUEq(const int& numbTransports, const int& numbPhases) :
  Model(NAME, numbTransports)
{
  fluxBuff = new FluxPTUEq(numbPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxPTUEq(numbPhases));
  }
}

//***********************************************************************

ModPTUEq::~ModPTUEq()
{
  delete fluxBuff;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//***********************************************************************

void ModPTUEq::allocateCons(Flux** cons)
{
  *cons = new FluxPTUEq(numberPhases);
}

//***********************************************************************

void ModPTUEq::allocatePhase(Phase** phase)
{
  *phase = new PhasePTUEq;
}

//***********************************************************************

void ModPTUEq::allocateMixture(Mixture** mixture)
{
  *mixture = new MixPTUEq;
}

//***********************************************************************

void ModPTUEq::allocatePhaseGradient(GradPhase** phase)
{
  *phase = new GradPhasePTUEq;
}

//***********************************************************************

void ModPTUEq::allocateMixtureGradient(GradMixture** mixture)
{
  *mixture = new GradMixPTUEq;
}

//***********************************************************************

void ModPTUEq::fulfillState(Phase** phases, Mixture* mixture)
{
  //Complete phases and mixture states from : alphak, pressure and temperature
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->setPressure(mixture->getPressure());
    phases[k]->setDensity(phases[k]->getEos()->computeDensity(mixture->getPressure(), mixture->getTemperature()));
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  mixture->computeMixtureVariables(phases);
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModPTUEq::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double> &boundData) const
{
  Phase* vecPhase;
  double sL, sR;
  double pStar(0.), rhoStar(0.), EStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getMixSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uR = cellRight.getMixture()->getVelocity().getX(), cR = cellRight.getMixture()->getMixSoundSpeed(), pR = cellRight.getMixture()->getPressure(), rhoR = cellRight.getMixture()->getDensity();

  //Davis
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
      static_cast<FluxPTUEq*> (fluxBuff)->m_mass[k] = alpha*density*uL;
    }
    double vitY = cellLeft.getMixture()->getVelocity().getY(); double vitZ = cellLeft.getMixture()->getVelocity().getZ();
    double totalEnergy = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setX(rhoL*uL*uL + pL);
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setY(rhoL*vitY*uL);
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setZ(rhoL*vitZ*uL);
    static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture = (rhoL*totalEnergy + pL)*uL;

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
      static_cast<FluxPTUEq*> (fluxBuff)->m_mass[k] = alpha*density*uR;
    }
    double vitY = cellRight.getMixture()->getVelocity().getY(); double vitZ = cellRight.getMixture()->getVelocity().getZ();
    double totalEnergy = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setX(rhoR*uR*uR + pR);
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setY(rhoR*vitY*uR);
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setZ(rhoR*vitZ*uR);
    static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture = (rhoR*totalEnergy + pR)*uR;

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
      mkL = alpha*density*(sL - uL);
      static_cast<FluxPTUEq*> (fluxBuff)->m_mass[k] = mkL / (sL - sM) * sM;
    }
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setY(rhoStar*vitY*sM);
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setZ(rhoStar*vitZ*sM);
    static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*sM;

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
      mkR = alpha*density*(sR - uR);
      static_cast<FluxPTUEq*> (fluxBuff)->m_mass[k] = mkR / (sR - sM) * sM;
    }
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setY(rhoStar*vitY*sM);
    static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setZ(rhoStar*vitZ*sM);
    static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*sM;

    // Boundary data for output
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vitY;
    boundData[VarBoundary::velW] = vitZ;
  }

  //Contact discontinuity velocity
  static_cast<FluxPTUEq*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions ****************
//****************************************************************************

void ModPTUEq::solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& boundData) const
{
  double sL;
  double pStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getMixSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();

  sL = std::min(uL - cL, -uL - cL);
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));

  pStar = rhoL*(uL - sL)*uL + pL;

  for (int k = 0; k < numberPhases; k++)
  {
    static_cast<FluxPTUEq*> (fluxBuff)->m_mass[k] = 0.;
  }
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setX(pStar);
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture = 0.;

  //Contact discontinuity velocity
  static_cast<FluxPTUEq*> (fluxBuff)->m_sM = 0.;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = 0.;
  boundData[VarBoundary::velU] = 0.;
  boundData[VarBoundary::velV] = 0.;
  boundData[VarBoundary::velW] = 0.;
}

//****************************************************************************

void ModPTUEq::solveRiemannInletTank(Cell& cellLeft, const double& dxLeft, double& dtMax, const double* ak0, const double* rhok0, const double& p0, const double& T0, std::vector<double> &boundData) const
{
  double sL, zL, sM, vmv0, mL;
  double pStar(0.), uStar(0.), rhoStar(0.), uyStar(0.), uzStar(0.), EStar(0.), vStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getMixSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity(), TL = cellLeft.getMixture()->getTemperature();
  double uyL = cellLeft.getMixture()->getVelocity().getY(), uzL = cellLeft.getMixture()->getVelocity().getZ(), EL = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();

  zL = rhoL*cL;
  for (int k = 0; k < numberPhases; k++) {
    TB->Yk[k] = cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() / rhoL;
  }

  //1) Left wave velocity estimation using pStar = p0
  //-------------------------------------------------
  pStar = p0;
  vStar = cellLeft.getMixture()->computeVolumeIsentrope(TB->Yk, pL, TL, pStar);
  vmv0 = vStar - 1. / rhoL;
  if (std::fabs(vmv0) > 1e-10) { mL = sqrt((pL - pStar) / vmv0); }
  else { mL = zL; }
  sL = uL - mL / rhoL;
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  sM = uL + mL * vmv0;

  //2) Check for pathologic cases
  //-----------------------------
  if (sL >= 0.) { //supersonic outflow => left state solution
    uStar = uL;
    pStar = pL;
    for (int k = 0; k < numberPhases; k++) {
      TB->YkStar[k] = cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() / rhoL;
    }
    rhoStar = rhoL;
    uyStar = uyL;
    uzStar = uzL;
    EStar = EL;
  }
  else if (sM >= -1e-3) { //subsonic outflow => star left state solution
    uStar = sM;
    pStar = p0; //approximation
    for (int k = 0; k < numberPhases; k++) {
      TB->YkStar[k] = cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() / rhoL;
    }
    rhoStar = 1. / vStar;
    uyStar = uyL;
    uzStar = uzL;
    EStar = EL + (uStar - uL)*(uStar + pL / mL);
  }

  //3) Tank
  //-------
  else { //tank inflow => star right state solution
    //Total enthalpy and entropy in tank state
    double H0(0.);
    double rho0 = cellLeft.getMixture()->computeDensity(ak0, rhok0);
    for (int k = 0;k < numberPhases;k++) {
      TB->Yk0[k] = ak0[k] * rhok0[k] / rho0;
      H0 += TB->Yk0[k] * TB->eos[k]->computeTotalEnthalpy(rhok0[k], p0, 0.);  //default zero velocity in tank
    }
    //ITERATIVE PROCESS FOR PRESSURE DETERMINATION 
    //--------------------------------------------
    int iteration(0);
    double p(0.5*p0);
    double f(0.), df(1.);
    double dmL;
    double uStarR(0.), duStarR(0.), uStarL(0.), duStarL(0.);
    double TStarR(0.);
    double hStarR(0.), dhStarR(0.);
    double vStarL(0.), dvStarL(0.);
    do {
      p -= f / df; iteration++;
      if (iteration > 50) Errors::errorMessage("solveRiemannInletTank not converged in ModPTUEq");
      //Physical pressure ?
      for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(p); }
      if (p > p0) { p = p0 - 1e-6; }
      //R) Tank rekations in the right (H=cte et s=cste)
      //mixture entropy constant brings the relation between Tr* andd p*
      TStarR = cellLeft.getMixture()->computeTemperatureIsentrope(TB->Yk0, p0, T0, p, &TStarR);
      hStarR = cellLeft.getMixture()->computeEnthalpyIsentrope(TB->Yk0, p0, T0, p, &dhStarR);
      uStarR = -sqrt(2.*(H0 - hStarR));
      duStarR = -dhStarR / uStarR;
      //L) Left relations s=cste (could be R-H if necessary)
      vStarL = cellLeft.getMixture()->computeVolumeIsentrope(TB->Yk, pL, TL, p, &dvStarL);
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
      if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
      uStarL = uL + mL*vmv0;
      duStarL = dmL*vmv0 + mL* dvStarL;
      //solved function
      f = uStarR - uStarL;
      df = duStarR - duStarL;
    } while (std::fabs(f)>1e-3); //End iterative loop
    pStar = p;
    uStar = 0.5*(uStarL + uStarR);
    rhoStar = 0.;
    for (int k = 0; k < numberPhases; k++) { 
      TB->YkStar[k] = TB->Yk0[k];
      rhoStar += TB->YkStar[k] / TB->eos[k]->computeDensity(pStar,TStarR);
    }
    rhoStar = 1. / rhoStar;
    uyStar = 0.;
    uzStar = 0.;
    EStar = H0 - pStar/rhoStar;
  } //End tank case

  //4) Flux completion
  //------------------
  for (int k = 0; k < numberPhases; k++) {
    static_cast<FluxPTUEq*> (fluxBuff)->m_mass[k] = rhoStar* TB->YkStar[k] * uStar;
  }
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setY(rhoStar*uStar*uyStar);
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setZ(rhoStar*uStar*uzStar);
  static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxPTUEq*> (fluxBuff)->m_sM = sM;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = uyStar;
  boundData[VarBoundary::velW] = uzStar;
}

//****************************************************************************

void ModPTUEq::solveRiemannOutletPressure(Cell& cellLeft, const double& dxLeft, double& dtMax, const double p0, std::vector<double> &boundData) const
{
  double sL, zL;
  double pStar(p0);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getMixSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity(), TL = cellLeft.getMixture()->getTemperature();;
  double vL = cellLeft.getMixture()->getVelocity().getY(), wL = cellLeft.getMixture()->getVelocity().getZ();

  for (int k = 0; k < numberPhases; k++) {
    TB->Yk[k] = cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() / rhoL;
  }

  //Left wave acoustic estimation
  //-----------------------------
  zL = rhoL*cL;
  double rhoStar(0.), vmv0, mL, uStar;
  rhoStar = 1./cellLeft.getMixture()->computeVolumeIsentrope(TB->Yk, pL, TL, p0);
  vmv0 = 1./ rhoStar - 1. / rhoL;
  if (std::fabs(vmv0) > 1e-10) {
    mL = sqrt((pL - p0) / vmv0);
  }
  else {
    mL = zL;
  }
  sL = uL - mL / rhoL;
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  uStar = uL + mL*vmv0;

  //Pathologic case sL>0            //FP//Q// Look for special case u<0
  if (sL >= 0.) { //Supersonic outflow => Left state solution
    uStar = uL;
    pStar = pL;
    rhoStar = rhoL;
  }

  //Flux completion
  double totalEnergy = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
  double EStar(totalEnergy + (uStar - uL)*(uStar - pL / mL));
  for (int k = 0; k < numberPhases; k++) {
    static_cast<FluxPTUEq*> (fluxBuff)->m_mass[k] = rhoStar * TB->Yk[k] * uStar ;
  }
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setX(uStar*uStar*rhoStar + pStar);
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setY(uStar*vL*rhoStar);
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.setZ(uStar*wL*rhoStar);
  static_cast<FluxPTUEq*> (fluxBuff)->m_energMixture = (EStar*rhoStar + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxPTUEq*> (fluxBuff)->m_sM = uStar;

  // Boundary data for output
  boundData[VarBoundary::p] = pStar;
  boundData[VarBoundary::rho] = rhoStar;
  boundData[VarBoundary::velU] = uStar;
  boundData[VarBoundary::velV] = vL;
  boundData[VarBoundary::velW] = wL;
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModPTUEq::selectScalar(Phase** phases, Mixture* mixture, Transport* transports, Variable nameVariable, int num) const
{
  switch (nameVariable) {
    case Variable::pressure:
      return mixture->getPressure();
      break;
    case Variable::temperature:
      return mixture->getTemperature();
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
    case Variable::density:
      return mixture->getDensity();
      break;
    default:
      Errors::errorMessage("nameVariable unknown in selectScalar"); return 0;
      break;
  }
}

//****************************************************************************

const double& ModPTUEq::getSM()
{
  return static_cast<FluxPTUEq*> (fluxBuff)->m_sM;
}

//****************************************************************************
//***************************** others methods *******************************
//****************************************************************************

void ModPTUEq::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  static_cast<FluxPTUEq*> (fluxBuff)->m_momentum.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************