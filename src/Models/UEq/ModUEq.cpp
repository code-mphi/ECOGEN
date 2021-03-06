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
#include "ModUEq.h"
#include "PhaseUEq.h"

const std::string ModUEq::NAME = "VELOCITYEQ";

//***********************************************************************

ModUEq::ModUEq(int& numberTransports, const int& numberPhases) : Model(NAME,numberTransports)
{
  fluxBuff = new FluxUEq(numberPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxUEq(numberPhases));
  }
}

//***********************************************************************

ModUEq::ModUEq(const std::string& name, const int& numberTransports) : Model(name,numberTransports){}

//***********************************************************************

ModUEq::~ModUEq()
{
  delete fluxBuff;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//***********************************************************************

void ModUEq::allocateCons(Flux** cons, const int& numberPhases)
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

void ModUEq::fulfillState(Phase** phases, Mixture* mixture, const int& numberPhases, Prim type)
{
  //Specific to restart simulation
  if (type == restart) {
    for (int k = 0; k < numberPhases; k++) { phases[k]->setPressure(mixture->getPressure()); } //KS// Only works for PUEq, to update and maybe improve to avoid this test (put an additional method for restart)
  }
  //Complete phases state
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  //Complete mixture variables using phases variable
  mixture->computeMixtureVariables(phases, numberPhases);
}

//***********************************************************************

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModUEq::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const int& numberPhases, const double& dxLeft, const double& dxRight, double& dtMax, double& massflow, double& powerFlux) const
{
  Phase* vecPhase;
  double sL, sR;
  double pStar(0.), rhoStar(0.), EStar(0.);

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
      double energie = vecPhase->getEnergy();
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = alpha*density*uL;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha*density*energie*uL;
    }
    double vitY = cellLeft.getMixture()->getVelocity().getY(); double vitZ = cellLeft.getMixture()->getVelocity().getZ();
    double totalEnergy = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setX(rhoL*uL*uL + pL);
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setY(rhoL*vitY*uL);
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setZ(rhoL*vitZ*uL);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoL*totalEnergy + pL)*uL;

    //Specific power flux through interface (W.m-2)
    powerFlux = rhoL * uL * (totalEnergy + pL / rhoL);
  }
  else if (sR <= 0.){
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellRight.getPhase(k);
      double alpha = vecPhase->getAlpha();
      double density = vecPhase->getDensity();
      double energie = vecPhase->getEnergy();
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = alpha*density*uR;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha*density*energie*uR;
    }
    double vitY = cellRight.getMixture()->getVelocity().getY(); double vitZ = cellRight.getMixture()->getVelocity().getZ();
    double totalEnergy = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setX(rhoR*uR*uR + pR);
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setY(rhoR*vitY*uR);
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setZ(rhoR*vitZ*uR);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoR*totalEnergy + pR)*uR;

    //Specific power flux through interface (W.m-2)
    powerFlux = rhoR * uR * (totalEnergy + pR / rhoR);
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
      TB->pkStar[k] = TB->eos[k]->computePressureIsentropic(pressure, density, TB->rhokStar[k]);
      // TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
      TB->ekStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]);
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = alpha* TB->rhokStar[k] * sM;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha* TB->rhokStar[k] * TB->ekStar[k] * sM;
    }
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setY(rhoStar*vitY*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setZ(rhoStar*vitZ*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*sM;

    //Specific power flux through interface (W.m-2)
    powerFlux = rhoStar * sM * (EStar + pStar / rhoStar);
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
      TB->pkStar[k] = TB->eos[k]->computePressureIsentropic(pressure, density, TB->rhokStar[k]);
      // TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
      TB->ekStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]);
      static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = alpha*sM;
      static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = alpha* TB->rhokStar[k] * sM;
      static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = alpha* TB->rhokStar[k] * TB->ekStar[k] * sM;
    }
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setY(rhoStar*vitY*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_qdm.setZ(rhoStar*vitZ*sM);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*sM;

    //Specific power flux through interface (W.m-2)
    powerFlux = rhoStar * sM * (EStar + pStar / rhoStar);
  }

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = sM;

  //Specific mass flow through interface (kg.s-1.m-2)
  massflow = 0.;
  for (int k = 0; k < numberPhases; k++) {
    massflow += static_cast<FluxUEq*> (fluxBuff)->m_masse[k];
  }
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions ****************
//****************************************************************************

void ModUEq::solveRiemannWall(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax) const
{
  double sL;
  double pStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getFrozenSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();

  sL = std::min(uL - cL, -uL - cL);
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));

  pStar = rhoL*(uL - sL)*uL + pL;

  for (int k = 0; k < numberPhases; k++)
  {
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setX(pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setY(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setZ(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = 0.;
}

//****************************************************************************

void ModUEq::solveRiemannInflow(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax, const double m0, const double* ak0, const double* rhok0, const double* pk0, double& massflow, double& powerFlux) const
{
  double sL, zL;
  double pStar(0.), uStar(0.), rhoStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getFrozenSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double vL = cellLeft.getMixture()->getVelocity().getY(), wL = cellLeft.getMixture()->getVelocity().getZ();

  //Compute total enthalpy of injected fluid and speed of sound
  double rho0 = cellLeft.getMixture()->computeDensity(ak0, rhok0, numberPhases);
  double u0 = m0 / rho0;
  double c0(0.), p0(0.);
  for (int k = 0;k < numberPhases;k++) {
    p0 += ak0[k] * pk0[k];
    TB->Hk0[k] = TB->eos[k]->computeTotalEnthalpy(rhok0[k], pk0[k], u0);
    TB->Yk0[k] = ak0[k] * rhok0[k] / rho0;
    double ck = cellLeft.getPhase(k)->getEos()->computeSoundSpeed(rhok0[k], pk0[k]);
    c0 += ak0[k]/ std::max((rhok0[k] * ck * ck), epsilonAlphaNull);
  }
  c0 = 1./ sqrt(rho0 * c0);

  //Estimates for acoustic wave sL
  sL = uL - cL;
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  zL = rhoL*cL;

  //Null Mass flow
  //--------------
  if (fabs(u0) < 1.e-6) {
    uStar = 0.;
    pStar = pL;
    rhoStar = rhoL;
    for (int k = 0; k < numberPhases; k++) {
      TB->vkStar[k] = 1.;
    }
  }
  //Supersonic inflow
  //-----------------
  else if (u0 < -c0) {
    uStar = u0;
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

    do {
      pStar -= f / df; iteration++;
      if (iteration > 50) Errors::errorMessage("solveRiemannInflow not converged in ModUEq");
      //Physical pressure ?
      for (int k = 0; k < numberPhases; k++) {
        TB->eos[k]->verifyAndModifyPressure(pStar);
      }
      //Left acoustic relations
      u = uL + (pL - pStar) / zL;
      if (u >= -1e-6) u = -1e-6;
      du = -1. / zL;
      f = u / u0; df = du / u0;
      //Compute from m0, Hk0, Yk0 on the right
      for (int k = 0; k < numberPhases; k++) {
        hk = TB->Hk0[k] - 0.5 * u * u;
        TB->vkStar[k] = TB->eos[k]->vfpfh(pStar, hk);
        double dvk = TB->eos[k]->dvdpch(pStar, hk) - TB->eos[k]->dvdhcp(pStar) * u * du;
        f -= TB->Yk0[k] * TB->vkStar[k] * rho0;
        df -= TB->Yk0[k] * dvk * rho0;
      }
    } while (std::fabs(f) > 1e-8 && iteration <= 50);
    uStar = u;
    rhoStar = m0 / uStar;
  }

  //Flux completion
  double Estar(0.5*(uStar * uStar + vL*vL + wL*wL)), ek, rhok;
  for (int k = 0; k<numberPhases; k++) {
    rhok = 1. / TB->vkStar[k];
    ek = TB->eos[k]->computeEnergy(rhok, pStar); Estar += TB->Yk0[k] * ek;
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = TB->Yk0[k] * TB->vkStar[k] * rhoStar * uStar;
    static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] * rhok;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] * rhok*ek;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setX(uStar * uStar * rhoStar + pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setY(uStar *vL * rhoStar);
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setZ(uStar *wL * rhoStar);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (Estar * rhoStar + pStar)* uStar;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = uStar;

  //Specific mass flow through interface (kg.s-1.m-2)
  massflow = rhoStar * uStar;

  //Specific power flux through interface (W.m-2)
  powerFlux = massflow * (Estar + pStar / rhoStar);
}

//****************************************************************************

void ModUEq::solveRiemannTank(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax, const double* ak0, const double* rhok0, const double& p0, const double& /*T0*/, double& massflow, double& powerFlux) const
{
  double tabp[50], tabf[50];
  double sL, zL, sM, vmv0, mL;
  double pStar(0.), uStar(0.), rhoStar(0.), vStar(0.), uyStar(0.), uzStar(0.);
  Phase* vecPhase;

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getFrozenSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uyL = cellLeft.getMixture()->getVelocity().getY(), uzL = cellLeft.getMixture()->getVelocity().getZ();

  zL = rhoL*cL;

  //1) Left wave velocity estimation using pStar = p0
  //-------------------------------------------------
  pStar = p0; vStar = 0.;
  for (int k = 0; k < numberPhases; k++) {
    vecPhase = cellLeft.getPhase(k);
    //TB->rhokStar[k] = TB->eos[k]->computeDensityIsentropic(vecPhase->getPressure(), vecPhase->getDensity(), pStar); //other possiblity
    TB->rhokStar[k] = TB->eos[k]->computeDensityHugoniot(vecPhase->getPressure(), vecPhase->getDensity(), pStar);
    vStar += vecPhase->getAlpha()*vecPhase->getDensity() / rhoL / std::max(TB->rhokStar[k], epsilonAlphaNull);
  }
  vmv0 = vStar - 1. / rhoL;
  if (std::fabs(vmv0) > 1e-10) { mL = sqrt((pL - pStar) / vmv0); }
  else { mL = zL; }
  sL = uL - mL / rhoL;
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  sM = uL + mL*vmv0;

  //2) Check for pathologic cases
  //-----------------------------
  if (sL >= 0.) { //supersonic outflow => left state solution
    uStar = uL;
    pStar = pL;
    for (int k = 0; k < numberPhases; k++) {
      vecPhase = cellLeft.getPhase(k);
      TB->rhokStar[k] = vecPhase->getDensity();
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
    double rho0 = cellLeft.getMixture()->computeDensity(ak0, rhok0, numberPhases);
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
        Errors::errorMessage("solveRiemannTank not converged in ModUEq");
      }
      //R) Tank rekations in the right (H=cte et sk=cste)
      uStarR = H0; duStarR = 0.;
      for (int k = 0; k < numberPhases; k++) {
        TB->rhokStar[k] = TB->eos[k]->computeDensityIsentropic(p0, rhok0[k], p);
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
      if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
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

  } //End tank case

  //4) Flux completion
  //------------------
  double EStar(0.5*(uStar*uStar + uyStar*uyStar + uzStar*uzStar)), ek;
  for (int k = 0; k < numberPhases; k++) {
    ek = TB->eos[k]->computeEnergy(TB->rhokStar[k], pStar); EStar += TB->YkStar[k] * ek;
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = TB->YkStar[k] * rhoStar / std::max(TB->rhokStar[k], epsilonAlphaNull) * uStar;
    static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] * TB->rhokStar[k];
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->m_masse[k] * ek;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setY(rhoStar*uStar*uyStar);
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setZ(rhoStar*uStar*uzStar);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (rhoStar*EStar + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = uStar;

  //Specific mass flow rate output (kg.s-1.m-2)
  massflow = rhoStar * uStar;

  //Specific power flux through interface (W.m-2)
  powerFlux = massflow * (EStar + pStar / rhoStar);
}

//****************************************************************************

void ModUEq::solveRiemannOutflow(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax, const double p0, double& massflow, double& powerFlux) const
{
  double sL, sM, zL;
  double pStar(p0), EStar(0.), vStar(0.), uStar(0.);

  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = cellLeft.getMixture()->getFrozenSoundSpeed(), pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double uyL = cellLeft.getMixture()->getVelocity().getY(), uzL = cellLeft.getMixture()->getVelocity().getZ();
  zL = rhoL*cL;

  //Left wave : isentropic wave assumption
  //--------------------------------------
  double vSmvL, mL(zL);
  for (int k = 0; k < numberPhases; k++) {
    //TB->rhokStar[k] = TB->eos[k]->computeDensityIsentropic(pL, vecPhase->getDensity(), pStar); //other possiblity
    TB->rhokStar[k] = TB->eos[k]->computeDensityHugoniot(pL, cellLeft.getPhase(k)->getDensity(), pStar);
    vStar += cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() /rhoL / std::max(TB->rhokStar[k], epsilonAlphaNull);
  }
  vSmvL = vStar - 1. / rhoL;
  if (abs(vSmvL) > 1e-10) { mL = sqrt((pL - pStar) / vSmvL); }
  sL = uL - mL / rhoL;
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
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
  double ekStar;
  EStar = 0.5*(uStar*uStar + uyL * uyL + uzL * uzL);
  for (int k = 0; k < numberPhases; k++) {
    double YkL = cellLeft.getPhase(k)->getAlpha()*cellLeft.getPhase(k)->getDensity() / rhoL;
    ekStar = TB->eos[k]->computeEnergy(TB->rhokStar[k], pStar);
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = YkL / std::max(TB->rhokStar[k], epsilonAlphaNull) / vStar * uStar;
    static_cast<FluxUEq*> (fluxBuff)->m_masse[k] = static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] * TB->rhokStar[k];
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = static_cast<FluxUEq*> (fluxBuff)->m_masse[k] * ekStar;
    EStar += YkL * ekStar;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setX(uStar*uStar / vStar + pStar);
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setY(uStar*uyL / vStar);
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setZ(uStar*uzL / vStar);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = (EStar / vStar + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxUEq*> (fluxBuff)->m_sM = uStar;

  //Specific mass flow rate output (kg.s-1.m-2)
  massflow = uStar / vStar;
  
  //Specific power flux through interface (W.m-2)
  powerFlux = massflow * (EStar + pStar * vStar);
}

//****************************************************************************
//********************** Transport Riemann solvers ***************************
//****************************************************************************

void ModUEq::solveRiemannTransportIntern(Cell& cellLeft, Cell& cellRight, const int& numberTransports)
{
	for (int k = 0; k < numberTransports; k++) {
		fluxBufferTransport[k].solveRiemann(cellLeft.getTransport(k).getValue(), cellRight.getTransport(k).getValue(), static_cast<FluxUEq*> (fluxBuff)->m_sM);
	}
}

//****************************************************************************

void ModUEq::solveRiemannTransportWall(const int& numberTransports)
{
	for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannWall();
	}
}

//****************************************************************************

void ModUEq::solveRiemannTransportInflow(Cell& cellLeft, const int& numberTransports, double* valueTransports)
{
	for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannInflow(cellLeft.getTransport(k).getValue(), static_cast<FluxUEq*> (fluxBuff)->m_sM, valueTransports[k]);
	}
}

//****************************************************************************

void ModUEq::solveRiemannTransportTank(Cell& cellLeft, const int& numberTransports, double* valueTransports)
{
	for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannTank(cellLeft.getTransport(k).getValue(), static_cast<FluxUEq*> (fluxBuff)->m_sM, valueTransports[k]);
	}
}

//****************************************************************************

void ModUEq::solveRiemannTransportOutflow(Cell& cellLeft, const int& numberTransports, double* valueTransports)
{
	for (int k = 0; k < numberTransports; k++) {
    fluxBufferTransport[k].solveRiemannOutflow(cellLeft.getTransport(k).getValue(), static_cast<FluxUEq*> (fluxBuff)->m_sM, valueTransports[k]);
	}
}

//****************************************************************************

const double& ModUEq::getSM()
{
  return static_cast<FluxUEq*> (fluxBuff)->m_sM;
}

//****************************************************************************
//***************************** others methods *******************************
//****************************************************************************

void ModUEq::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  Coord fluxProjected;
  fluxProjected.setX(normal.getX()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getX() + tangent.getX()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getY() + binormal.getX()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getZ());
  fluxProjected.setY(normal.getY()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getX() + tangent.getY()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getY() + binormal.getY()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getZ());
  fluxProjected.setZ(normal.getZ()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getX() + tangent.getZ()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getY() + binormal.getZ()*static_cast<FluxUEq*> (fluxBuff)->m_qdm.getZ());
  static_cast<FluxUEq*> (fluxBuff)->m_qdm.setXYZ(fluxProjected.getX(), fluxProjected.getY(), fluxProjected.getZ());
}

//****************************************************************************