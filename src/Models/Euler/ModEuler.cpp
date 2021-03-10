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
#include "ModEuler.h"
#include "PhaseEuler.h"

const std::string ModEuler::NAME = "EULER";

//****************************************************************************

ModEuler::ModEuler(const int& numberTransports) :
  Model(NAME, numberTransports)
{
  fluxBuff = new FluxEuler();
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxEuler());
  }
}

//****************************************************************************

ModEuler::~ModEuler()
{
  delete fluxBuff;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//****************************************************************************

void ModEuler::allocateCons(Flux** cons, const int& /*numberPhases*/)
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

void ModEuler::fulfillState(Phase** phases, Mixture* /*mixture*/, const int& /*numberPhases*/, Prim /*type*/)
{
  phases[0]->extendedCalculusPhase(phases[0]->getVelocity());
}

//***********************************************************************

//****************************************************************************
//********************* Cell to cell Riemann solvers** ***********************
//****************************************************************************

void ModEuler::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const int& /*numberPhases*/, const double& dxLeft, const double& dxRight, double& dtMax, double& massflow, double& powerFlux) const
{
  double cL, cR, sL, sR;
  double uL, uR, vL, vR, wL, wR, pL, pR, rhoL, rhoR, EL, ER;

  Phase* phaseLeft(0), *phaseRight(0);
  phaseLeft = cellLeft.getPhase(0);
  phaseRight = cellRight.getPhase(0);
  Eos* eos(0);
  eos = phaseLeft->getEos();

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

  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);

  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoL*(sL - uL)), mR(rhoR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));
  if (std::fabs(sM)<1.e-8) sM = 0.;

  if (sL > 0.){
    static_cast<FluxEuler*> (fluxBuff)->m_masse = rhoL*uL;
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(rhoL*uL*uL + pL);
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(rhoL*vL*uL);
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(rhoL*wL*uL);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoL*EL + pL)*uL;

    //Specific power flux through interface (W.m-2)
    powerFlux = rhoL * uL * eos->computeTotalEnthalpy(rhoL, pL, uL);
  }
  else if (sR < 0.){
    static_cast<FluxEuler*> (fluxBuff)->m_masse = rhoR*uR;
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(rhoR*uR*uR + pR);
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(rhoR*vR*uR);
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(rhoR*wR*uR);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoR*ER + pR)*uR;

    //Specific power flux through interface (W.m-2)
    powerFlux = rhoR * uR * eos->computeTotalEnthalpy(rhoR, pR, uR);
  }

  ////1) Option HLL
  //else if (std::fabs(sR - sL)>1.e-3)
  //{
  //  static_cast<FluxEuler*> (fluxBuff)->m_masse = (rhoR*uR*sL - rhoL*uL*sR + sL*sR*(rhoL - rhoR)) / (sL - sR);
  //  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(((rhoR*uR*uR + pR)*sL - (rhoL*uL*uL + pL)*sR + sL*sR*(rhoL*uL - rhoR*uR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY((rhoR*uR*vR*sL - rhoL*uL*vL*sR + sL*sR*(rhoL*vL - rhoR*vR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ((rhoR*uR*wR*sL - rhoL*uL*wL*sR + sL*sR*(rhoL*wL - rhoR*wR)) / (sL - sR));
  //  static_cast<FluxEuler*> (fluxBuff)->m_energ = ((rhoR*ER + pR)*uR*sL - (rhoL*EL + pL)*uL*sR + sL*sR*(rhoL*EL - rhoR*ER)) / (sL - sR);
  //}

  //2) Option HLLC
  else if (sM >= 0.) {
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    static_cast<FluxEuler*> (fluxBuff)->m_masse = rhoStar*sM;
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM+pStar);
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(rhoStar*sM*vL);
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(rhoStar*sM*wL);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;

    //Specific power flux through interface (W.m-2)
    powerFlux = rhoStar * sM * eos->computeTotalEnthalpy(rhoStar, pStar, sM);
  }
  else {
    double pStar = mR*(sM - uR) + pR;
    double rhoStar = mR / (sR - sM);
    double Estar = ER + (sM - uR)*(sM + pR / mR);
    static_cast<FluxEuler*> (fluxBuff)->m_masse = rhoStar*sM;
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(rhoStar*sM*vR);
    static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(rhoStar*sM*wR);
    static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;

    //Specific power flux through interface (W.m-2)
    powerFlux = rhoStar * sM * eos->computeTotalEnthalpy(rhoStar, pStar, sM);
  }

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = sM;

  //Specific mass flow through interface (kg.s-1.m-2)
  massflow = static_cast<FluxEuler*> (fluxBuff)->m_masse;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions** **************
//****************************************************************************

void ModEuler::solveRiemannWall(Cell& cellLeft, const int& /*numberPhases*/, const double& dxLeft, double& dtMax) const
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

  sL = std::min(uL - cL, -uL - cL);
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));

  pStar = rhoL*uL*(uL - sL) + pL;

  static_cast<FluxEuler*> (fluxBuff)->m_masse = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(0.);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(0.);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = 0.;

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = 0.;
}

//****************************************************************************

void ModEuler::solveRiemannInflow(Cell& cellLeft, const int& /*numberPhases*/, const double& dxLeft, double& dtMax, const double m0, const double* /*ak0*/, const double* rhok0, const double* pk0, double& massflow, double& powerFlux) const
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
  zL = rhoL*cL;

  sL = uL - cL;
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));

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
      if (iteration > 50) Errors::errorMessage("solveRimannInflow not converged in modEuler");
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

  static_cast<FluxEuler*> (fluxBuff)->m_masse = rhoStar*uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(rhoStar*uStar*vL);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(rhoStar*uStar*wL);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vL*vL + wL*wL)) + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;
  
  //Specific mass flow through interface (kg.s-1.m-2)
  massflow = static_cast<FluxEuler*> (fluxBuff)->m_masse;

  //Specific power flux through interface (W.m-2)
  powerFlux = massflow * eos->computeTotalEnthalpy(rhoStar, pStar, uStar);
}

//****************************************************************************

void ModEuler::solveRiemannSubInj(Cell& cellLeft, const int& /*numberPhases*/, const double& dxLeft, double& dtMax, const double m0, const double T0, double& massflow, double& powerFlux) const
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
  zL = rhoL * cL;

  sL = uL - cL;
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));

  double uRef(m0 / rhoL); // N-R reference velocity (don't use directly uL to avoid /0)

  int it(0);
  double p(pL), f(0.), df(1.);
  double mL(zL), dmL(0.), vSmvL(0.);
  double rhoStarR(0.), uStarR(0.), drhoStarR(0.), duStarR(0.);
  double rhoStarL(0.), uStarL(0.), drhoStarL(0.), duStarL(0.), vStarL(0.), dvStarL(0.);
  double eStar(0.), rhoStar(0.), pStar(0.), uStar(0.);

  // ITERATIVE PROCESS FOR PRESSURE DETERMINATION
  // --------------------------------------------
  do {
    p -= f / df;
    it++;
    if (it > 50) Errors::errorMessage("solveRiemannSubInj not converged in modEuler");
    eos->verifyAndModifyPressure(p);

    // Right intermediate state 
    rhoStarR = eos->computeDensity(p, T0);
    drhoStarR = eos->drhodpcT(p, T0);
    uStarR = m0 / rhoStarR;
    duStarR = - m0 * drhoStarR / (rhoStarR * rhoStarR);

    // Left intermediate state
    rhoStarL = eos->computeDensityHugoniot(pL, rhoL, p, &drhoStarL);
    vStarL = 1. / rhoStarL;
    dvStarL = - drhoStarL / (rhoStarL * rhoStarL);
    vSmvL = vStarL - 1. / rhoL;
    if (std::fabs(vSmvL) > 1.e-10) { // Rankine-Hugoniot
      mL = std::sqrt((pL - p) / vSmvL);
      dmL = 0.5 * (- vSmvL + (p - pL) * dvStarL) / (vSmvL * vSmvL) / mL;
      duStarL = dmL * vSmvL + mL * dvStarL;
    }
    else { // Acoustic
      mL = zL;
      duStarL = zL * dvStarL;
    }
    uStarL = uL + mL * vSmvL;

    f = (uStarL - uStarR);
    df = (duStarL - duStarR);
  } while (std::fabs(f/uRef) > 1.e-5);

  pStar = p;
  rhoStar = rhoStarR;
  eStar = eos->computeEnergy(rhoStar, pStar);
  uStar = 0.5 * (uStarL + uStarR);

  static_cast<FluxEuler*> (fluxBuff)->m_masse = rhoStar * uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(rhoStar * uStar * uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(rhoStar * uStar * vL);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(rhoStar * uStar * wL);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar * (eStar + 0.5 * (uStar * uStar + vL * vL + wL * wL)) + pStar) * uStar;

  // Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;
  
  //Specific mass flow through interface (kg.s-1.m-2)
  massflow = static_cast<FluxEuler*> (fluxBuff)->m_masse;

  //Specific power flux through interface (W.m-2)
  powerFlux = massflow * eos->computeTotalEnthalpy(rhoStar, pStar, uStar);
}

//****************************************************************************

void ModEuler::solveRiemannTank(Cell& cellLeft, const int& /*numberPhases*/, const double& dxLeft, double& dtMax, const double* /*ak0*/, const double* rhok0, const double& p0, const double& /*T0*/, double& massflow, double& powerFlux) const
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
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
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
      if (iteration > 50) Errors::errorMessage("solveRiemannTank not converged in modEuler");
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
      if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
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

  static_cast<FluxEuler*> (fluxBuff)->m_masse = rhoStar*uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(rhoStar*uStar*vStar);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(rhoStar*uStar*wStar);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vStar*vStar + wStar*wStar)) + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;

  //Specific mass flow rate output (kg.s-1.m-2)
  massflow = static_cast<FluxEuler*> (fluxBuff)->m_masse;

  //Specific power flux through interface (W.m-2)
  powerFlux = massflow * eos->computeTotalEnthalpy(rhoStar, pStar, uStar);
}

//****************************************************************************

void ModEuler::solveRiemannOutflow(Cell& cellLeft, const int& /*numberPhases*/, const double& dxLeft, double& dtMax, const double p0, double& massflow, double& powerFlux) const
{
  double cL, sL, zL;
  double uL, pL, rhoL, vL, wL;
  double uStar(0.), rhoStar(0.), pStar(p0), eStar(0.);

  Phase* phaseLeft(0);
  phaseLeft = cellLeft.getPhase(0);
  Eos* eos(0);
  eos = phaseLeft->getEos();
  
  uL = phaseLeft->getU();
  vL = phaseLeft->getV();
  wL = phaseLeft->getW();
  pL = phaseLeft->getPressure();
  rhoL = phaseLeft->getDensity();
  cL = phaseLeft->getSoundSpeed();
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
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
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
  static_cast<FluxEuler*> (fluxBuff)->m_masse = rhoStar*uStar;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setX(rhoStar*uStar*uStar + pStar);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setY(rhoStar*uStar*vL);
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setZ(rhoStar*uStar*wL);
  static_cast<FluxEuler*> (fluxBuff)->m_energ = (rhoStar*(eStar + 0.5*(uStar*uStar + vL*vL + wL*wL)) + pStar)*uStar;

  //Contact discontinuity velocity
  static_cast<FluxEuler*> (fluxBuff)->m_sM = uStar;

  //Specific mass flow rate output (kg.s-1.m-2)
  massflow = static_cast<FluxEuler*> (fluxBuff)->m_masse;

  //Specific power flux through interface (W.m-2)
  powerFlux = massflow * eos->computeTotalEnthalpy(rhoStar, pStar, uStar);
}

//****************************************************************************

const double& ModEuler::getSM()
{
  return static_cast<FluxEuler*> (fluxBuff)->m_sM;
}

//****************************************************************************

void ModEuler::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  Coord fluxProjete;
  fluxProjete.setX(normal.getX()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getX() + tangent.getX()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getY() + binormal.getX()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getZ());
  fluxProjete.setY(normal.getY()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getX() + tangent.getY()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getY() + binormal.getY()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getZ());
  fluxProjete.setZ(normal.getZ()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getX() + tangent.getZ()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getY() + binormal.getZ()*static_cast<FluxEuler*> (fluxBuff)->m_qdm.getZ());
  static_cast<FluxEuler*> (fluxBuff)->m_qdm.setXYZ(fluxProjete.getX(), fluxProjete.getY(), fluxProjete.getZ());
}

//****************************************************************************