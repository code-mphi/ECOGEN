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

#include "ModNEq.h"

const std::string ModNEq::NAME = "NONEQ";

//***********************************************************************

ModNEq::ModNEq(const int& numbTransports, const int& numbPhases) : Model(NAME, numbTransports)
{
  fluxBuff    = new FluxNEq(numbPhases);
  fluxBuffMRF = new FluxNEq(numbPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxNEq(numbPhases));
  }
}

//***********************************************************************

ModNEq::ModNEq(const std::string& name, const int& numbTransports) : Model(name, numbTransports) {}

//***********************************************************************

ModNEq::~ModNEq()
{
  delete fluxBuff;
  delete fluxBuffMRF;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//***********************************************************************

void ModNEq::allocateCons(Flux** cons) { *cons = new FluxNEq(numberPhases); }

//***********************************************************************

void ModNEq::allocatePhase(Phase** phase) { *phase = new PhaseNEq; }

//***********************************************************************

void ModNEq::allocateMixture(Mixture** mixture) { *mixture = new MixNEq; }

//***********************************************************************

void ModNEq::allocatePhaseGradient(GradPhase** phase) { *phase = new GradPhaseNEq; }

//***********************************************************************

void ModNEq::allocateMixtureGradient(GradMixture** mixture) { *mixture = new GradMixNEq; }

//***********************************************************************

void ModNEq::fulfillState(Phase** phases, Mixture* mixture)
{
  //Complete phases state
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase();
  }
  //Complete mixture variables using phases variable
  mixture->computeMixtureVariables(phases);
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModNEq::solveRiemannIntern(
  Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& boundData) const
{
  this->solver(cellLeft, cellRight, dxLeft, dxRight, dtMax, boundData);
}

// ****************************************************************************

//****************************************************************************
//*** Half Riemann solver for MRF interface between static/rotating region ***
//****************************************************************************

//FP//TODO

//****************************************************************************
//************** Half Riemann solvers for boundary conditions ****************
//****************************************************************************

//FP//TODO

//****************************************************************************
//********************** Transport Riemann solvers ***************************
//****************************************************************************

//FP//TODO

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModNEq::selectScalar(Phase** phases, Mixture* mixture, Transport* transports, Variable nameVariable, int num) const
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
    if (num < 0) {
      return mixture->getVelocity().getX();
    }
    else {
      return phases[num]->getU();
    }
    break;
  case Variable::velocityV:
    if (num < 0) {
      return mixture->getVelocity().getY();
    }
    else {
      return phases[num]->getV();
    }
    break;
  case Variable::velocityW:
    if (num < 0) {
      return mixture->getVelocity().getZ();
    }
    else {
      return phases[num]->getW();
    }
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
    Errors::errorMessage("nameVariable unknown in selectScalar");
    return 0;
    break;
  }
}

//***********************************************************************

const double& ModNEq::getSM() { return static_cast<FluxNEq*>(fluxBuff)->m_sM; }

//****************************************************************************
//***************************** Others methods *******************************
//****************************************************************************

void ModNEq::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  for (int k = 0; k < numberPhases; k++) {
    static_cast<FluxNEq*>(fluxBuff)->m_momentum[k].reverseProjection(normal, tangent, binormal);
  }
}

//****************************************************************************
