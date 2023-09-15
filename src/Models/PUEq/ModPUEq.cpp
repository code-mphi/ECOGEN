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

#include "ModPUEq.h"
#include "PhasePUEq.h"
#include "../../Relaxations/RelaxationPInfinite.h"

const std::string ModPUEq::NAME = "PRESSUREVELOCITYEQ";

//***********************************************************************

ModPUEq::ModPUEq(const int& numbTransports, const int& numbPhases) : ModUEq(NAME, numbTransports)
{
  fluxBuff = new FluxPUEq(numbPhases);
  fluxBuffMRF = new FluxPUEq(numbPhases);
  m_relaxations.push_back(new RelaxationPInfinite); //Pressure relaxation imposed in this model
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxPUEq(numbPhases));
  }
}

//***********************************************************************

ModPUEq::~ModPUEq(){}

//***********************************************************************

void ModPUEq::allocateCons(Flux** cons)
{
  *cons = new FluxPUEq(numberPhases);
}

//***********************************************************************

void ModPUEq::allocatePhase(Phase** phase)
{
  *phase = new PhasePUEq;
}

//***********************************************************************

void ModPUEq::allocateMixture(Mixture** mixture)
{
  *mixture = new MixPUEq;
}

//***********************************************************************

void ModPUEq::fulfillStateRestart(Phase** phases, Mixture* mixture)
{
  for (int k = 0; k < numberPhases; k++) { 
    phases[k]->setPressure(mixture->getPressure());
    phases[k]->verifyAndCorrectPhase();
  }
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModPUEq::selectScalar(Phase** phases, Mixture* mixture, Transport* transports, Variable nameVariable, int num) const
{
  switch (nameVariable) {
    case Variable::pressure:
      return mixture->getPressure();
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
      //double psi(0.), coeff(0.75);
      //psi = std::pow(transports[num].getValue(), coeff) / (std::pow(transports[num].getValue(), coeff) + std::pow((1 - transports[num].getValue()), coeff));
      //return psi;
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