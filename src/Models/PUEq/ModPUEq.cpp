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
#include "ModPUEq.h"
#include "PhasePUEq.h"
#include "../../Relaxations/RelaxationPInfinite.h"

const std::string ModPUEq::NAME = "PRESSUREVELOCITYEQ";

//***********************************************************************

ModPUEq::ModPUEq(int& numberTransports, const int& numberPhases) : ModUEq(NAME,numberTransports)
{
  fluxBuff = new FluxPUEq(numberPhases);
  m_relaxations.push_back(new RelaxationPInfinite); //Pressure relaxation imposed in this model
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxPUEq(numberPhases));
  }
}

//***********************************************************************

ModPUEq::~ModPUEq(){}

//***********************************************************************

void ModPUEq::allocateCons(Flux** cons, const int& numberPhases)
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