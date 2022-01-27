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

#include "ModNonLinearSchrodinger.h"
#include "PhaseNonLinearSchrodinger.h"

const std::string ModNonLinearSchrodinger::NAME = "NONLINEARSCHRODINGER";

//****************************************************************************

ModNonLinearSchrodinger::ModNonLinearSchrodinger(const int& numbTransports, const double& alpha, const double &beta) :
  ModEulerKorteweg(NAME, numbTransports, alpha, beta, 0., 0.)
{
  fluxBuff = new FluxNonLinearSchrodinger();
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxNonLinearSchrodinger());
  }
}

//****************************************************************************

ModNonLinearSchrodinger::~ModNonLinearSchrodinger() {}

//****************************************************************************

void ModNonLinearSchrodinger::allocateCons(Flux** cons)
{
  *cons = new FluxNonLinearSchrodinger;
}

//****************************************************************************

void ModNonLinearSchrodinger::allocatePhase(Phase** phase)
{
  *phase = new PhaseNonLinearSchrodinger;
}

//****************************************************************************

void ModNonLinearSchrodinger::allocateMixture(Mixture** mixture)
{
  *mixture = new MixNonLinearSchrodinger;
}

//****************************************************************************

double ModNonLinearSchrodinger::kappa(const double& density) const
{
  return 1./(4.*density);
}

//****************************************************************************

double ModNonLinearSchrodinger::kappaPrime(const double& density) const
{
  return -1./(4.*density*density);
}

//****************************************************************************

double ModNonLinearSchrodinger::kappaSecond(const double& density) const
{
  return 1./(2.*density*density*density);
}

//****************************************************************************

double ModNonLinearSchrodinger::epsilonPrime(Cell& /*cell*/, const double& /*density*/) const
{
  return 0.5;
}

//****************************************************************************

double ModNonLinearSchrodinger::epsilonSecond(Cell& /*cell*/, const double& /*density*/) const
{
  return 0.;
}

//****************************************************************************