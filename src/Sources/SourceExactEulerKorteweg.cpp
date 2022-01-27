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

#include "SourceExactEulerKorteweg.h"

//***********************************************************************

SourceExactEulerKorteweg::SourceExactEulerKorteweg(int physicalEntity) : SourceExact(physicalEntity)
{}

//***********************************************************************

SourceExactEulerKorteweg::~SourceExactEulerKorteweg()
{}

//***********************************************************************

void SourceExactEulerKorteweg::integrationExactSolution(Cell* cell, const double& dt)
{
  double rho(0.), eta(0.), omega(0.), root(0.);
  rho = cell->getPhase(0)->getDensity();
  eta = cell->getPhase(0)->getEta();
  omega = cell->getPhase(0)->getOmega();
  root = std::sqrt(1. / (alphaEK*betaEK*rho*rho));
  cell->getPhase(0)->setEta(rho + (eta - rho) * std::cos(dt * root) + omega / root * std::sin(dt * root));
  cell->getPhase(0)->setOmega(root * (rho - eta) * std::sin(dt * root) + omega * std::cos(dt * root));
}

//***********************************************************************