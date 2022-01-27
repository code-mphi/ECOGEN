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

#include "BoundCondNullFlux.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondNullflux::BoundCondNullflux(const BoundCondNullflux& Source, const int& lvl) : BoundCond(Source, lvl)
{
}

//****************************************************************************

BoundCondNullflux::BoundCondNullflux(int numPhysique) : 
  BoundCond(numPhysique)
{
}

//****************************************************************************

BoundCondNullflux::~BoundCondNullflux() {}

//****************************************************************************

void BoundCondNullflux::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondNullflux(*(this)));
}

//****************************************************************************

void BoundCondNullflux::solveRiemannBoundary(Cell& /*cellLeft*/, const double& /*dxLeft*/, double& /*dtMax*/)
{
  model->solveRiemannNullFlux();
}

//****************************************************************************
//******************************AMR Method***********************************
//****************************************************************************

void BoundCondNullflux::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondNullflux(*this, m_lvl + 1));
}

//****************************************************************************