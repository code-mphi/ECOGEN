//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

//! \file      BoundCondWall.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCondWall.h"

//****************************************************************************

BoundCondWall::BoundCondWall() {}

//****************************************************************************

BoundCondWall::BoundCondWall(const BoundCondWall& Source, const int lvl) : BoundCond(Source)
{
  m_lvl = lvl;
}

//****************************************************************************

BoundCondWall::BoundCondWall(int numPhysique) : BoundCond(numPhysique)
{}

//****************************************************************************

BoundCondWall::~BoundCondWall() {}

//****************************************************************************

void BoundCondWall::creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondWall(*(this)));
}

//****************************************************************************

void BoundCondWall::solveRiemannLimite(Cell &cellLeft, const int & numberPhases, const double & dxLeft, double & dtMax)
{
  m_mod->solveRiemannWall(cellLeft, numberPhases, dxLeft, dtMax);
}

//****************************************************************************

void BoundCondWall::solveRiemannTransportLimite(Cell &cellLeft, const int & numberTransports) const
{
  m_mod->solveRiemannTransportWall(numberTransports);
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BoundCondWall::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondWall(*this, m_lvl + 1));
}

//****************************************************************************