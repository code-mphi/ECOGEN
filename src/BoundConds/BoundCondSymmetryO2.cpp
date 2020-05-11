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

//! \file      BoundCondSymmetryO2.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCondSymmetryO2.h"

//****************************************************************************

BoundCondSymmetryO2::BoundCondSymmetryO2() {}

//****************************************************************************

BoundCondSymmetryO2::BoundCondSymmetryO2(const BoundCondSymmetryO2& Source, const int lvl) : BoundCondWallO2(Source, lvl)
{}

//****************************************************************************

BoundCondSymmetryO2::BoundCondSymmetryO2(int numPhysique) : BoundCondWallO2(numPhysique)
{}

//****************************************************************************

BoundCondSymmetryO2::~BoundCondSymmetryO2()
{}

//****************************************************************************

void BoundCondSymmetryO2::creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondSymmetryO2(*(this)));
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BoundCondSymmetryO2::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondSymmetryO2(*this, m_lvl + 1));
}

//****************************************************************************