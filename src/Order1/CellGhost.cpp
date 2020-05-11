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

//! \file      CellGhost.cpp
//! \author    K. Schmidmayer, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include "CellGhost.h"

//***********************************************************************

CellGhost::CellGhost() : Cell(), m_rankOfNeighborCPU(0) {}

//***********************************************************************

CellGhost::CellGhost(int lvl) : Cell(lvl), m_rankOfNeighborCPU(0) {}

//***********************************************************************

CellGhost::~CellGhost() {}

//***********************************************************************

void CellGhost::createChildCell(const int &lvl)
{
  m_childrenCells.push_back(new CellGhost(lvl + 1));
  m_childrenCells.back()->setRankOfNeighborCPU(m_rankOfNeighborCPU);
}

//***************************************************************************

int CellGhost::getRankOfNeighborCPU() const
{
  return m_rankOfNeighborCPU;
}

//***************************************************************************

void CellGhost::setRankOfNeighborCPU(int rank)
{
  m_rankOfNeighborCPU = rank;
}

//***************************************************************************
