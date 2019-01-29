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

//! \file      QAPConductivity.cpp
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "QAPConductivity.h"
#include <iostream>

using namespace std;

//***********************************************************************

QAPConductivity::QAPConductivity(){}

//***********************************************************************

QAPConductivity::QAPConductivity(AddPhys* addPhys, const int &numberPhases) : QuantitiesAddPhys(addPhys)
{
  m_gradTk = new Coord[numberPhases];
  for (int k = 0; k < numberPhases; k++) {
    m_gradTk[k] = 0.;
  }
}

//***********************************************************************

QAPConductivity::~QAPConductivity(){  delete[] m_gradTk; }

//***********************************************************************

void QAPConductivity::computeQuantities(Cell* cell)
{
  for (int k = 0; k < cell->getNumberPhases(); k++) {
    m_gradTk[k] = cell->computeGradient("T", k);
  }
}

//***********************************************************************

void QAPConductivity::setGrad(const Coord &grad, int num)
{
  m_gradTk[num] = grad;
}

//***********************************************************************

Coord QAPConductivity::getGrad(int num) const
{
  return m_gradTk[num];
}

//***********************************************************************
