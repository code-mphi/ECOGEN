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

//! \file      QAPViscosity.cpp
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "QAPViscosity.h"
#include <iostream>

using namespace std;

//***********************************************************************

QAPViscosity::QAPViscosity(){}

//***********************************************************************

QAPViscosity::QAPViscosity(AddPhys* addPhys) : QuantitiesAddPhys(addPhys), m_gradU(0.), m_gradV(0.), m_gradW(0.)
{}

//***********************************************************************

QAPViscosity::~QAPViscosity(){}

//***********************************************************************

void QAPViscosity::computeQuantities(Cell* cell)
{
  m_gradU = cell->computeGradient("u");
  m_gradV = cell->computeGradient("v");
  m_gradW = cell->computeGradient("w");
}

//***********************************************************************

void QAPViscosity::setGrad(const Coord &grad, int num)
{
  switch (num) {
  case 1: m_gradU = grad; break;
  case 2: m_gradV = grad; break;
  case 3: m_gradW = grad; break;
  default: Errors::errorMessage("Error in QAPViscosity::setGrad value of num non defined"); break;
  }
}

//***********************************************************************

Coord QAPViscosity::getGrad(int num) const
{
  switch (num) {
  case 1: return m_gradU; break;
  case 2: return m_gradV; break;
  case 3: return m_gradW; break;
  default: Errors::errorMessage("Error in QAPViscosity::getGrad value of num non defined"); break;
  }
  return 0;
}

//***********************************************************************