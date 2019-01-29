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

//! \file      QAPSurfaceTension.cpp
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "QAPSurfaceTension.h"
#include <iostream>

using namespace std;

//***********************************************************************

QAPSurfaceTension::QAPSurfaceTension(){}


//***********************************************************************

QAPSurfaceTension::QAPSurfaceTension(AddPhys* addPhys) : QuantitiesAddPhys(addPhys), m_gradC(0.,0.,0.)
{}

//***********************************************************************

QAPSurfaceTension::~QAPSurfaceTension(){}

//***********************************************************************

void QAPSurfaceTension::computeQuantities(Cell* cell)
{
  m_gradC = cell->computeGradient("TR", m_addPhys->getNumTransportAssociated());
}

//***********************************************************************

void QAPSurfaceTension::setGrad(const Coord &grad, int num)
{
  m_gradC = grad;
}

//***********************************************************************

Coord QAPSurfaceTension::getGrad(int num) const
{
  return m_gradC;
}

//***********************************************************************

