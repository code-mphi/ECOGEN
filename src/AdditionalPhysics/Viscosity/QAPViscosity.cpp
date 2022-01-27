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

#include "QAPViscosity.h"

//***********************************************************************

QAPViscosity::QAPViscosity(AddPhys* addPhys) : QuantitiesAddPhys(addPhys), m_grads(3)
{
  variableNamesVisc.resize(3);
  numPhasesVisc.resize(3);
  for (int i = 0; i < 3; ++i) {
    m_grads[i] = 0.;
    numPhasesVisc[i] = 0;
  }
  variableNamesVisc[0] = velocityU;
  variableNamesVisc[1] = velocityV;
  variableNamesVisc[2] = velocityW;
}

//***********************************************************************

QAPViscosity::~QAPViscosity(){}

//***********************************************************************

void QAPViscosity::computeQuantities(Cell* cell)
{
  cell->computeGradient(m_grads, variableNamesVisc, numPhasesVisc);
}

//***********************************************************************

void QAPViscosity::setGrad(const Coord& grad, const int& num)
{
  m_grads[num-1] = grad;
}

//***********************************************************************