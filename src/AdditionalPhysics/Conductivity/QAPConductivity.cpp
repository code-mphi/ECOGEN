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

#include "QAPConductivity.h"

//***********************************************************************

QAPConductivity::QAPConductivity(AddPhys* addPhys) : QuantitiesAddPhys(addPhys), m_gradTk(numberPhases)
{
  variableNamesCond.resize(numberPhases);
  numPhasesCond.resize(numberPhases);
  for (int k = 0; k < numberPhases; ++k) {
    m_gradTk[k] = 0.;
    variableNamesCond[k] = temperature;
    numPhasesCond[k] = k;
  }
}

//***********************************************************************

QAPConductivity::~QAPConductivity(){}

//***********************************************************************

void QAPConductivity::computeQuantities(Cell* cell)
{
  cell->computeGradient(m_gradTk, variableNamesCond, numPhasesCond);
}

//***********************************************************************

void QAPConductivity::setGrad(const Coord& grad, const int& num)
{
  m_gradTk[num] = grad;
}

//***********************************************************************
