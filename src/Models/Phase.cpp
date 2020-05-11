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

//! \file      Phase.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      December 21 2017

#include "Phase.h"

//***************************************************************************

Phase::Phase(){}

//***************************************************************************

Phase::~Phase(){}

//***************************************************************************

void Phase::printPhase(std::ofstream &fileStream) const
{
  //Scalar variables
  for (int var = 1; var <= this->getNumberScalars(); var++) {
    fileStream << this->returnScalar(var) << " ";
  }
  //Vector variables
  for (int var = 1; var <= this->getNumberVectors(); var++) {
    fileStream << this->returnVector(var).norm() << " ";
  }
}

//***************************************************************************
