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

//! \file      GOLine.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      January 5 2018

#include "GOLine.h"
#include <iostream>

//***********************************************************************

GOLine::GOLine(){}

//***********************************************************************

GOLine::GOLine(const Coord &vertex, const Coord &vecDir) :
  GeometricObject(LINE), m_point(vertex), m_vecDir(vecDir)
{
  if (vecDir.norm() < 1e-6) {
    Errors::errorMessage("GOLine::GOLine impossible to create line, director vector null");
  }
  else {
    m_vecDir.normalized();
  }
}

//***********************************************************************

GOLine::~GOLine(){}

//***********************************************************************

double GOLine::distancePoint(const Coord &vertex) const
{
  Coord vec; vec.setFromSubtractedVectors(m_point, vertex);
  if (std::fabs(vec.scalar(m_vecDir)) < 1e-6) { return 0.; }
  else { return  (vec.cross(m_vecDir)).norm(); }
}

//***********************************************************************

Coord GOLine::projectionPoint(const Coord &vertex) const
{
  Coord projection;
  Coord vec; vec.setFromSubtractedVectors(m_point, vertex);
  projection.setXYZ(m_vecDir.scalar(vec),0.,0.);
  return projection;
}

//***********************************************************************