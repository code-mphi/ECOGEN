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

//! \file      Face.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      December 20 2017

#include "Face.h"

using namespace std;

//***********************************************************************

Face::Face(){}

//***********************************************************************

Face::~Face(){}

//***********************************************************************

Coord Face::getNormal() const
{
  return m_normal;
}

//***********************************************************************

Coord Face::getTangent() const
{
  return m_tangent;
}

//***********************************************************************

Coord Face::getBinormal() const
{
  return m_binormal;
}

//***********************************************************************

double Face::getSurface() const
{
  return m_surface;
}

//***********************************************************************

Coord Face::getPos() const
{
  return m_position;
}

//***********************************************************************

Coord Face::vecteur(Element *e)
{
  Coord vec;
  vec.setFromSubtractedVectors(m_position, e->getPosition());
  return vec;
}

//***********************************************************************

double Face::distance(Element *e)
{
  Coord vec(m_position - e->getPosition());
  return vec.norm();
}

//***********************************************************************