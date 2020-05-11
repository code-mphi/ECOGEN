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

//! \file      Element.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include "Element.h"

//***********************************************************************

Element::Element() : m_position(0), m_volume(0.), m_lCFL(0.), m_numCellAssociee(0) {}

//***********************************************************************

Element::~Element(){}

//***********************************************************************

void Element::ecritPos(std::ofstream &fileStream, Axis axis)
{
  switch (axis) {
  case X:
    fileStream << m_position.getX() << " "; break;
  case Y:
    fileStream << m_position.getY() << " "; break;
  case Z:
    fileStream << m_position.getZ() << " "; break;
  default:
    Errors::errorMessage("Element::ecritPos : Axis unknown"); break;
  }
}

//***********************************************************************

Coord Element::vecteur(const Element *e)
{
  Coord vec;
  vec.setFromSubtractedVectors(m_position, e->m_position);
  return vec;
}

//***********************************************************************

Coord Element::vecteur(const Face *f)
{
  Coord vec;
  vec.setFromSubtractedVectors(m_position, f->getPos());
  return vec;
}

//***********************************************************************

double Element::distance(const Element *e) 
{
  Coord vec(m_position - e->m_position);
  return vec.norm();
}

//***********************************************************************

double Element::distanceX(const Element *e)
{
	Coord vec(m_position - e->m_position);
	return vec.getX();
}

//***********************************************************************

double Element::distanceY(const Element *e)
{
	Coord vec(m_position - e->m_position);
	return vec.getY();
}

//***********************************************************************

double Element::distanceZ(const Element *e)
{
	Coord vec(m_position - e->m_position);
	return vec.getZ();
}

//***********************************************************************

double Element::distance(const Face *f)
{
  Coord vec(m_position - f->getPos());
  return vec.norm();
}

//***********************************************************************

double Element::distanceX(const Face *f)
{
  Coord vec(m_position - f->getPos());
  return vec.getX();
}

//***********************************************************************

double Element::distanceY(const Face *f)
{
  Coord vec(m_position - f->getPos());
  return vec.getY();
}

//***********************************************************************

double Element::distanceZ(const Face *f)
{
  Coord vec(m_position - f->getPos());
  return vec.getZ();
}

//***********************************************************************

bool Element::traverseObjet(const GeometricObject &objet) const
{
  if (objet.distancePoint(m_position) < m_lCFL) return true;
  return false;
}

//***********************************************************************

//***********************************************************************
//*********************** Parallel load balancing ***********************
//***********************************************************************

void Element::setKey(const decomposition::Key<3> &key)
{
  m_key = key;
}

//***********************************************************************