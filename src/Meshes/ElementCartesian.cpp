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

//! \file      ElementCartesian.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include "ElementCartesian.h"

//***********************************************************************

ElementCartesian::ElementCartesian() : m_elementsChildren(0) {}

//***********************************************************************

ElementCartesian::~ElementCartesian()
{
  m_elementsChildren.clear();
}

//***********************************************************************

void ElementCartesian::setVolume(const double &volume)
{
  m_volume = volume;
}

//***********************************************************************

void ElementCartesian::setLCFL(const double &lCFL)
{
  m_lCFL = lCFL;
}

//***********************************************************************

void ElementCartesian::setPos(const double &X, const double &Y, const double &Z)
{
  m_position.setXYZ(X, Y, Z);
}

//***********************************************************************


void ElementCartesian::setPos(const Coord &pos)
{
  m_position = pos;
}

//***********************************************************************

void ElementCartesian::setPosX(const double &X)
{
  m_position.setX(X);
}

//***********************************************************************

void ElementCartesian::setPosY(const double &Y)
{
  m_position.setY(Y);
}

//***********************************************************************

void ElementCartesian::setPosZ(const double &Z)
{
  m_position.setZ(Z);
}

//****************************************************************************

void ElementCartesian::setSize(const double &sizeX, const double &sizeY, const double &sizeZ)
{
  m_size.setXYZ(sizeX, sizeY, sizeZ);
}

//****************************************************************************

void ElementCartesian::setSize(const Coord &size)
{
  m_size = size;
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

void ElementCartesian::creerElementChild()
{
  m_elementsChildren.push_back(new ElementCartesian);
}

//****************************************************************************

Element* ElementCartesian::getElementChild(const int &numberChild)
{
  return m_elementsChildren[numberChild];
}

//****************************************************************************

Element* ElementCartesian::getElementChildBack()
{
  return m_elementsChildren.back();
}
//****************************************************************************

void ElementCartesian::finalizeElementsChildren()
{
  m_elementsChildren.clear();
}
