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
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      December 20 2017

#include "ElementCartesian.h"

using namespace std;

//***********************************************************************

ElementCartesian::ElementCartesian() : m_elementsChildren(0) {}

//***********************************************************************

ElementCartesian::~ElementCartesian()
{
  for (unsigned int i = 0; i < m_elementsChildren.size(); i++) {
    delete m_elementsChildren[i];
  }
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

double ElementCartesian::getSizeX()
{
  return m_size.getX();
}

//****************************************************************************

double ElementCartesian::getSizeY()
{
  return m_size.getY();
}

//****************************************************************************

double ElementCartesian::getSizeZ()
{
  return m_size.getZ();
}

//****************************************************************************

Coord ElementCartesian::getSize()
{
  return m_size;
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

void ElementCartesian::finalizeElementsChildren()
{
  for (unsigned int i = 0; i < m_elementsChildren.size(); i++) {
    delete m_elementsChildren[i];
  }
  m_elementsChildren.clear();
}
