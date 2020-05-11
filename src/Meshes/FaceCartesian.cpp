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

//! \file      FaceCartesian.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "FaceCartesian.h"

//***********************************************************************

FaceCartesian::FaceCartesian(){}

//***********************************************************************

FaceCartesian::~FaceCartesian(){}

//***********************************************************************

void FaceCartesian::setSurface(const double &surface)
{
  m_surface = surface;
}

//***********************************************************************

void FaceCartesian::initializeAutres(const double &surface, const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  m_surface = surface;
  m_normal = normal;
  m_tangent = tangent;
  m_binormal = binormal;
}

//***********************************************************************

void FaceCartesian::setPos(const double &X, const double &Y, const double &Z)
{
  m_position.setXYZ(X, Y, Z);
}

//***********************************************************************

void FaceCartesian::setNormal(const double &X, const double &Y, const double &Z)
{
  m_normal.setXYZ(X, Y, Z);
}

//***********************************************************************

void FaceCartesian::setTangent(const double &X, const double &Y, const double &Z)
{
  m_tangent.setXYZ(X, Y, Z);
}

//***********************************************************************

void FaceCartesian::setBinormal(const double &X, const double &Y, const double &Z)
{
  m_binormal.setXYZ(X, Y, Z);
}

//****************************************************************************

void FaceCartesian::setSize(const double &sizeX, const double &sizeY, const double &sizeZ)
{
  m_size.setXYZ(sizeX, sizeY, sizeZ);
}

//****************************************************************************

void FaceCartesian::setSize(const Coord &size)
{
  m_size = size;
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

Face* FaceCartesian::creerNouvelleFace()
{
  return new FaceCartesian;
}

