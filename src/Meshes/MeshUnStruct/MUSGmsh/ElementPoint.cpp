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

//! \file      ElementPoint.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include "ElementPoint.h"

const int ElementPoint::TYPEGMSH = 15;
const int ElementPoint::NOMBRENOEUDS = 1;
const int ElementPoint::NOMBREFACES = 0;
const int ElementPoint::TYPEVTK = 1;

//***********************************************************************

ElementPoint::ElementPoint() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementPoint::~ElementPoint(){}

//***********************************************************************

void ElementPoint::computeVolume(const Coord *noeuds)
{
  m_volume = 1.0; //sans unite inutile
}

//***********************************************************************

void ElementPoint::computeLCFL(const Coord *noeuds)
{
  m_lCFL = 1.0; //inutile
}

//***********************************************************************

void ElementPoint::attributFaceLimite(const Coord *noeuds, FaceNS **faces, const int &indexMaxFaces)
{
  int indexFaceExiste(0);
  FacePoint face(m_numNoeuds[0]);
  if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
  {
    faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
  }
  else
  {
    Errors::errorMessage("Probleme attribution des faces limites element Segment");
  }
}

//**************************************