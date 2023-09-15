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

#include "ElementPoint.h"

const int ElementPoint::TYPEGMSH = 15;
const int ElementPoint::NUMBERNODES = 1;
const int ElementPoint::NUMBERFACES = 0;
const int ElementPoint::TYPEVTK = 1;

//***********************************************************************

ElementPoint::ElementPoint() :
ElementNS(TYPEGMSH, NUMBERNODES, NUMBERFACES, TYPEVTK)
{}

//***********************************************************************

ElementPoint::~ElementPoint(){}

//***********************************************************************

void ElementPoint::computeVolume(const Coord* /*nodes*/)
{
  m_volume = 1.0; // Without unit, useless
}

//***********************************************************************

void ElementPoint::computeLCFL(const Coord* /*nodes*/)
{
  m_lCFL = 1.0; // Useless
}

//***********************************************************************

void ElementPoint::attributFaceLimite(FaceNS** faces, const int& indexMaxFaces)
{
  int indexFaceExiste(0);
  FacePoint face(m_numNodes[0]);
  if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
  {
    faces[indexFaceExiste]->addElementNeighborLimite(this);
  }
  else
  {
    Errors::errorMessage("Probleme attribution des faces limites element Segment");
  }
}

//**************************************