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

#include "FacePoint.h"

const int FacePoint::NUMBERNODES = 1;

//***********************************************************************

FacePoint::FacePoint(const int& numNode1) :
FaceNS(NUMBERNODES)
{
  m_numNodes[0] = numNode1;
  m_sumNumNodes = m_numNodes[0];
}

//***********************************************************************

FacePoint::~FacePoint(){}

//***********************************************************************

void FacePoint::computeSurface(const Coord* /*nodes*/)
{
  m_surface = 1.0; //unite
}

//***********************************************************************

void FacePoint::computeRepere(const Coord* nodes, const int& numNodeOther, ElementNS *elementNeighbor)
{
  Coord v1; v1.setFromSubtractedVectors(nodes[m_numNodes[0]], nodes[numNodeOther]);
  m_normal = v1 / v1.norm();
  m_binormal.setXYZ(0., 0., 1.);
  m_tangent = Coord::crossProduct(m_binormal, m_normal);

  m_elementDroite = elementNeighbor;
  m_elementGauche = 0;
}

//***********************************************************************