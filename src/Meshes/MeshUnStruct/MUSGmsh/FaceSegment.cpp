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

#include "FaceSegment.h"

const int FaceSegment::NUMBERNODES = 2;

//***********************************************************************

FaceSegment::FaceSegment(const int& numNode1, const int& numNode2, int tri) :
FaceNS(NUMBERNODES)
{
  m_numNodes[0] = numNode1;
  m_numNodes[1] = numNode2;
  if(tri) std::sort(m_numNodes, m_numNodes+2);
  m_sumNumNodes = m_numNodes[0] + m_numNodes[1];
}

//***********************************************************************

FaceSegment::~FaceSegment(){}

//***********************************************************************

void FaceSegment::computeSurface(const Coord* nodes)
{
  m_surface = (nodes[m_numNodes[1]] - nodes[m_numNodes[0]]).norm(); //Longeur du segment
}

//***********************************************************************

void FaceSegment::computeRepere(const Coord* nodes, const int& numNodeOther, ElementNS *elementNeighbor)
{
  Coord v1; v1.setFromSubtractedVectors(nodes[m_numNodes[0]], nodes[m_numNodes[1]]);
  m_tangent = v1 / v1.norm();
  m_binormal.setXYZ(0., 0., 1.);
  m_normal = Coord::crossProduct(m_tangent, m_binormal);

  Coord v2; v2.setFromSubtractedVectors(nodes[m_numNodes[0]], nodes[numNodeOther]);
  if (v2.scalar(m_normal) > 0.)
  {
    m_elementDroite = elementNeighbor;
    m_elementGauche = 0;
  }
  else
  {
    m_elementGauche = elementNeighbor;
    m_elementDroite = 0;
  }
}

//***********************************************************************