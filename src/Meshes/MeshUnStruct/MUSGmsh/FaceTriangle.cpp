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

#include "FaceTriangle.h"

const int FaceTriangle::NUMBERNODES=3;

//***********************************************************************

FaceTriangle::FaceTriangle(const int& numNode1, const int& numNode2, const int& numNode3, int tri) :
FaceNS(NUMBERNODES)
{
  m_numNodes[0] = numNode1;
  m_numNodes[1] = numNode2;
  m_numNodes[2] = numNode3;
  if(tri) std::sort(m_numNodes, m_numNodes+3);
  m_sumNumNodes = m_numNodes[0] + m_numNodes[1] + m_numNodes[2];
}

//***********************************************************************

FaceTriangle::~FaceTriangle(){}

//***********************************************************************

void FaceTriangle::computeSurface(const Coord* nodes)
{
  Coord v1(nodes[m_numNodes[1]] - nodes[m_numNodes[0]]);
  Coord v2(nodes[m_numNodes[2]] - nodes[m_numNodes[1]]);
  Coord v3(nodes[m_numNodes[0]] - nodes[m_numNodes[2]]);
  double a(v1.norm()), b(v2.norm()), c(v3.norm());
  double dp = 0.5*(a + b + c);
  m_surface = sqrt(dp*(dp - a)*(dp - b)*(dp - c));
}

//***********************************************************************

void FaceTriangle::computeRepere(const Coord* nodes, const int& numNodeOther, ElementNS *elementNeighbor)
{
  Coord v1; v1.setFromSubtractedVectors(nodes[m_numNodes[0]], nodes[m_numNodes[1]]);
  Coord v2; v2.setFromSubtractedVectors(nodes[m_numNodes[0]], nodes[m_numNodes[2]]);

  m_tangent = v1 / v1.norm();
  Coord v1v2; v1v2 = Coord::crossProduct(v1, v2);
  m_normal = v1v2 / v1v2.norm();
  m_binormal = Coord::crossProduct(m_normal, m_tangent);
  
  Coord v3; v3.setFromSubtractedVectors(nodes[m_numNodes[0]], nodes[numNodeOther]);
  if (v3.scalar(m_normal) > 0.)
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