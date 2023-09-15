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

#include "FaceQuadrangle.h"

const int FaceQuadrangle::NUMBERNODES = 4;

//***********************************************************************

FaceQuadrangle::FaceQuadrangle(const int& numNode1, const int& numNode2, const int& numNode3, const int& numNode4, int tri) :
FaceNS(NUMBERNODES)
{
  //Sauvegarde ordre initial
  m_numNodesOrigine = new int[4];
  m_numNodesOrigine[0] = numNode1;
  m_numNodesOrigine[1] = numNode2;
  m_numNodesOrigine[2] = numNode3;
  m_numNodesOrigine[3] = numNode4;

  m_numNodes[0] = numNode1;
  m_numNodes[1] = numNode2;
  m_numNodes[2] = numNode3;
  m_numNodes[3] = numNode4;
  if(tri) std::sort(m_numNodes, m_numNodes+4);
  m_sumNumNodes = m_numNodes[0] + m_numNodes[1] + m_numNodes[2] + m_numNodes[3];
}

//***********************************************************************

FaceQuadrangle::~FaceQuadrangle()
{
  delete[] m_numNodesOrigine;
}

//***********************************************************************

void FaceQuadrangle::computeSurface(const Coord* nodes)
{
  //Atention utilisation des numbering de node d origin pour assurer le compute des surfaces
  //une diagonale :
  Coord v0(nodes[m_numNodesOrigine[2]] - nodes[m_numNodesOrigine[0]]);
  double diagonale(v0.norm());
  //Les 4 cotes :
  Coord v1(nodes[m_numNodesOrigine[1]] - nodes[m_numNodesOrigine[0]]);
  Coord v2(nodes[m_numNodesOrigine[2]] - nodes[m_numNodesOrigine[1]]);
  Coord v3(nodes[m_numNodesOrigine[3]] - nodes[m_numNodesOrigine[2]]);
  Coord v4(nodes[m_numNodesOrigine[0]] - nodes[m_numNodesOrigine[3]]);
  double a(v1.norm()); double b(v2.norm()); double c(v3.norm()); double d(v4.norm());
  //Aire premier triangle
  double dp1 = 0.5*(a + b + diagonale); 
  double surf1 = sqrt(dp1*(dp1 - a)*(dp1 - b)*(dp1 - diagonale));
  //Aire second triangle
  double dp2 = 0.5*(c + d + diagonale);
  double surf2 = sqrt(dp2*(dp2 - c)*(dp2 - d)*(dp2 - diagonale));
  //Air quadrangle
  m_surface = surf1 + surf2;
}

//***********************************************************************

void FaceQuadrangle::computeRepere(const Coord* nodes, const int& numNodeOther, ElementNS *elementNeighbor)
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