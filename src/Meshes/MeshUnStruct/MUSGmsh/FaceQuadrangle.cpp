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

//! \file      FaceQuadrangle.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include "FaceQuadrangle.h"

const int FaceQuadrangle::NOMBRENOEUDS = 4;

//***********************************************************************

FaceQuadrangle::FaceQuadrangle(const int &numNoeud1, const int &numNoeud2, const int &numNoeud3, const int &numNoeud4, int tri) :
FaceNS(NOMBRENOEUDS)
{
  //Sauvegarde ordre initial
  m_numNoeudsOrigine = new int[4];
  m_numNoeudsOrigine[0] = numNoeud1;
  m_numNoeudsOrigine[1] = numNoeud2;
  m_numNoeudsOrigine[2] = numNoeud3;
  m_numNoeudsOrigine[3] = numNoeud4;

  m_numNoeuds[0] = numNoeud1;
  m_numNoeuds[1] = numNoeud2;
  m_numNoeuds[2] = numNoeud3;
  m_numNoeuds[3] = numNoeud4;
  if(tri) std::sort(m_numNoeuds, m_numNoeuds+4);
  m_sommeNumNoeuds = m_numNoeuds[0] + m_numNoeuds[1] + m_numNoeuds[2] + m_numNoeuds[3];
}

//***********************************************************************

FaceQuadrangle::~FaceQuadrangle()
{
  delete[] m_numNoeudsOrigine;
}

//***********************************************************************

void FaceQuadrangle::computeSurface(const Coord *noeuds)
{
  //Atention utilisation des numbering de noeud d origin pour assurer le compute des surfaces
  //une diagonale :
  Coord v0(noeuds[m_numNoeudsOrigine[2]] - noeuds[m_numNoeudsOrigine[0]]);
  double diagonale(v0.norm());
  //Les 4 cotes :
  Coord v1(noeuds[m_numNoeudsOrigine[1]] - noeuds[m_numNoeudsOrigine[0]]);
  Coord v2(noeuds[m_numNoeudsOrigine[2]] - noeuds[m_numNoeudsOrigine[1]]);
  Coord v3(noeuds[m_numNoeudsOrigine[3]] - noeuds[m_numNoeudsOrigine[2]]);
  Coord v4(noeuds[m_numNoeudsOrigine[0]] - noeuds[m_numNoeudsOrigine[3]]);
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

void FaceQuadrangle::computeRepere(const Coord *noeuds, const int &numNoeudAutre, ElementNS *elementVoisin)
{
  Coord v1; v1.setFromSubtractedVectors(noeuds[m_numNoeuds[0]], noeuds[m_numNoeuds[1]]);
  Coord v2; v2.setFromSubtractedVectors(noeuds[m_numNoeuds[0]], noeuds[m_numNoeuds[2]]);

  m_tangent = v1 / v1.norm();
  Coord v1v2; v1v2 = Coord::crossProduct(v1, v2);
  m_normal = v1v2 / v1v2.norm();
  m_binormal = Coord::crossProduct(m_normal, m_tangent);

  Coord v3; v3.setFromSubtractedVectors(noeuds[m_numNoeuds[0]], noeuds[numNoeudAutre]);
  if (v3.scalar(m_normal) > 0.)
  {
    m_elementDroite = elementVoisin;
    m_elementGauche = 0;
  }
  else
  {
    m_elementGauche = elementVoisin;
    m_elementDroite = 0;
  }
}

//***********************************************************************