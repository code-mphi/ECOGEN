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

//! \file      ElementPyramid.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include "ElementPyramid.h"

const int ElementPyramid::TYPEGMSH = 7;
const int ElementPyramid::NOMBRENOEUDS = 5;
const int ElementPyramid::NOMBREFACES = 5; /* ici il s'agit d'1 quadrangle et de 4 triangles*/
const int ElementPyramid::TYPEVTK = 14;

//***********************************************************************

ElementPyramid::ElementPyramid() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementPyramid::~ElementPyramid(){}

//***********************************************************************

void ElementPyramid::computeVolume(const Coord *noeuds)
{
  ////Aire du Quadrangle a la base
  ////----------------------------
  ////une diagonale :
  //Coord v0(noeuds[2] - noeuds[0]);
  //double diagonale(v0.norm());
  ////Les 4 cotes :
  //Coord v1(noeuds[1] - noeuds[0]);
  //Coord v2(noeuds[2] - noeuds[1]);
  //Coord v3(noeuds[3] - noeuds[2]);
  //Coord v4(noeuds[0] - noeuds[3]);
  //double a(v1.norm()); double b(v2.norm()); double c(v3.norm()); double d(v4.norm());
  ////Aire premier triangle
  //double dp = 0.5*(a + b + diagonale);
  //double surf1 = sqrt(dp*(dp - a)*(dp - b)*(dp - diagonale));
  ////Aire second triangle
  //dp = 0.5*(c + d + diagonale);
  //double surf2 = sqrt(dp*(dp - c)*(dp - d)*(dp - diagonale));
  ////Air quadrangle
  //double surfaceBase = surf1 + surf2;

  ////Hauteur
  ////-------
  //v1.setFromSubtractedVectors(noeuds[0], noeuds[1]);
  //v2.setFromSubtractedVectors(noeuds[0], noeuds[2]);
  //Coord v1v2;
  //v1v2 = Coord::crossProduct(v1, v2);
  //v1v2 = v1v2 / v1v2.norm();
  //Coord vh(noeuds[4] - noeuds[0]);
  //double hauteur = std::fabs(vh.scalar(v1v2));

  ////Volume de la pyramide a base quadrangle
  ////---------------------------------------
  //m_volume = (surfaceBase*hauteur) / 3.;

  //On va computeer le volume des 2 tetraedres inclus dans la pyramide
  Coord v1, v2, v3;
  v1.setFromSubtractedVectors(noeuds[4], noeuds[1]); v2.setFromSubtractedVectors(noeuds[4], noeuds[0]); v3.setFromSubtractedVectors(noeuds[4], noeuds[2]);
  double volumeT1 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  v1.setFromSubtractedVectors(noeuds[4], noeuds[2]); v2.setFromSubtractedVectors(noeuds[4], noeuds[3]); v3.setFromSubtractedVectors(noeuds[4], noeuds[0]);
  double volumeT2 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  m_volume = volumeT1 + volumeT2; //volume de la pyramide
  
}

//***********************************************************************

void ElementPyramid::computeLCFL(const Coord *noeuds)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((noeuds[0] + noeuds[1] + noeuds[2] + noeuds[3]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[0] + noeuds[1] + noeuds[4]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[1] + noeuds[2] + noeuds[4]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[2] + noeuds[3] + noeuds[4]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[3] + noeuds[0] + noeuds[4]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************

void ElementPyramid::construitFaces(const Coord *noeuds, FaceNS **faces, int &indexMaxFaces)
{
  //1 face de type quadrangle et 4 faces a traiter de type triangle
  int indexFaceExiste(0);
  int noeud1, noeud2, noeud3, noeud4, noeudAutre;

  //Construction de la base de la pyramide
  //--------------------------------------
  noeud1 = 0; noeud2 = 1; noeud3 = 2; noeud4 = 3; noeudAutre = 4;
  FaceQuadrangle face(m_numNoeuds[noeud1], m_numNoeuds[noeud2], m_numNoeuds[noeud3], m_numNoeuds[noeud4]);
  if (!face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
  {
    faces[indexMaxFaces] = new FaceQuadrangle(m_numNoeuds[noeud1], m_numNoeuds[noeud2], m_numNoeuds[noeud3], m_numNoeuds[noeud4]);
    faces[indexMaxFaces]->construitFace(noeuds, m_numNoeuds[noeudAutre], this);
    indexMaxFaces++;
  }
  else
  {
    faces[indexFaceExiste]->ajouteElementVoisin(this);
  }

  //Construction des 4 faces triangle
  //---------------------------------
  for (int i = 1; i < NOMBREFACES; i++)
  {
    switch (i)
    {
    case 1: noeud1 = 0; noeud2 = 1; noeud3 = 4; noeudAutre = 2; break;
    case 2: noeud1 = 1; noeud2 = 2; noeud3 = 4; noeudAutre = 3; break;
    case 3: noeud1 = 2; noeud2 = 3; noeud3 = 4; noeudAutre = 0; break;
    case 4: noeud1 = 3; noeud2 = 0; noeud3 = 4; noeudAutre = 1; break;
    }
    FaceTriangle face(m_numNoeuds[noeud1], m_numNoeuds[noeud2], m_numNoeuds[noeud3]);
    if (!face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexMaxFaces] = new FaceTriangle(m_numNoeuds[noeud1], m_numNoeuds[noeud2], m_numNoeuds[noeud3]);
      faces[indexMaxFaces]->construitFace(noeuds, m_numNoeuds[noeudAutre], this);
      indexMaxFaces++;
    }
    else
    {
      faces[indexFaceExiste]->ajouteElementVoisin(this);
    }
  }

}

//***********************************************************************