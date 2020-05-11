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

//! \file      ElementTetrahedron.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include "ElementTetrahedron.h"

const int ElementTetrahedron::TYPEGMSH = 4;
const int ElementTetrahedron::NOMBRENOEUDS = 4;
const int ElementTetrahedron::NOMBREFACES = 4; /* ici il s'agit de triangles*/
const int ElementTetrahedron::TYPEVTK = 10;

//***********************************************************************

ElementTetrahedron::ElementTetrahedron() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementTetrahedron::~ElementTetrahedron(){}

//***********************************************************************

void ElementTetrahedron::computeVolume(const Coord *noeuds)
{
  Coord v1(noeuds[1] - noeuds[0]), v2(noeuds[2] - noeuds[0]), v3(noeuds[3] - noeuds[0]);
  m_volume = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
}

//***********************************************************************

void ElementTetrahedron::computeLCFL(const Coord *noeuds)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((noeuds[0] + noeuds[1] + noeuds[2]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[1] + noeuds[2] + noeuds[3]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[2] + noeuds[3] + noeuds[0]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[3] + noeuds[0] + noeuds[1]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementTetrahedron::construitFaces(const Coord *noeuds, FaceNS **faces, int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //4 faces a traiter de type triangle
  int indexFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[2]; noeudAutre = 3; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[3]; noeudAutre = 0; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[3]; facesTemp[iMax][2] = m_numNoeuds[0]; noeudAutre = 1; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[0]; facesTemp[iMax][2] = m_numNoeuds[1]; noeudAutre = 2; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2];
    std::sort(facesTemp[iMax],facesTemp[iMax]+3);  //Tri des noeuds
    //Existance face ?
    indexFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,3,iMax);
    //Creation face ou rattachement
    if (indexFaceExiste==-1)
    {
      faces[iMax] = new FaceTriangle(facesTemp[iMax][0], facesTemp[iMax][1], facesTemp[iMax][2], 0); //pas besoin du tri ici
      faces[iMax]->construitFace(noeuds, m_numNoeuds[noeudAutre], this);
      iMax++;
    }
    else
    {
      faces[indexFaceExiste]->ajouteElementVoisin(this);
    }
  }
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementTetrahedron::construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //4 faces a traiter de type triangle
  int indexFaceExiste(-1);
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[2]; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[3]; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[3]; facesTemp[iMax][2] = m_numNoeuds[0]; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[0]; facesTemp[iMax][2] = m_numNoeuds[1]; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2];
    std::sort(facesTemp[iMax],facesTemp[iMax]+3);  //Tri des noeuds
    //Existance face ?
    indexFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,3,iMax);
    //Creation face ou rattachement
    if (indexFaceExiste==-1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementTetrahedron::attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indexMaxFaces, const int &numberNoeudsInternes)
{
  int indexFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[3]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[0] < numberNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[2], m_numNoeuds[3], m_numNoeuds[0]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[3], m_numNoeuds[0], m_numNoeuds[1]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementTetrahedron::compteFaceCommunicante(std::vector<int*> &facesTemp, std::vector<int> &sommeNoeudsTemp)
{
  //4 faces a traiter de type triangle
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[3], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[3]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; face[2] = m_numNoeuds[0]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; face[2] = m_numNoeuds[1]; break;     
    }
    int iMax = sommeNoeudsTemp.size();
    sommeNoeuds = face[0]+face[1]+face[2];
    std::sort(face, face+3);
    //Recherche existance faces
    indexFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,3,iMax);
    if (indexFaceExiste!=-1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementTetrahedron::compteFaceCommunicante(int &iMax, int **facesTemp, int *sommeNoeudsTemp)
{
  //4 faces a traiter de type triangle
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[3], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[3]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; face[2] = m_numNoeuds[0]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; face[2] = m_numNoeuds[1]; break;     
    }
    sommeNoeuds = face[0]+face[1]+face[2];
    std::sort(face, face+3);
    //Recherche existance faces
    indexFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,3,iMax);
    if (indexFaceExiste!=-1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************