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

//! \file      ElementQuadrangle.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include "ElementQuadrangle.h"

const int ElementQuadrangle::TYPEGMSH = 3;
const int ElementQuadrangle::NOMBRENOEUDS = 4;
const int ElementQuadrangle::NOMBREFACES = 4; /* ici il s'agit du number de segments*/
const int ElementQuadrangle::TYPEVTK = 9;

//***********************************************************************

ElementQuadrangle::ElementQuadrangle() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementQuadrangle::~ElementQuadrangle(){}

//***********************************************************************

void ElementQuadrangle::computeVolume(const Coord *noeuds)
{
  //une diagonale :
  Coord v0(noeuds[2] - noeuds[0]);
  double diagonale(v0.norm());
  //Les 4 cotes :
  Coord v1(noeuds[1] - noeuds[0]);
  Coord v2(noeuds[2] - noeuds[1]);
  Coord v3(noeuds[3] - noeuds[2]);
  Coord v4(noeuds[0] - noeuds[3]);
  double a(v1.norm()); double b(v2.norm()); double c(v3.norm()); double d(v4.norm());
  //Aire premier triangle
  double dp1 = 0.5*(a + b + diagonale);
  double surf1 = sqrt(dp1*(dp1 - a)*(dp1 - b)*(dp1 - diagonale));
  //Aire second triangle
  double dp2 = 0.5*(c + d + diagonale);
  double surf2 = sqrt(dp2*(dp2 - c)*(dp2 - d)*(dp2 - diagonale));
  //Air quadrangle
  m_volume = surf1 + surf2; //aire du quadrangle
}

//***********************************************************************

void ElementQuadrangle::computeLCFL(const Coord *noeuds)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((noeuds[0] + noeuds[1]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[1] + noeuds[2]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[2] + noeuds[3]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[3] + noeuds[0]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementQuadrangle::construitFaces(const Coord *noeuds, FaceNS **faces, int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //4 faces a traiter de type segment
  int indexFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; noeudAutre = 2; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; noeudAutre = 3; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[3]; noeudAutre = 0; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[0]; noeudAutre = 1; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1];
    std::sort(facesTemp[iMax],facesTemp[iMax]+2);  //Tri des noeuds
    //Existance face ?
    indexFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,2,iMax);
    //Creation face ou rattachement
    if (indexFaceExiste==-1)
    {
      faces[iMax] = new FaceSegment(facesTemp[iMax][0], facesTemp[iMax][1], 0); //pas besoin du tri ici
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
void ElementQuadrangle::construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //4 faces a traiter de type segment
  int indexFaceExiste(-1);
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[3]; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[0]; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1];
    std::sort(facesTemp[iMax],facesTemp[iMax]+2);  //Tri des noeuds
    //Existance face ?
    indexFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,2,iMax);
    //Creation face ou rattachement
    if (indexFaceExiste==-1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementQuadrangle::attributFaceLimite(const Coord *noeuds, FaceNS **faces, const int &indexMaxFaces)
{
  int indexFaceExiste(0);
  FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[3]);
  if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
  {
    faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
  }
  else
  {
    Errors::errorMessage("Probleme attribution des faces limites element Quadrangle");
  }
}

//***********************************************************************

void ElementQuadrangle::attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indexMaxFaces, const int &numberNoeudsInternes)
{
  int indexFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[0], m_numNoeuds[1]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[1], m_numNoeuds[2]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[2], m_numNoeuds[3]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[0] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[3], m_numNoeuds[0]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementQuadrangle::compteFaceCommunicante(std::vector<int*> &facesTemp, std::vector<int> &sommeNoeudsTemp)
{
  //4 faces a traiter de type segment
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[2], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; break;     
    }
    int iMax = sommeNoeudsTemp.size();
    sommeNoeuds = face[0]+face[1];
    std::sort(face, face+2);
    //Recherche existance faces
    indexFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,2,iMax);
    if (indexFaceExiste!=-1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementQuadrangle::compteFaceCommunicante(int &iMax, int **facesTemp, int *sommeNoeudsTemp)
{
  //4 faces a traiter de type segment
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[2], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; break;     
    }
    sommeNoeuds = face[0]+face[1];
    std::sort(face, face+2);
    //Recherche existance faces
    indexFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,2,iMax);
    if (indexFaceExiste!=-1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************