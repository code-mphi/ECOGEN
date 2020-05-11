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

//! \file      ElementHexahedron.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include "ElementHexahedron.h"

const int ElementHexahedron::TYPEGMSH = 5;
const int ElementHexahedron::NOMBRENOEUDS = 8;
const int ElementHexahedron::NOMBREFACES = 6; /* ici il s'agit de quadrangles*/
const int ElementHexahedron::TYPEVTK = 12;

//***********************************************************************

ElementHexahedron::ElementHexahedron() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementHexahedron::~ElementHexahedron(){}

//***********************************************************************

void ElementHexahedron::computeVolume(const Coord *noeuds)
{
  //On va computeer le volume des 6 tetraedres inclus dans l hexaerdre
  Coord v1, v2, v3;
  v1.setFromSubtractedVectors(noeuds[0], noeuds[1]); v2.setFromSubtractedVectors(noeuds[0], noeuds[3]); v3.setFromSubtractedVectors(noeuds[0], noeuds[4]);
  double volumeT1 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  v1.setFromSubtractedVectors(noeuds[3], noeuds[1]); v2.setFromSubtractedVectors(noeuds[3], noeuds[7]); v3.setFromSubtractedVectors(noeuds[3], noeuds[4]);
  double volumeT2 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  v1.setFromSubtractedVectors(noeuds[3], noeuds[1]); v2.setFromSubtractedVectors(noeuds[3], noeuds[7]); v3.setFromSubtractedVectors(noeuds[3], noeuds[2]);
  double volumeT3 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  v1.setFromSubtractedVectors(noeuds[4], noeuds[1]); v2.setFromSubtractedVectors(noeuds[4], noeuds[7]); v3.setFromSubtractedVectors(noeuds[4], noeuds[5]);
  double volumeT4 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  v1.setFromSubtractedVectors(noeuds[7], noeuds[2]); v2.setFromSubtractedVectors(noeuds[7], noeuds[6]); v3.setFromSubtractedVectors(noeuds[7], noeuds[5]);
  double volumeT5 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  v1.setFromSubtractedVectors(noeuds[7], noeuds[2]); v2.setFromSubtractedVectors(noeuds[7], noeuds[5]); v3.setFromSubtractedVectors(noeuds[7], noeuds[1]);
  double volumeT6 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  m_volume = volumeT1 + volumeT2 + volumeT3 + volumeT4 + volumeT5 + volumeT6; //volume de l Hexahedron
}

//***********************************************************************

void ElementHexahedron::computeLCFL(const Coord *noeuds)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((noeuds[0] + noeuds[1] + noeuds[2] + noeuds[3]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[4] + noeuds[5] + noeuds[6] + noeuds[7]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[0] + noeuds[4] + noeuds[7] + noeuds[3]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[1] + noeuds[5] + noeuds[6] + noeuds[2]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[0] + noeuds[1] + noeuds[5] + noeuds[4]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((noeuds[3] + noeuds[2] + noeuds[6] + noeuds[7]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementHexahedron::construitFaces(const Coord *noeuds, FaceNS **faces, int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //6 faces a traiter de type quadrangle
  int indexFaceExiste(-1);
  int noeudAutre;  
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[2]; facesTemp[iMax][3] = m_numNoeuds[3]; noeudAutre = 4; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[4]; facesTemp[iMax][1] = m_numNoeuds[5]; facesTemp[iMax][2] = m_numNoeuds[6]; facesTemp[iMax][3] = m_numNoeuds[7]; noeudAutre = 0; break;
      case 2: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[3]; facesTemp[iMax][2] = m_numNoeuds[7]; facesTemp[iMax][3] = m_numNoeuds[4]; noeudAutre = 1; break;
      case 3: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[6]; facesTemp[iMax][3] = m_numNoeuds[5]; noeudAutre = 0; break;
      case 4: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[5]; facesTemp[iMax][3] = m_numNoeuds[4]; noeudAutre = 2; break;
      case 5: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[6]; facesTemp[iMax][3] = m_numNoeuds[7]; noeudAutre = 0; break;
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2] + facesTemp[iMax][3];
    std::sort(facesTemp[iMax],facesTemp[iMax]+4);  //Tri des noeuds
    //Existance face ?
    indexFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,4,iMax);
    //Creation face ou rattachement
    if (indexFaceExiste==-1) // on fill simultanement le tableau faces et le tableau facesTemp
    {
      faces[iMax] = new FaceQuadrangle(facesTemp[iMax][0], facesTemp[iMax][1], facesTemp[iMax][2], facesTemp[iMax][3], 0); //pas besoin du tri ici
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
void ElementHexahedron::construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //6 faces a traiter de type quadrangle
  int indexFaceExiste(-1);
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
    case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[2]; facesTemp[iMax][3] = m_numNoeuds[3]; break;
    case 1: facesTemp[iMax][0] = m_numNoeuds[4]; facesTemp[iMax][1] = m_numNoeuds[5]; facesTemp[iMax][2] = m_numNoeuds[6]; facesTemp[iMax][3] = m_numNoeuds[7]; break;
    case 2: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[3]; facesTemp[iMax][2] = m_numNoeuds[7]; facesTemp[iMax][3] = m_numNoeuds[4]; break;
    case 3: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[6]; facesTemp[iMax][3] = m_numNoeuds[5]; break;
    case 4: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[5]; facesTemp[iMax][3] = m_numNoeuds[4]; break;
    case 5: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[6]; facesTemp[iMax][3] = m_numNoeuds[7]; break;
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2] + facesTemp[iMax][3];
    std::sort(facesTemp[iMax], facesTemp[iMax] + 4);  //Tri des noeuds
    //Existance face ?
    indexFaceExiste = FaceNS::rechercheFace(facesTemp[iMax], sommeNoeudsTemp[iMax], facesTemp, sommeNoeudsTemp, 4, iMax);
    //Creation face ou rattachement
    if (indexFaceExiste == -1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementHexahedron::attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indexMaxFaces, const int &numberNoeudsInternes)
{
  int indexFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[3]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[4] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes && m_numNoeuds[6] < numberNoeudsInternes && m_numNoeuds[7] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[4], m_numNoeuds[5], m_numNoeuds[6], m_numNoeuds[7]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[7] < numberNoeudsInternes && m_numNoeuds[4] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[3], m_numNoeuds[7], m_numNoeuds[4]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[6] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[6], m_numNoeuds[5]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 5 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes && m_numNoeuds[4] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[5], m_numNoeuds[4]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 6 :
  if (m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[6] < numberNoeudsInternes && m_numNoeuds[7] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[3], m_numNoeuds[2], m_numNoeuds[6], m_numNoeuds[7]);
    if (face.faceExiste(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementHexahedron::compteFaceCommunicante(std::vector<int*> &facesTemp, std::vector<int> &sommeNoeudsTemp)
{
  //6 faces a traiter de type quadrangle
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[4], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
    case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; face[3] = m_numNoeuds[3]; break;
    case 1: face[0] = m_numNoeuds[4]; face[1] = m_numNoeuds[5]; face[2] = m_numNoeuds[6]; face[3] = m_numNoeuds[7]; break;
    case 2: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[3]; face[2] = m_numNoeuds[7]; face[3] = m_numNoeuds[4]; break;
    case 3: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[6]; face[3] = m_numNoeuds[5]; break;
    case 4: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[4]; break;
    case 5: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[6]; face[3] = m_numNoeuds[7]; break;
    }
    int iMax = sommeNoeudsTemp.size();
    sommeNoeuds = face[0] + face[1] + face[2] + face[3];
    std::sort(face, face + 4);
    //Recherche existance faces
    indexFaceExiste = FaceNS::rechercheFace(face, sommeNoeuds, facesTemp, sommeNoeudsTemp, 4, iMax);
    if (indexFaceExiste != -1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementHexahedron::compteFaceCommunicante(int &iMax, int **facesTemp, int *sommeNoeudsTemp)
{
  //6 faces a traiter de type quadrangle
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[4], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
    case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; face[3] = m_numNoeuds[3]; break;
    case 1: face[0] = m_numNoeuds[4]; face[1] = m_numNoeuds[5]; face[2] = m_numNoeuds[6]; face[3] = m_numNoeuds[7]; break;
    case 2: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[3]; face[2] = m_numNoeuds[7]; face[3] = m_numNoeuds[4]; break;
    case 3: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[6]; face[3] = m_numNoeuds[5]; break;
    case 4: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[4]; break;
    case 5: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[6]; face[3] = m_numNoeuds[7]; break;
    }
    sommeNoeuds = face[0] + face[1] + face[2] + face[3];
    std::sort(face, face + 4);
    //Recherche existance faces
    indexFaceExiste = FaceNS::rechercheFace(face, sommeNoeuds, facesTemp, sommeNoeudsTemp, 4, iMax);
    if (indexFaceExiste != -1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************