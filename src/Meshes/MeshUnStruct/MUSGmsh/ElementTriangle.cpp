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

#include "ElementTriangle.h"

const int ElementTriangle::TYPEGMSH = 2;
const int ElementTriangle::NUMBERNODES = 3;
const int ElementTriangle::NUMBERFACES = 3; /* Here it is the number of segments */
const int ElementTriangle::TYPEVTK = 5;

//***********************************************************************

ElementTriangle::ElementTriangle() :
ElementNS(TYPEGMSH, NUMBERNODES, NUMBERFACES, TYPEVTK)
{}

//***********************************************************************

ElementTriangle::~ElementTriangle(){}

//***********************************************************************

void ElementTriangle::computeVolume(const Coord* nodes)
{
  Coord v1(nodes[1] - nodes[0]), v2(nodes[2] - nodes[1]), v3(nodes[0] - nodes[2]);
  double a(v1.norm()), b(v2.norm()), c(v3.norm());
  double dp = 0.5*(a + b + c);
  m_volume = sqrt(dp*(dp - a)*(dp - b)*(dp - c)); // Triangle area using Heron's formula

  // Kahan's formula helps to prevent numerical instability
  // especially for needle like triangles
  // It is required to have a >= b >= c 
  // and to use additionnal parenthesis to preserve computation order
  // if (a < b) Tools::swap(a, b);
  // if (a < c) Tools::swap(a, c);
  // if (b < c) Tools::swap(b, c);
  // m_volume = sqrt( (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)) ) * 0.25; 
}

//***********************************************************************

void ElementTriangle::computeLCFL(const Coord* nodes)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((nodes[0] + nodes[1]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[0] + nodes[2]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[1] + nodes[2]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************

void ElementTriangle::construitFaces(const Coord* nodes, FaceNS** faces, int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //3 faces a traiter de type segment
  int indexFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; noeudAutre = 2; break;
      case 1: facesBuff[iMax][0] = m_numNoeuds[1]; facesBuff[iMax][1] = m_numNoeuds[2]; noeudAutre = 0; break;      
      case 2: facesBuff[iMax][0] = m_numNoeuds[2]; facesBuff[iMax][1] = m_numNoeuds[0]; noeudAutre = 1; break;      
    }
    sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1];
    std::sort(facesBuff[iMax],facesBuff[iMax]+2);  //Tri des nodes
    // Checking face existence
    indexFaceExiste = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,2,iMax);
    //Creation face ou rattachement
    if (indexFaceExiste==-1)
    {
      faces[iMax] = new FaceSegment(facesBuff[iMax][0], facesBuff[iMax][1], 0); //pas besoin du tri ici
      faces[iMax]->construitFace(nodes, m_numNoeuds[noeudAutre], this);
      iMax++;
    }
    else
    {
      faces[indexFaceExiste]->ajouteElementVoisin(this);
    }
  }
}

//***********************************************************************

void ElementTriangle::construitFacesSimplifie(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //3 faces a traiter de type segment
  int indexFaceExiste(-1);
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; break;
      case 1: facesBuff[iMax][0] = m_numNoeuds[1]; facesBuff[iMax][1] = m_numNoeuds[2]; break;      
      case 2: facesBuff[iMax][0] = m_numNoeuds[2]; facesBuff[iMax][1] = m_numNoeuds[0]; break;      
    }
    sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1];
    std::sort(facesBuff[iMax],facesBuff[iMax]+2);  //Tri des nodes
    // Checking face existence
    indexFaceExiste = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,2,iMax);
    //Creation face ou rattachement
    if (indexFaceExiste==-1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementTriangle::attributFaceLimite(FaceNS** faces, const int& indexMaxFaces)
{
  int indexFaceExiste(0);
  FaceTriangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2]);
  if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
  {
    faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
  }
  else
  {
    Errors::errorMessage("Probleme attribution des faces limites element Triangle");
  }
}

//***********************************************************************

void ElementTriangle::attributFaceCommunicante(FaceNS** faces, const int& indexMaxFaces, const int& numberNoeudsInternes)
{
  int indexFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[0], m_numNoeuds[1]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[1], m_numNoeuds[2]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[0] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[2], m_numNoeuds[0]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementTriangle::compteFaceCommunicante(std::vector<int*>& facesBuff, std::vector<int>& sumNodesBuff)
{
  //4 faces a traiter de type segment
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[2], sumNodes;
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[0]; break;     
    }
    int iMax = sumNodesBuff.size();
    sumNodes = face[0]+face[1];
    std::sort(face, face+2);
    //Recherche existance faces
    indexFaceExiste = FaceNS::searchFace(face,sumNodes,facesBuff,sumNodesBuff,2,iMax);
    if (indexFaceExiste!=-1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementTriangle::compteFaceCommunicante(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //4 faces a traiter de type segment
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[2], sumNodes;
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[0]; break;       
    }
    sumNodes = face[0]+face[1];
    std::sort(face, face+2);
    //Recherche existance faces
    indexFaceExiste = FaceNS::searchFace(face,sumNodes,facesBuff,sumNodesBuff,2,iMax);
    if (indexFaceExiste!=-1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************