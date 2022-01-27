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

#include "ElementQuadrangle.h"

const int ElementQuadrangle::TYPEGMSH = 3;
const int ElementQuadrangle::NUMBERNODES = 4;
const int ElementQuadrangle::NUMBERFACES = 4; /* Here it is the number of segments */
const int ElementQuadrangle::TYPEVTK = 9;

//***********************************************************************

ElementQuadrangle::ElementQuadrangle() :
ElementNS(TYPEGMSH, NUMBERNODES, NUMBERFACES, TYPEVTK)
{}

//***********************************************************************

ElementQuadrangle::~ElementQuadrangle(){}

//***********************************************************************

void ElementQuadrangle::computeVolume(const Coord* nodes)
{
  //A diagonal
  Coord v0(nodes[2] - nodes[0]);
  double diagonale(v0.norm());
  //The 4 sides
  Coord v1(nodes[1] - nodes[0]);
  Coord v2(nodes[2] - nodes[1]);
  Coord v3(nodes[3] - nodes[2]);
  Coord v4(nodes[0] - nodes[3]);
  double a(v1.norm()); double b(v2.norm()); double c(v3.norm()); double d(v4.norm());
  //Area 1st triangle
  double dp1 = 0.5*(a + b + diagonale);
  double surf1 = sqrt(dp1*(dp1 - a)*(dp1 - b)*(dp1 - diagonale));
  //Area 2nd triangle
  double dp2 = 0.5*(c + d + diagonale);
  double surf2 = sqrt(dp2*(dp2 - c)*(dp2 - d)*(dp2 - diagonale));
  //Area quadrangle
  m_volume = surf1 + surf2;
}

//***********************************************************************

void ElementQuadrangle::computeLCFL(const Coord* nodes)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((nodes[0] + nodes[1]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[1] + nodes[2]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[2] + nodes[3]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[3] + nodes[0]) / 2.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************

void ElementQuadrangle::construitFaces(const Coord* nodes, FaceNS** faces, int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //4 faces a traiter de type segment
  int indexFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; noeudAutre = 2; break;
      case 1: facesBuff[iMax][0] = m_numNoeuds[1]; facesBuff[iMax][1] = m_numNoeuds[2]; noeudAutre = 3; break;      
      case 2: facesBuff[iMax][0] = m_numNoeuds[2]; facesBuff[iMax][1] = m_numNoeuds[3]; noeudAutre = 0; break;      
      case 3: facesBuff[iMax][0] = m_numNoeuds[3]; facesBuff[iMax][1] = m_numNoeuds[0]; noeudAutre = 1; break;      
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

void ElementQuadrangle::construitFacesSimplifie(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //4 faces a traiter de type segment
  int indexFaceExiste(-1);
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; break;
      case 1: facesBuff[iMax][0] = m_numNoeuds[1]; facesBuff[iMax][1] = m_numNoeuds[2]; break;      
      case 2: facesBuff[iMax][0] = m_numNoeuds[2]; facesBuff[iMax][1] = m_numNoeuds[3]; break;      
      case 3: facesBuff[iMax][0] = m_numNoeuds[3]; facesBuff[iMax][1] = m_numNoeuds[0]; break;      
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

void ElementQuadrangle::attributFaceLimite(FaceNS** faces, const int& indexMaxFaces)
{
  int indexFaceExiste(0);
  FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[3]);
  if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
  {
    faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
  }
  else
  {
    Errors::errorMessage("Probleme attribution des faces limites element Quadrangle");
  }
}

//***********************************************************************

void ElementQuadrangle::attributFaceCommunicante(FaceNS** faces, const int& indexMaxFaces, const int& numberNoeudsInternes)
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
  if (m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[2], m_numNoeuds[3]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[0] < numberNoeudsInternes)
  {
    FaceSegment face(m_numNoeuds[3], m_numNoeuds[0]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementQuadrangle::compteFaceCommunicante(std::vector<int*>& facesBuff, std::vector<int>& sumNodesBuff)
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
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; break;     
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
int ElementQuadrangle::compteFaceCommunicante(int& iMax, int** facesBuff, int* sumNodesBuff)
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
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; break;     
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