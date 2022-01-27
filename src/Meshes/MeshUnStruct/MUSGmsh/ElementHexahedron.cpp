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

#include "ElementHexahedron.h"

const int ElementHexahedron::TYPEGMSH = 5;
const int ElementHexahedron::NUMBERNODES = 8;
const int ElementHexahedron::NUMBERFACES = 6; /* ici il s'agit de quadrangles*/
const int ElementHexahedron::TYPEVTK = 12;

//***********************************************************************

ElementHexahedron::ElementHexahedron() :
ElementNS(TYPEGMSH, NUMBERNODES, NUMBERFACES, TYPEVTK)
{}

//***********************************************************************

ElementHexahedron::~ElementHexahedron(){}

//***********************************************************************

void ElementHexahedron::computeVolume(const Coord* nodes)
{
  // Volume of hexahedron = volume of 6 tetrahedron
  //JC can be done with only 5 tetrahedron
  Coord v1, v2, v3;
  m_volume = 0.; // Volume of hexahedron
  // v1.setFromSubtractedVectors(nodes[0], nodes[1]); v2.setFromSubtractedVectors(nodes[0], nodes[3]); v3.setFromSubtractedVectors(nodes[0], nodes[4]);
  // m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 1st tetrahedron
  // v1.setFromSubtractedVectors(nodes[3], nodes[1]); v2.setFromSubtractedVectors(nodes[3], nodes[7]); v3.setFromSubtractedVectors(nodes[3], nodes[4]);
  // m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 2nd tetrahedron
  // v1.setFromSubtractedVectors(nodes[3], nodes[1]); v2.setFromSubtractedVectors(nodes[3], nodes[7]); v3.setFromSubtractedVectors(nodes[3], nodes[2]);
  // m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 3rd tetrahedron
  // v1.setFromSubtractedVectors(nodes[4], nodes[1]); v2.setFromSubtractedVectors(nodes[4], nodes[7]); v3.setFromSubtractedVectors(nodes[4], nodes[5]);
  // m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 4th tetrahedron
  // v1.setFromSubtractedVectors(nodes[7], nodes[2]); v2.setFromSubtractedVectors(nodes[7], nodes[6]); v3.setFromSubtractedVectors(nodes[7], nodes[5]);
  // m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 5th tetrahedron
  // v1.setFromSubtractedVectors(nodes[7], nodes[2]); v2.setFromSubtractedVectors(nodes[7], nodes[5]); v3.setFromSubtractedVectors(nodes[7], nodes[1]);
  // m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 6th tetrahedron

  // Joris' version with 6 tetrahedron
  v1.setFromSubtractedVectors(nodes[0], nodes[3]); v2.setFromSubtractedVectors(nodes[0], nodes[2]); v3.setFromSubtractedVectors(nodes[0], nodes[4]);
  m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 1st tetrahedron
  v1.setFromSubtractedVectors(nodes[3], nodes[7]); v2.setFromSubtractedVectors(nodes[3], nodes[2]); v3.setFromSubtractedVectors(nodes[3], nodes[4]);
  m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 2nd tetrahedron
  v1.setFromSubtractedVectors(nodes[0], nodes[1]); v2.setFromSubtractedVectors(nodes[0], nodes[4]); v3.setFromSubtractedVectors(nodes[0], nodes[2]);
  m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 3rd tetrahedron
  v1.setFromSubtractedVectors(nodes[7], nodes[2]); v2.setFromSubtractedVectors(nodes[7], nodes[6]); v3.setFromSubtractedVectors(nodes[7], nodes[4]);
  m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 4th tetrahedron
  v1.setFromSubtractedVectors(nodes[4], nodes[6]); v2.setFromSubtractedVectors(nodes[4], nodes[5]); v3.setFromSubtractedVectors(nodes[4], nodes[2]);
  m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 5th tetrahedron
  v1.setFromSubtractedVectors(nodes[4], nodes[5]); v2.setFromSubtractedVectors(nodes[4], nodes[2]); v3.setFromSubtractedVectors(nodes[4], nodes[1]);
  m_volume += std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Add volume of 6th tetrahedron
}

//***********************************************************************

void ElementHexahedron::computeLCFL(const Coord* nodes)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((nodes[0] + nodes[1] + nodes[2] + nodes[3]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[4] + nodes[5] + nodes[6] + nodes[7]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[0] + nodes[4] + nodes[7] + nodes[3]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[1] + nodes[5] + nodes[6] + nodes[2]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[0] + nodes[1] + nodes[5] + nodes[4]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[3] + nodes[2] + nodes[6] + nodes[7]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************

void ElementHexahedron::construitFaces(const Coord* nodes, FaceNS** faces, int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //6 faces a traiter de type quadrangle
  int indexFaceExiste(-1);
  int noeudAutre;
  int currentFaceNodes[4]; // buffer array of node used to create current face
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: currentFaceNodes[0] = m_numNoeuds[0]; currentFaceNodes[1] = m_numNoeuds[1]; currentFaceNodes[2] = m_numNoeuds[2]; currentFaceNodes[3] = m_numNoeuds[3]; noeudAutre = 4; break;
      case 1: currentFaceNodes[0] = m_numNoeuds[4]; currentFaceNodes[1] = m_numNoeuds[5]; currentFaceNodes[2] = m_numNoeuds[6]; currentFaceNodes[3] = m_numNoeuds[7]; noeudAutre = 0; break;
      case 2: currentFaceNodes[0] = m_numNoeuds[0]; currentFaceNodes[1] = m_numNoeuds[3]; currentFaceNodes[2] = m_numNoeuds[7]; currentFaceNodes[3] = m_numNoeuds[4]; noeudAutre = 1; break;
      case 3: currentFaceNodes[0] = m_numNoeuds[1]; currentFaceNodes[1] = m_numNoeuds[2]; currentFaceNodes[2] = m_numNoeuds[6]; currentFaceNodes[3] = m_numNoeuds[5]; noeudAutre = 0; break;
      case 4: currentFaceNodes[0] = m_numNoeuds[0]; currentFaceNodes[1] = m_numNoeuds[1]; currentFaceNodes[2] = m_numNoeuds[5]; currentFaceNodes[3] = m_numNoeuds[4]; noeudAutre = 2; break;
      case 5: currentFaceNodes[0] = m_numNoeuds[3]; currentFaceNodes[1] = m_numNoeuds[2]; currentFaceNodes[2] = m_numNoeuds[6]; currentFaceNodes[3] = m_numNoeuds[7]; noeudAutre = 0; break;
    }
    for(int n=0; n<4; n++){ facesBuff[iMax][n] = currentFaceNodes[n]; } // Filling facesBuff array before sorting
    sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2] + facesBuff[iMax][3];
    std::sort(facesBuff[iMax],facesBuff[iMax]+4);  //Tri des nodes
    // Checking face existence
    indexFaceExiste = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,4,iMax);
    //Creation face ou rattachement
    if (indexFaceExiste==-1) // on fill simultanement le tableau faces et le tableau facesBuff
    {
      faces[iMax] = new FaceQuadrangle(currentFaceNodes[0], currentFaceNodes[1], currentFaceNodes[2], currentFaceNodes[3], 1); // Nodes ordering of quadrangle matters
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

void ElementHexahedron::construitFacesSimplifie(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //6 faces a traiter de type quadrangle
  int indexFaceExiste(-1);
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
    case 0: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; facesBuff[iMax][2] = m_numNoeuds[2]; facesBuff[iMax][3] = m_numNoeuds[3]; break;
    case 1: facesBuff[iMax][0] = m_numNoeuds[4]; facesBuff[iMax][1] = m_numNoeuds[5]; facesBuff[iMax][2] = m_numNoeuds[6]; facesBuff[iMax][3] = m_numNoeuds[7]; break;
    case 2: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[3]; facesBuff[iMax][2] = m_numNoeuds[7]; facesBuff[iMax][3] = m_numNoeuds[4]; break;
    case 3: facesBuff[iMax][0] = m_numNoeuds[1]; facesBuff[iMax][1] = m_numNoeuds[2]; facesBuff[iMax][2] = m_numNoeuds[6]; facesBuff[iMax][3] = m_numNoeuds[5]; break;
    case 4: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; facesBuff[iMax][2] = m_numNoeuds[5]; facesBuff[iMax][3] = m_numNoeuds[4]; break;
    case 5: facesBuff[iMax][0] = m_numNoeuds[3]; facesBuff[iMax][1] = m_numNoeuds[2]; facesBuff[iMax][2] = m_numNoeuds[6]; facesBuff[iMax][3] = m_numNoeuds[7]; break;
    }
    sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2] + facesBuff[iMax][3];
    std::sort(facesBuff[iMax], facesBuff[iMax] + 4);  //Tri des nodes
    // Checking face existence
    indexFaceExiste = FaceNS::searchFace(facesBuff[iMax], sumNodesBuff[iMax], facesBuff, sumNodesBuff, 4, iMax);
    //Creation face ou rattachement
    if (indexFaceExiste == -1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementHexahedron::attributFaceCommunicante(FaceNS** faces, const int& indexMaxFaces, const int& numberNoeudsInternes)
{
  int indexFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[3]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[4] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes && m_numNoeuds[6] < numberNoeudsInternes && m_numNoeuds[7] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[4], m_numNoeuds[5], m_numNoeuds[6], m_numNoeuds[7]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[7] < numberNoeudsInternes && m_numNoeuds[4] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[3], m_numNoeuds[7], m_numNoeuds[4]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[6] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[6], m_numNoeuds[5]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 5 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes && m_numNoeuds[4] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[5], m_numNoeuds[4]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 6 :
  if (m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[6] < numberNoeudsInternes && m_numNoeuds[7] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[3], m_numNoeuds[2], m_numNoeuds[6], m_numNoeuds[7]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementHexahedron::compteFaceCommunicante(std::vector<int*>& facesBuff, std::vector<int>& sumNodesBuff)
{
  //6 faces a traiter de type quadrangle
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[4], sumNodes;
  for (int i = 0; i < NUMBERFACES; i++)
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
    int iMax = sumNodesBuff.size();
    sumNodes = face[0] + face[1] + face[2] + face[3];
    std::sort(face, face + 4);
    //Recherche existance faces
    indexFaceExiste = FaceNS::searchFace(face, sumNodes, facesBuff, sumNodesBuff, 4, iMax);
    if (indexFaceExiste != -1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementHexahedron::compteFaceCommunicante(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //6 faces a traiter de type quadrangle
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[4], sumNodes;
  for (int i = 0; i < NUMBERFACES; i++)
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
    sumNodes = face[0] + face[1] + face[2] + face[3];
    std::sort(face, face + 4);
    //Recherche existance faces
    indexFaceExiste = FaceNS::searchFace(face, sumNodes, facesBuff, sumNodesBuff, 4, iMax);
    if (indexFaceExiste != -1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************