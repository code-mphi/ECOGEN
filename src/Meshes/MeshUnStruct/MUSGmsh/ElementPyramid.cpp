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

#include "ElementPyramid.h"

const int ElementPyramid::TYPEGMSH = 7;
const int ElementPyramid::NUMBERNODES = 5;
const int ElementPyramid::NUMBERFACES = 5; /* Here there are 1 quadrangle and 4 triangles*/
const int ElementPyramid::TYPEVTK = 14;

//***********************************************************************

ElementPyramid::ElementPyramid() :
ElementNS(TYPEGMSH, NUMBERNODES, NUMBERFACES, TYPEVTK)
{}

//***********************************************************************

ElementPyramid::~ElementPyramid(){}

//***********************************************************************

void ElementPyramid::computeVolume(const Coord* nodes)
{
  // The volume is computed using the two tetrahedrons included in the pyramid
  Coord v1, v2, v3;
  v1.setFromSubtractedVectors(nodes[4], nodes[1]); v2.setFromSubtractedVectors(nodes[4], nodes[0]); v3.setFromSubtractedVectors(nodes[4], nodes[2]);
  double volumeT1 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Tetrahedron's volume
  v1.setFromSubtractedVectors(nodes[4], nodes[2]); v2.setFromSubtractedVectors(nodes[4], nodes[3]); v3.setFromSubtractedVectors(nodes[4], nodes[0]);
  double volumeT2 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Tetrahedron's volume
  m_volume = volumeT1 + volumeT2; // Pyramid's volume
}

//***********************************************************************

void ElementPyramid::computeLCFL(const Coord* nodes)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((nodes[0] + nodes[1] + nodes[2] + nodes[3]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[0] + nodes[1] + nodes[4]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[1] + nodes[2] + nodes[4]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[2] + nodes[3] + nodes[4]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[3] + nodes[0] + nodes[4]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************

void ElementPyramid::construitFaces(const Coord* nodes, FaceNS** faces, int& iMax, int** facesBuff, int* sumNodesBuff)
{
  // 1 face is quadrangle and 4 faces are triangles
  int indexFaceExiste(-1);
  int otherNode;

  // Building quadrangle as basis of pyramid
  // ---------------------------------------
  int currentFaceNodes[4]; // Buffer array of node used to create current face (nodes ordering matters for quadrangle's face)
  currentFaceNodes[0] = m_numNodes[0];
  currentFaceNodes[1] = m_numNodes[1];
  currentFaceNodes[2] = m_numNodes[2];
  currentFaceNodes[3] = m_numNodes[3];
  otherNode = 4;

  for(int n=0; n<4; n++){ facesBuff[iMax][n] = currentFaceNodes[n]; } // Filling facesBuff array before sorting
  sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2] + facesBuff[iMax][3];
  std::sort(facesBuff[iMax],facesBuff[iMax]+4); // Nodes ordering to check if face exists
  // Checking face existence
  indexFaceExiste = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,4,iMax);
  // Create face or attach it to element if already existing
  if (indexFaceExiste==-1) // faces and facesBuff arrays are filled simultaneously
  {
    faces[iMax] = new FaceQuadrangle(currentFaceNodes[0], currentFaceNodes[1], currentFaceNodes[2], currentFaceNodes[3], 1); // Nodes ordering of quadrangle matters
    faces[iMax]->construitFace(nodes, m_numNodes[otherNode], this);
    iMax++;
  }
  else
  {
    faces[indexFaceExiste]->addElementNeighbor(this);
  }

  // Building the 4 triangular faces
  // -------------------------------
  for (int i = 1; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 1: facesBuff[iMax][0] = m_numNodes[0]; facesBuff[iMax][1] = m_numNodes[1]; facesBuff[iMax][2] = m_numNodes[4]; otherNode = 2; break;
      case 2: facesBuff[iMax][0] = m_numNodes[1]; facesBuff[iMax][1] = m_numNodes[2]; facesBuff[iMax][2] = m_numNodes[4]; otherNode = 3; break;
      case 3: facesBuff[iMax][0] = m_numNodes[2]; facesBuff[iMax][1] = m_numNodes[3]; facesBuff[iMax][2] = m_numNodes[4]; otherNode = 0; break;
      case 4: facesBuff[iMax][0] = m_numNodes[3]; facesBuff[iMax][1] = m_numNodes[0]; facesBuff[iMax][2] = m_numNodes[4]; otherNode = 1; break;
    }
    sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2];
    std::sort(facesBuff[iMax],facesBuff[iMax]+3); // Nodes ordering
    // Checking face existence
    indexFaceExiste = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,3,iMax);
    // Create face or attach it to element if already existing
    if (indexFaceExiste==-1)
    {
      faces[iMax] = new FaceTriangle(facesBuff[iMax][0], facesBuff[iMax][1], facesBuff[iMax][2], 0); // No need to sort nodes to build triangle's face
      faces[iMax]->construitFace(nodes, m_numNodes[otherNode], this);
      iMax++;
    }
    else
    {
      faces[indexFaceExiste]->addElementNeighbor(this);
    }
  }

}

//***********************************************************************

void ElementPyramid::construitFacesSimplifie(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  // 1 face is quadrangle and 4 faces are triangles
  int indexFaceExiste(-1);

  // Building quadrangle as basis of pyramid
  // ---------------------------------------
  facesBuff[iMax][0] = m_numNodes[0];
  facesBuff[iMax][1] = m_numNodes[1];
  facesBuff[iMax][2] = m_numNodes[2];
  facesBuff[iMax][3] = m_numNodes[3];
  sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2] + facesBuff[iMax][3];
  std::sort(facesBuff[iMax],facesBuff[iMax]+4); // Nodes ordering
  // Checking face existence
  indexFaceExiste = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,4,iMax);
  // Create face or attach it to element if already existing
  if (indexFaceExiste==-1) // faces anf facesBuff arrays are filled simultaneously
  {
    iMax++;
  }

  // Building the 4 triangular faces
  // -------------------------------
  for (int i = 1; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 1: facesBuff[iMax][0] = m_numNodes[0]; facesBuff[iMax][1] = m_numNodes[1]; facesBuff[iMax][2] = m_numNodes[4]; break;
      case 2: facesBuff[iMax][0] = m_numNodes[1]; facesBuff[iMax][1] = m_numNodes[2]; facesBuff[iMax][2] = m_numNodes[4]; break;
      case 3: facesBuff[iMax][0] = m_numNodes[2]; facesBuff[iMax][1] = m_numNodes[3]; facesBuff[iMax][2] = m_numNodes[4]; break;
      case 4: facesBuff[iMax][0] = m_numNodes[3]; facesBuff[iMax][1] = m_numNodes[0]; facesBuff[iMax][2] = m_numNodes[4]; break;
    }
    sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2];
    std::sort(facesBuff[iMax],facesBuff[iMax]+3); // Nodes ordering
    // Checking face existence
    indexFaceExiste = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,3,iMax);
    // Create face or attach it to element if already existing
    if (indexFaceExiste==-1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementPyramid::attributFaceCommunicante(FaceNS** faces, const int& indexMaxFaces, const int& numberNodesInternal)
{
  int indexFaceExiste(0);
  // Verification face Quadrangle:
  if (m_numNodes[0] < numberNodesInternal && m_numNodes[1] < numberNodesInternal && m_numNodes[2] < numberNodesInternal && m_numNodes[3] < numberNodesInternal)
  {
    FaceQuadrangle face(m_numNodes[0], m_numNodes[1], m_numNodes[2], m_numNodes[3]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->addElementNeighborLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  // Verification face Triangle 1:
  if (m_numNodes[0] < numberNodesInternal && m_numNodes[1] < numberNodesInternal && m_numNodes[4] < numberNodesInternal)
  {
    FaceTriangle face(m_numNodes[0], m_numNodes[1], m_numNodes[4]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->addElementNeighborLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  // Verification face Triangle 2:
  if (m_numNodes[1] < numberNodesInternal && m_numNodes[2] < numberNodesInternal && m_numNodes[4] < numberNodesInternal)
  {
    FaceTriangle face(m_numNodes[1], m_numNodes[2], m_numNodes[4]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->addElementNeighborLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  // Verification face Triangle 3:
  if (m_numNodes[2] < numberNodesInternal && m_numNodes[3] < numberNodesInternal && m_numNodes[4] < numberNodesInternal)
  {
    FaceTriangle face(m_numNodes[2], m_numNodes[3], m_numNodes[4]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->addElementNeighborLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }
  // Verification face Triangle 4:
  if (m_numNodes[3] < numberNodesInternal && m_numNodes[0] < numberNodesInternal && m_numNodes[4] < numberNodesInternal)
  {
    FaceTriangle face(m_numNodes[3], m_numNodes[0], m_numNodes[4]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExiste))
    {
      faces[indexFaceExiste]->addElementNeighborLimite(this);
      faces[indexFaceExiste]->setEstComm(true);
    }
  }    
}

//***********************************************************************


int ElementPyramid::compteFaceCommunicante(std::vector<int*>& facesBuff, std::vector<int>& sumNodesBuff)
{
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[4], sumNodes, iMax;

  // 1 quadrangular face
  face[0] = m_numNodes[0]; 
  face[1] = m_numNodes[1]; 
  face[2] = m_numNodes[2]; 
  face[3] = m_numNodes[3];
  iMax = sumNodesBuff.size();
  sumNodes = face[0] + face[1] + face[2] + face[3];
  std::sort(face, face + 4);
  // Checking face existence
  indexFaceExiste = FaceNS::searchFace(face, sumNodes, facesBuff, sumNodesBuff, 4, iMax);
  if (indexFaceExiste != -1)
  {
    numberFacesCommunicante++;
  }

  // 4 triangular faces
  for (int i = 1; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 1: face[0] = m_numNodes[0]; face[1] = m_numNodes[1]; face[2] = m_numNodes[4]; break;
      case 2: face[0] = m_numNodes[1]; face[1] = m_numNodes[2]; face[2] = m_numNodes[4]; break;
      case 3: face[0] = m_numNodes[2]; face[1] = m_numNodes[3]; face[2] = m_numNodes[4]; break;     
      case 4: face[0] = m_numNodes[3]; face[1] = m_numNodes[0]; face[2] = m_numNodes[4]; break;     
    }
    iMax = sumNodesBuff.size();
    sumNodes = face[0]+face[1]+face[2];
    std::sort(face, face+3);
    //Recherche existance faces
    indexFaceExiste = FaceNS::searchFace(face,sumNodes,facesBuff,sumNodesBuff,3,iMax);
    if (indexFaceExiste!=-1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementPyramid::compteFaceCommunicante(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  int indexFaceExiste(-1), numberFacesCommunicante(0);
  int face[4], sumNodes;

  //1 quadrangle face
  face[0] = m_numNodes[0]; 
  face[1] = m_numNodes[1]; 
  face[2] = m_numNodes[2]; 
  face[3] = m_numNodes[3];
  sumNodes = face[0] + face[1] + face[2] + face[3];
  std::sort(face, face + 4);
  // Checking face existence
  indexFaceExiste = FaceNS::searchFace(face, sumNodes, facesBuff, sumNodesBuff, 4, iMax);
  if (indexFaceExiste != -1)
  {
    numberFacesCommunicante++;
  }

  //4 faces a traiter de type triangle
  for (int i = 1; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 1: face[0] = m_numNodes[0]; face[1] = m_numNodes[1]; face[2] = m_numNodes[4]; break;
      case 2: face[0] = m_numNodes[1]; face[1] = m_numNodes[2]; face[2] = m_numNodes[4]; break;
      case 3: face[0] = m_numNodes[2]; face[1] = m_numNodes[3]; face[2] = m_numNodes[4]; break;     
      case 4: face[0] = m_numNodes[3]; face[1] = m_numNodes[0]; face[2] = m_numNodes[4]; break;     
    }
    sumNodes = face[0]+face[1]+face[2];
    std::sort(face, face+3);
    // Checking face existence
    indexFaceExiste = FaceNS::searchFace(face,sumNodes,facesBuff,sumNodesBuff,3,iMax);
    if (indexFaceExiste!=-1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************