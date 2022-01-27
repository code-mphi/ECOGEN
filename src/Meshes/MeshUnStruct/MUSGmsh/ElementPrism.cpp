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

#include "ElementPrism.h"

const int ElementPrism::TYPEGMSH = 6;
const int ElementPrism::NUMBERNODES = 6;
const int ElementPrism::NUMBERFACES = 5; /* Here there are 3 quadrangles and 2 triangles */
const int ElementPrism::TYPEVTK = 13;

//***********************************************************************

ElementPrism::ElementPrism() :
ElementNS(TYPEGMSH, NUMBERNODES, NUMBERFACES, TYPEVTK)
{}

//***********************************************************************

ElementPrism::~ElementPrism(){}

//***********************************************************************

void ElementPrism::computeVolume(const Coord* nodes)
{
  // Prism volume is computed using 3 tetrahedron
  Coord v1, v2, v3;
  v1.setFromSubtractedVectors(nodes[0], nodes[1]); v2.setFromSubtractedVectors(nodes[0], nodes[2]); v3.setFromSubtractedVectors(nodes[0], nodes[3]);
  double volumeT1 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Volume tetrahedron
  v1.setFromSubtractedVectors(nodes[3], nodes[4]); v2.setFromSubtractedVectors(nodes[3], nodes[5]); v3.setFromSubtractedVectors(nodes[3], nodes[2]);
  double volumeT2 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Volume tetrahedron
  v1.setFromSubtractedVectors(nodes[3], nodes[4]); v2.setFromSubtractedVectors(nodes[3], nodes[2]); v3.setFromSubtractedVectors(nodes[3], nodes[1]);
  double volumeT3 = std::fabs(Coord::determinant(v1, v2, v3)) / 6.; // Volume tetrahedron
  m_volume = volumeT1 + volumeT2 + volumeT3; // Volume prism
}

//***********************************************************************

void ElementPrism::computeLCFL(const Coord* nodes)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((nodes[0] + nodes[1] + nodes[4] + nodes[3]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[0] + nodes[2] + nodes[5] + nodes[3]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[1] + nodes[2] + nodes[5] + nodes[4]) / 4.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[0] + nodes[1] + nodes[2]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
  vec = ((nodes[3] + nodes[4] + nodes[5]) / 3.) - m_position;
  m_lCFL = std::min(m_lCFL, vec.norm());
}

//***********************************************************************

void ElementPrism::construitFaces(const Coord* nodes, FaceNS** faces, int& iMax, int** facesBuff, int* sumNodesBuff)
{
  // 3 faces are quadrangles and 2 faces are triangles
  int indexFaceExists(-1);
  int noeudAutre;
  int currentFaceNodes[4]; // Buffer array of nodes used to create current face
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; facesBuff[iMax][2] = m_numNoeuds[4]; facesBuff[iMax][3] = m_numNoeuds[3]; noeudAutre = 2; break;
      case 1: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[2]; facesBuff[iMax][2] = m_numNoeuds[5]; facesBuff[iMax][3] = m_numNoeuds[3]; noeudAutre = 1; break;      
      case 2: facesBuff[iMax][0] = m_numNoeuds[1]; facesBuff[iMax][1] = m_numNoeuds[2]; facesBuff[iMax][2] = m_numNoeuds[5]; facesBuff[iMax][3] = m_numNoeuds[4]; noeudAutre = 0; break;      
      case 3: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; facesBuff[iMax][2] = m_numNoeuds[2]; noeudAutre = 3; break;      
      case 4: facesBuff[iMax][0] = m_numNoeuds[3]; facesBuff[iMax][1] = m_numNoeuds[4]; facesBuff[iMax][2] = m_numNoeuds[5]; noeudAutre = 0; break;      
    }
    if (i < 3) // Faces Quadrangles
    {
      for (int n = 0; n < 4; n++) { currentFaceNodes[n] = facesBuff[iMax][n]; }
      sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2] + facesBuff[iMax][3];
      std::sort(facesBuff[iMax],facesBuff[iMax]+4);  // Nodes ordering
      // Checking face existence
      indexFaceExists = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,4,iMax);
      // Create face or attach it to element if already existing
      if (indexFaceExists==-1) // faces and facesBuff arrays are filled simultaneously
      {
        faces[iMax] = new FaceQuadrangle(currentFaceNodes[0], currentFaceNodes[1], currentFaceNodes[2], currentFaceNodes[3], 1); // Nodes ordering of quadrangle matters
        faces[iMax]->construitFace(nodes, m_numNoeuds[noeudAutre], this);
        iMax++;
      }
      else
      {
        faces[indexFaceExists]->ajouteElementVoisin(this);
      }
    }
    else // Faces triangles
    {
      sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2];
      std::sort(facesBuff[iMax],facesBuff[iMax]+3); // Nodes ordering 
      // Checking face existence
      indexFaceExists = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,3,iMax);
      // Create face or attach it to element if already existing
      if (indexFaceExists==-1)
      {
        faces[iMax] = new FaceTriangle(facesBuff[iMax][0], facesBuff[iMax][1], facesBuff[iMax][2], 0); // Nodes ordering does not matter for triangle
        faces[iMax]->construitFace(nodes, m_numNoeuds[noeudAutre], this);
        iMax++;
      }
      else
      {
        faces[indexFaceExists]->ajouteElementVoisin(this);
      }
    } // End if
  } // End loop faces
}

//***********************************************************************

void ElementPrism::construitFacesSimplifie(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //3 faces a traiter de type quadrangle et 2 faces triangle
  int indexFaceExists(-1);
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; facesBuff[iMax][2] = m_numNoeuds[4]; facesBuff[iMax][3] = m_numNoeuds[3]; break;
      case 1: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[2]; facesBuff[iMax][2] = m_numNoeuds[5]; facesBuff[iMax][3] = m_numNoeuds[3]; break;      
      case 2: facesBuff[iMax][0] = m_numNoeuds[1]; facesBuff[iMax][1] = m_numNoeuds[2]; facesBuff[iMax][2] = m_numNoeuds[5]; facesBuff[iMax][3] = m_numNoeuds[4]; break;      
      case 3: facesBuff[iMax][0] = m_numNoeuds[0]; facesBuff[iMax][1] = m_numNoeuds[1]; facesBuff[iMax][2] = m_numNoeuds[2]; break;      
      case 4: facesBuff[iMax][0] = m_numNoeuds[3]; facesBuff[iMax][1] = m_numNoeuds[4]; facesBuff[iMax][2] = m_numNoeuds[5]; break;      
    }
    if(i<3)
    {
      sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2] + facesBuff[iMax][3];
      std::sort(facesBuff[iMax], facesBuff[iMax] + 4);  //Tri des nodes
      // Checking face existence
      indexFaceExists = FaceNS::searchFace(facesBuff[iMax], sumNodesBuff[iMax], facesBuff, sumNodesBuff, 4, iMax);

    }
    else
    {
      sumNodesBuff[iMax] = facesBuff[iMax][0] + facesBuff[iMax][1] + facesBuff[iMax][2];
      std::sort(facesBuff[iMax],facesBuff[iMax]+3);  //Tri des nodes
      // Checking face existence
      indexFaceExists = FaceNS::searchFace(facesBuff[iMax],sumNodesBuff[iMax],facesBuff,sumNodesBuff,3,iMax);
    } //Fin if
    if (indexFaceExists == -1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementPrism::attributFaceCommunicante(FaceNS** faces, const int& indexMaxFaces, const int& numberNoeudsInternes)
{
  int indexFaceExists(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[4] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[4], m_numNoeuds[3]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExists))
    {
      faces[indexFaceExists]->ajouteElementVoisinLimite(this);
      faces[indexFaceExists]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes && m_numNoeuds[3] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[2], m_numNoeuds[5], m_numNoeuds[3]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExists))
    {
      faces[indexFaceExists]->ajouteElementVoisinLimite(this);
      faces[indexFaceExists]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes && m_numNoeuds[4] < numberNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[5], m_numNoeuds[4]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExists))
    {
      faces[indexFaceExists]->ajouteElementVoisinLimite(this);
      faces[indexFaceExists]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[0] < numberNoeudsInternes && m_numNoeuds[1] < numberNoeudsInternes && m_numNoeuds[2] < numberNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExists))
    {
      faces[indexFaceExists]->ajouteElementVoisinLimite(this);
      faces[indexFaceExists]->setEstComm(true);
    }
  }
  //Verification face 5 :
  if (m_numNoeuds[3] < numberNoeudsInternes && m_numNoeuds[4] < numberNoeudsInternes && m_numNoeuds[5] < numberNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[3], m_numNoeuds[4], m_numNoeuds[5]);
    if (face.faceExists(faces, indexMaxFaces, indexFaceExists))
    {
      faces[indexFaceExists]->ajouteElementVoisinLimite(this);
      faces[indexFaceExists]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementPrism::compteFaceCommunicante(std::vector<int*>& facesBuff, std::vector<int>& sumNodesBuff)
{
  //3 faces a traiter de type quadrangle et 2 faces triangle
  int indexFaceExists(-1), numberFacesCommunicante(0);
  int face[4], sumNodes;
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[4]; face[3] = m_numNoeuds[3]; break;
      case 1: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[3]; break;
      case 2: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[4]; break;     
      case 3: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; break;
      case 4: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[4]; face[2] = m_numNoeuds[5]; break;
    }
    int iMax = sumNodesBuff.size();
    if(i<3)
    {
      sumNodes = face[0] + face[1] + face[2] + face[3];
      std::sort(face, face+4);
      //Recherche existance faces
      indexFaceExists = FaceNS::searchFace(face, sumNodes, facesBuff, sumNodesBuff, 4, iMax);
    }
    else
    {
      sumNodes = face[0]+face[1]+face[2];
      std::sort(face, face+3);
      //Recherche existance faces
      indexFaceExists = FaceNS::searchFace(face,sumNodes,facesBuff,sumNodesBuff,3,iMax);
    }
    if (indexFaceExists != -1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementPrism::compteFaceCommunicante(int& iMax, int** facesBuff, int* sumNodesBuff)
{
  //3 faces a traiter de type quadrangle et 2 faces triangle
  int indexFaceExists(-1), numberFacesCommunicante(0);
  int face[4], sumNodes;
  for (int i = 0; i < NUMBERFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[4]; face[3] = m_numNoeuds[3]; break;
      case 1: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[3]; break;
      case 2: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[4]; break;     
      case 3: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; break;
      case 4: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[4]; face[2] = m_numNoeuds[5]; break;
    }
    if(i<3)
    {
      sumNodes = face[0] + face[1] + face[2] + face[3];
      std::sort(face, face + 4);
      //Recherche existance faces
      indexFaceExists = FaceNS::searchFace(face, sumNodes, facesBuff, sumNodesBuff, 4, iMax);
    }
    else
    {
      sumNodes = face[0]+face[1]+face[2];
      std::sort(face, face+3);
      //Recherche existance faces
      indexFaceExists = FaceNS::searchFace(face,sumNodes,facesBuff,sumNodesBuff,3,iMax);
    }
    if (indexFaceExists != -1)
    {
      numberFacesCommunicante++;
    }
  }
  return numberFacesCommunicante;
}

//***********************************************************************
