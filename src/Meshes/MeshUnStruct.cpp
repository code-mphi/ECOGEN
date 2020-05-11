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

//! \file      MeshUnStruct.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include <cmath>
#include <algorithm>
#include "MeshUnStruct.h"
#include "../Errors.h"

using namespace tinyxml2;

//***********************************************************************

MeshUnStruct::MeshUnStruct(const std::string &meshFile, const std::string &meshExtension) :
  Mesh(),
  m_meshFile(meshFile),
  m_nameMesh(meshFile),
  m_numberNodes(0),
  m_numberInnerNodes(0),
  m_numberInnerElements(0),
  m_numberGhostElements(0),
  m_numberCommunicatingElements(0),
  m_numberGhostCells(0),
  m_numberFacesParallel(0),
  m_numberElements0D(0),
  m_numberElements1D(0),
  m_numberElements2D(0),
  m_numberElements3D(0),
  m_numberSegments(0),
  m_numberTriangles(0),
  m_numberQuadrangles(0),
  m_numberTetrahedrons(0),
  m_numberPyramids(0),
  m_numberPoints(0),
  m_numberHexahedrons(0),
  m_totalSurface(0.),
  m_totalVolume(0.)
{
  m_meshFile = "./libMeshes/" + m_nameMesh; // meshFile with extension and full path
  m_nameMesh.resize(m_nameMesh.size() - meshExtension.size() - 1); // Remove mesh file extension
  m_type = UNS;

  if (rankCpu == 0) {
	  std::cout << "------------------------------------------------------" << std::endl;
	  std::cout << " B) SEARCHING MESH FILE " + m_meshFile << std::endl;
  }
  std::ifstream meshFileStream(m_meshFile.c_str(), std::ios::in);
  if (!meshFileStream){ throw ErrorECOGEN("mesh file not found : " + m_meshFile, __FILE__, __LINE__); }
  meshFileStream.close();
}

//***********************************************************************

MeshUnStruct::~MeshUnStruct(){
  delete[] m_faces;
  delete[] m_elements;
  delete[] m_nodes;
  for (unsigned int l = 0; l < m_bound.size(); l++) { delete m_bound[l]; }
}

//***********************************************************************

std::string MeshUnStruct::readMeshFileExtension(const std::string &meshFile)
{ 
	if (rankCpu == 0) {
		std::cout << "------------------------------------------------------" << std::endl;
		std::cout << " A) READING MESH FILE EXTENSION : " + meshFile << std::endl;
	}

  return meshFile.substr(meshFile.length() - 3, 3); //three last characters
}

//***********************************************************************

void MeshUnStruct::attributLimites(std::vector<BoundCond*> &boundCond)
{
//JC//Q// If the user specifies an high number for a defined boundary but not specifies numbers for other boundaries
		// you will fill other boundaries with non-reflecting.
		// But if you have a number of bc lower than the bc to imposed on the geometry 
		// you will not attribute non-reflecting to these bc

  // Look for highest boundary physic number
  int maxNumLim(0); 
  for (unsigned int i = 0; i < boundCond.size(); i++) {
    if (boundCond[i]->getNumPhys() > maxNumLim) maxNumLim = boundCond[i]->getNumPhys();
  }
  // Assignment of boundaries in m_bound array in sequence
  int limite(1);
  int limiteTrouvee(0);
  while (limite <= maxNumLim) {
    for (unsigned int i = 0; i < boundCond.size(); i++) {
      if (boundCond[i]->getNumPhys() == limite) {
        m_bound.push_back(boundCond[i]);
        limiteTrouvee = 1; break;
      }
    }
    if (!limiteTrouvee) m_bound.push_back(new BoundCondNonReflecting(limite));
    limite++; limiteTrouvee = 0;
  }
}

//***********************************************************************

int MeshUnStruct::initializeGeometrie(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces,
  const int &restartSimulation, bool pretraitementParallele, std::string ordreCalcul)
{
  try {
    if (Ncpu == 1) { this->initGeometryMonoCPU(cells, cellInterfaces, ordreCalcul); }
    else {
      // Preprocessing of mesh file with CPU 0
      if (pretraitementParallele) {
        if (rankCpu == 0) { this->preProcessMeshFileForParallel(); }
        MPI_Barrier(MPI_COMM_WORLD);
      }
      this->initGeometryParallel(cells, cellsGhost, cellInterfaces, ordreCalcul);
    }
    return m_geometrie;
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::writeMeshInfoData() const
{
	std::cout << "  --------------------------" << std::endl;
	std::cout << "    MESH INFORMATIONS :" << std::endl;
	std::cout << "  --------------------------" << std::endl;
	std::cout << "    mesh nodes number : " << m_numberNodes << std::endl;
	std::cout << "    ~~~~~~~~~~~~~~~~~~~~" << std::endl;
	if (m_numberElements0D != 0)
	{
		std::cout << std::endl << "    0D elements number : " << m_numberElements0D << std::endl;
		std::cout << "    ~~~~~~~~~~~~~~~~~~~~" << std::endl;
		if (m_numberPoints != 0) { std::cout << "      - number vertex : " << m_numberPoints << std::endl; }
	}
	if (m_numberElements1D != 0)
	{
		std::cout << std::endl << "    1D elements number : " << m_numberElements1D << std::endl;
		std::cout << "    ~~~~~~~~~~~~~~~~~~~~" << std::endl;
		if (m_numberSegments != 0) { std::cout << "      - number segments : " << m_numberSegments << std::endl; }
	}
	if (m_numberElements2D != 0)
	{
		std::cout << std::endl << "    2D elements number : " << m_numberElements2D << std::endl;
		std::cout << "    ~~~~~~~~~~~~~~~~~~~~" << std::endl;
		if (m_numberTriangles != 0) { std::cout << "      - number triangles   : " << m_numberTriangles << std::endl; }
		if (m_numberQuadrangles != 0) { std::cout << "      - number quadrangles : " << m_numberQuadrangles << std::endl; }
		std::cout << "      Total surface : " << m_totalSurface << " m2" << std::endl;
	}
	if (m_numberElements3D != 0)
	{
		std::cout << std::endl << "    3D elements number : " << m_numberElements3D << std::endl;
		std::cout << "    ~~~~~~~~~~~~~~~~~~~~" << std::endl;
		if (m_numberTetrahedrons != 0) { std::cout << "      - number tetraedres   : " << m_numberTetrahedrons << std::endl; }
		if (m_numberPyramids != 0) { std::cout << "      - number pyramides    : " << m_numberPyramids << std::endl; }
		if (m_numberHexahedrons != 0) { std::cout << "      - number hexaedres    : " << m_numberHexahedrons << std::endl; }
		std::cout << "      Total volume : " << m_totalVolume << " m3" << std::endl;
	}
	std::cout << std::endl << "... READING MESH FILE COMPLETE " << std::endl;
}

//***********************************************************************

// void MeshUnStruct::readGmshV4(std::vector<ElementNS*>** voisinsNoeuds, std::ifstream &meshFile)
// {
//   try {
//     std::string currentLine;

//     //1) Filling m_nodes array
//     //-------------------------
//     { 
//       std::cout << "  1/Mesh nodes reading ...";
//       while (currentLine != "$Nodes") {
//         getline(meshFile, currentLine);
//         if (meshFile.eof()) { throw ErrorECOGEN("Nodes block not found in mesh file", __FILE__, __LINE__); }
//       }
//       int numEntityBlocks;
//       { std::stringstream lineTotreat;
//       getline(meshFile, currentLine); lineTotreat << currentLine;
//       lineTotreat >> numEntityBlocks >> m_numberNodes; }
//       m_nodes = new Coord[m_numberNodes];
//       //Allocate node array
//       *voisinsNoeuds = new std::vector<ElementNS*>[m_numberNodes];
//       //Reading nodes
//       int tagEntity, dimEntity, parametric, numNodesEntity, tag; double x, y, z;
//       for (int e = 0; e < numEntityBlocks; e++) { //Reading entity
//         {std::stringstream lineTotreat;
//         getline(meshFile, currentLine); lineTotreat << currentLine;
//         lineTotreat >> tagEntity >> dimEntity >> parametric >> numNodesEntity; }
//         for (int i = 0; i < numNodesEntity; i++) {
//           std::stringstream lineTotreat;
//           getline(meshFile, currentLine); lineTotreat << currentLine;
//           lineTotreat >> tag >> x >> y >> z;
//           m_nodes[i].setXYZ(x, y, z);
//         }
//       }
//       getline(meshFile, currentLine);
//       if (currentLine != "$EndNodes") { throw ErrorECOGEN("Nodes block not completely read in mesh file", __FILE__, __LINE__); }
//       std::cout << "OK" << std::endl;
//     }

//     //2) 1D/2D/3D elements are stored in m_elements array / counting
//     //--------------------------------------------------------------
//     {
//       std::cout << "  2/0D/1D/2D/3D elements reading ...";
//       while (currentLine != "$Elements") {
//         getline(meshFile, currentLine);
//         if (meshFile.eof()) { throw ErrorECOGEN("Elements block not found in mesh file", __FILE__, __LINE__); }
//       }
//       int numEntityBlocks;
//       { std::stringstream lineTotreat;
//       getline(meshFile, currentLine); lineTotreat << currentLine;
//       lineTotreat >> numEntityBlocks >> m_numberElements; }
//       //Allocate elements array
//       m_elements = new ElementNS*[m_numberElements];
//       //Readin elements and geometrical properties attributions
//       m_numberElements1D = 0, m_numberElements2D = 0, m_numberElements3D = 0;
//       int posDebutElement = static_cast<int>(meshFile.tellg()); //reperage debut des elements pour retour rapide
//       int noeudG;
//       int tagEntity, dimEntity, typeEle, numElementsEntity;
//       int numElement(0);
//       for (int e = 0; e < numEntityBlocks; e++) { //Reading entity
//         {std::stringstream lineTotreat;
//         getline(meshFile, currentLine); lineTotreat << currentLine;
//         lineTotreat >> tagEntity >> dimEntity >> typeEle >> numElementsEntity; }
//         for (int i = 0; i < numElementsEntity; i++) {
//           this->lectureElementGmshV4(m_nodes, meshFile, &m_elements[numElement], typeEle,numElement,tagEntity);
//           //Les tags entity ne sont pas bon ni l'ordonancement des faces et cellules => a revoir
//           if (m_elements[numElement]->getTypeGmsh() == 15) { m_numberElements0D++; }
//           else if (m_elements[numElement]->getTypeGmsh() == 1) { m_numberElements1D++; }
//           else if (m_elements[numElement]->getTypeGmsh() <= 3) { m_numberElements2D++; m_totalSurface += m_elements[i]->getVolume(); }
//           else if (m_elements[numElement]->getTypeGmsh() <= 7) { m_numberElements3D++; m_totalVolume += m_elements[i]->getVolume(); }
//           else { throw ErrorECOGEN("Type element du .msh non gere dans ECOGEN", __FILE__, __LINE__); }
//           //Attribution element i voisin pour les noeuds concernes (Ordre 2 muiltislopes)
//           for (int n = 0; n < m_elements[numElement]->getNumberNoeuds(); n++) {
//             noeudG = m_elements[numElement]->getNumNoeud(n);
//             (*voisinsNoeuds)[noeudG].push_back(m_elements[numElement]);
//           }
//           numElement++;
//         }
//       }
//       m_numberInnerElements = m_numberElements;
//       getline(meshFile, currentLine);
//       if (currentLine != "$EndElements") { throw ErrorECOGEN("Elements block not completely read in mesh file", __FILE__, __LINE__); }
//       std::cout << "OK" << std::endl;
//     }

//   }
//   catch (ErrorECOGEN &) { throw; }
// }


//***********************************************************************

// void MeshUnStruct::lectureElementGmshV4(const Coord *TableauNoeuds, std::ifstream &fichierMesh, ElementNS **element, const int &typeElement, int &indiceElement, const int & physicalEntity)
// {
//   try {
//     //1)Number of vertex affectation
//     //------------------------------
//     switch (typeElement)
//     {
//     case 1: //segment (deux points)
//       *element = new ElementSegment;
//       m_numberSegments++;
//       break;
//     case 2: //triangle (trois points)
//       *element = new ElementTriangle;
//       m_numberTriangles++;
//       break;
//     case 3: //Quadrangle (quatre points)
//       *element = new ElementQuadrangle;
//       m_numberQuadrangles++;
//       break;
//     case 4: //Tetrahedron (quatre points)
//       *element = new ElementTetrahedron;
//       m_numberTetrahedrons++;
//       break;
//     case 7: //Pyramid quadrangulaire (cinq points)
//       *element = new ElementPyramid;
//       m_numberPyramids++;
//       break;
//     case 15: //Point (un vertex)
//       *element = new ElementPoint;
//       m_numberPoints++;
//       break;
//     case 5: //Hexahedron (huit points)
//       *element = new ElementHexahedron;
//       m_numberHexahedrons++;
//       break;
//     case 6: //Prism (six points)
//       *element = new ElementPrism;
//       m_numberHexahedrons++;
//       break;
//     default:
//       throw ErrorECOGEN("Element type unknown in mesh file", __FILE__, __LINE__);
//       break;
//     } //Fin switch typeElement

//     //2) Element building / properties filling
//     //----------------------------------------
//     int noeudCourant, tag;
//     int *numNoeud = new int[(*element)->getNumberNoeuds()];
//     Coord *noeud = new Coord[(*element)->getNumberNoeuds()];
//     std::string currentLine;
//     std::stringstream lineTotreat;
//     getline(fichierMesh, currentLine); lineTotreat << currentLine;
//     lineTotreat >> tag;
//     for (int i = 0; i < (*element)->getNumberNoeuds(); i++)
//     {
//       lineTotreat >> noeudCourant;
//       numNoeud[i] = noeudCourant - 1;         //decalage car tableau commencant a zero
//       noeud[i] = TableauNoeuds[noeudCourant - 1];
//     }
//     (*element)->construitElement(numNoeud, noeud, physicalEntity, 0, indiceElement);

//     delete[] noeud;
//     delete[] numNoeud;
//   }
//   catch (ErrorECOGEN &) { throw; }

// }


//**************************************************************************
//******************************** ECRITURE ********************************
//**************************************************************************

void MeshUnStruct::ecritHeaderPiece(std::ofstream &fileStream, TypeMeshContainer<Cell *> *cellsLvl) const
{
  fileStream << "    <Piece NumberOfPoints=\"" << m_numberNodes << "\" NumberOfCells=\"" << m_numberCellsCalcul - m_numberGhostCells << "\">" << std::endl;
}

//****************************************************************************

void MeshUnStruct::recupereNoeuds(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  for (int noeud = 0; noeud < m_numberNodes; noeud++)
  {
    jeuDonnees.push_back(m_nodes[noeud].getX());
    jeuDonnees.push_back(m_nodes[noeud].getY());
    jeuDonnees.push_back(m_nodes[noeud].getZ());
  }
}

//****************************************************************************

void MeshUnStruct::recupereConnectivite(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      for (int noeud = 0; noeud < m_elements[i]->getNumberNoeuds(); noeud++)
      {
        jeuDonnees.push_back(m_elements[i]->getNumNoeud(noeud));
      }
    }
  }
}

//****************************************************************************

void MeshUnStruct::recupereOffsets(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int offset(0);
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      offset += m_elements[i]->getNumberNoeuds();
      jeuDonnees.push_back(offset);
    }
  }
}

//****************************************************************************

void MeshUnStruct::recupereTypeCell(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      jeuDonnees.push_back(m_elements[i]->getTypeVTK());
    }
  }
}

//****************************************************************************

void MeshUnStruct::recupereDonnees(TypeMeshContainer<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, const int var, int phase) const
{
  jeuDonnees.clear();
  int numCell;
  double transport(0.);
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      if (var > 0) { //On veut recuperer les donnees scalars
        if (phase >= 0) { jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnScalar(var)); }      //Donnees de phases
        else if (phase == -1) { jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnScalar(var)); }   //Donnees de mixture
        else if (phase == -2) {
          transport = cellsLvl[0][numCell]->getTransport(var-1).getValue();
          if (transport < 1.e-20) { transport = 0.; }
          jeuDonnees.push_back(transport);
        }
        else if (phase == -3) { jeuDonnees.push_back(cellsLvl[0][numCell]->getXi()); }
        else if (phase == -4) { jeuDonnees.push_back(cellsLvl[0][numCell]->getDensityGradient()); }
        else { Errors::errorMessage("MeshUnStruct::recupereDonnees: unknown number of phase: ", phase); }
      }
      else { //On veut recuperer les donnees vectorielles
        if (phase >= 0) { //Phases data
          jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getX());
          jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getY());
          jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getZ());
        }
        else if(phase == -1){  //Mixture data
          jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getX());
          jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getY());
          jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getZ());
        }
        else { Errors::errorMessage("MeshUnStruct::recupereDonnees: unknown number of phase: ", phase); }
      } //Fin vecteur
    }
  }
}

//****************************************************************************

void MeshUnStruct::setDataSet(std::vector<double> &jeuDonnees, TypeMeshContainer<Cell *> *cellsLvl, const int var, int phase) const
{
  int iterDataSet(0);
  int numCell;
  Coord vec;
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      if (var > 0) { //Scalars data are first set
        if (phase >= 0) { cellsLvl[0][numCell]->getPhase(phase)->setScalar(var, jeuDonnees[iterDataSet++]); } //phases data
        else if (phase == -1) { cellsLvl[0][numCell]->getMixture()->setScalar(var, jeuDonnees[iterDataSet++]); }  //mixture data
        else if (phase == -2) { cellsLvl[0][numCell]->getTransport(var - 1).setValue(jeuDonnees[iterDataSet++]); } //transport data
        else { Errors::errorMessage("MeshUnStruct::setDataSet: unknown phase number: ", phase); }
      }
      else { //On veut recuperer les donnees vectorielles
        if (phase >= 0) { //Phases data
          vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet+1], jeuDonnees[iterDataSet+2]);
          cellsLvl[0][numCell]->getPhase(phase)->setVector(-var, vec);
          iterDataSet += 3;
        }
        else if (phase == -1) {  //Mixture data
          vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet+1], jeuDonnees[iterDataSet+2]);
          cellsLvl[0][numCell]->getMixture()->setVector(-var, vec);
          iterDataSet += 3;
        }
        else { Errors::errorMessage("MeshUnStruct::setDataSet: unknown phase number: ", phase); }
      } //Fin vecteur
    }
  }
}

//****************************************************************************
void MeshUnStruct::extractAbsVeloxityMRF(TypeMeshContainer<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, Source *sourceMRF) const
{
  jeuDonnees.clear();
  int numCell;
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      Coord absoluteVelocity = sourceMRF->computeAbsVelocity(cellsLvl[0][numCell]->getVelocity(), cellsLvl[0][numCell]->getPosition());
      jeuDonnees.push_back(absoluteVelocity.getX());
      jeuDonnees.push_back(absoluteVelocity.getY());
      jeuDonnees.push_back(absoluteVelocity.getZ());
    }
  }
}

//***********************************************************************

// void MeshUnStruct::rechercheElementsArrieres(ElementNS *element, FaceNS *face, CellInterface *cellInterface, std::vector<ElementNS *> voisins, Cell **cells) const
// {
//   //int index(element->getIndex() - m_numberBoundFaces);
//   //Coord vG = element->vecteur(face);
//   //ElementNS *eT; //Element a tester
//   ////1) Recherche du premier element arriere gauche BG1M
//   ////---------------------------------------------------
//   //double cos(0.);
//   //for (unsigned int v = 0; v < voisins.size(); v++) {
//   //  eT = voisins[v];
//   //  Coord vT = eT->vecteur(element);
//   //  double cosTemp(Coord::cos(vT, vG));
//   //  if (cosTemp >= cos) {
//   //    cos = cosTemp;
//   //    if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1M a NULL
//   //      cellInterface->setB(BG1M, 0);
//   //    }
//   //    else { //Sinon on met a jour avec la maille la plus loin
//   //      Cell *c(cells[eT->getNumCellAssociee()]);
//   //      cellInterface->setB(BG1M, c);
//   //    }
//   //  }
//   //}
//   ////2) Recherche du second element arriere gauche BG2M
//   ////---------------------------------------------------
//   //cos = 0.;
//   //if (cellInterface->getB(BG1M) != 0) {  //Cas non cellInterface
//   //  Element *e1(cellInterface->getB(BG1M)->getElement());
//   //  Coord v1 = e1->vecteur(element);
//   //  Coord sin1(Coord::sin(v1, vG));
//   //  if (std::fabs(sin1.norm()) <= 1e-8) {
//   //    cellInterface->setB(BG2M, cellInterface->getB(BG1M)); // Si le cos est nul, meme cell pour B1 et B2
//   //  }
//   //  else {
//   //    for (unsigned int v = 0; v < voisins.size(); v++) {
//   //      eT = voisins[v];
//   //      if (eT == e1) continue; // voisin suivant si voisin deja utilise 
//   //      Coord vT = eT->vecteur(element);
//   //      Coord sinTemp(Coord::sin(vT, vG));

//   //      if (sinTemp.scalar(sin1) <= 0.) {
//   //        double cosTemp(Coord::cos(vT, vG));
//   //        if (cosTemp >= cos) {
//   //          cos = cosTemp;
//   //          if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1M a NULL
//   //            cellInterface->setB(BG2M, 0);
//   //          }
//   //          else {  //Sinon on met a jour avec la 2 eme maille la plus loin
//   //            Cell *c(cells[eT->getNumCellAssociee()]);
//   //            cellInterface->setB(BG2M, c);
//   //          }
//   //        }
//   //      } //fin sin*sin <0
//   //    }  //fin boucle voisins
//   //  } //fin if 
//   //} //Fin if cellInterface

//   //  //Determination des ponderations arrieres
//   //Coord MB1; if (cellInterface->getB(BG1M) != 0) { MB1.setFromSubtractedVectors(face->getPos(), cellInterface->getB(BG1M)->getPosition()); }
//   //Coord MB2; if (cellInterface->getB(BG2M) != 0) { MB2.setFromSubtractedVectors(face->getPos(), cellInterface->getB(BG2M)->getPosition()); }
//   //Coord MB; MB.setFromSubtractedVectors(face->getPos(), element->getPosition());
//   //double a, b, beta1(1.), beta2(0.);
//   //a = (MB1.vectoriel(MB)).norm() / MB.norm();
//   //b = (MB2.vectoriel(MB)).norm() / MB.norm();
//   //if (std::fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
//   //cellInterface->setBeta(betaG1M, beta1);
//   //cellInterface->setBeta(betaG2M, beta2);
//   ////Calcul de la distance du vertex M (cell interface) au vertex H (barycenter)
//   //if (cellInterface->getB(BG1M) != 0) {
//   //  cos = Coord::cos(MB, MB1);
//   //  double d, e, c;
//   //  d = cos*MB1.norm();
//   //  Coord B1B2; B1B2 = MB1 - MB2;
//   //  e = B1B2.norm()*beta1;
//   //  c = sqrt(e*e - a*a);
//   //  d += c;
//   //  cellInterface->setDistanceH(distanceHGM, d);
//   //}
//   //cout << eG->getIndex() << " : " << endl;
//   //if (cellInterfaces[i]->getB(BG1M) != 0) cout << "BG1M = " << cellInterfaces[i]->getB(BG1M)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHGM) << endl;
//   //if (cellInterfaces[i]->getB(BG2M) != 0) cout << "BG2M = " << cellInterfaces[i]->getB(BG2M)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHGM) << endl;
// }

//***********************************************************************

// void MeshUnStruct::rechercheElementsAvants(ElementNS *element, FaceNS *face, CellInterface *cellInterface, std::vector<ElementNS *> voisins, Cell **cells) const
// {
//   int index(element->getIndex() - m_numberBoundFaces);
//   Coord vG = element->vecteur(face);
//   //ElementNS *eT; //Element a tester
//   ////1) Recherche du premier element arriere gauche BG1P
//   ////---------------------------------------------------
//   //double cos(0.);
//   //for (unsigned int v = 0; v < voisins.size(); v++) {
//   //  eT = voisins[v];
//   //  Coord vT = eT->vecteur(element);
//   //  double cosTemp(Coord::cos(vT, vG));
//   //  if (cosTemp <= cos) {
//   //    cos = cosTemp;
//   //    if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1P a NULL
//   //      cellInterface->setB(BG1P, 0);
//   //    }
//   //    else { //Sinon on met a jour avec la maille la plus loin
//   //      Cell *c(cells[eT->getNumCellAssociee()]);
//   //      cellInterface->setB(BG1P, c);
//   //    }
//   //  }
//   //}
//   ////2) Recherche du second element arriere gauche BG2P
//   ////---------------------------------------------------
//   //cos = 0.;
//   //if (cellInterface->getB(BG1P) != 0) {  //Cas non cellInterface
//   //  Element *e1(cellInterface->getB(BG1P)->getElement());
//   //  Coord v1 = e1->vecteur(element);
//   //  Coord sin1(Coord::sin(v1, vG));
//   //  if (std::fabs(sin1.norm()) <= 1e-8) {
//   //    cellInterface->setB(BG2P, cellInterface->getB(BG1P)); // Si le cos est nul, meme cell pour B1 et B2
//   //  }
//   //  else {
//   //    for (unsigned int v = 0; v < voisins.size(); v++) {
//   //      eT = voisins[v];
//   //      if (eT == e1) continue; // voisin suivant si voisin deja utilise 
//   //      Coord vT = eT->vecteur(element);
//   //      Coord sinTemp(Coord::sin(vT, vG));

//   //      if (sinTemp.scalar(sin1) <= 0.) {
//   //        double cosTemp(Coord::cos(vT, vG));
//   //        if (cosTemp <= cos) {
//   //          cos = cosTemp;
//   //          if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG2P a NULL
//   //            cellInterface->setB(BG2P, 0);
//   //          }
//   //          else {  //Sinon on met a jour avec la 2 eme maille la plus loin
//   //            Cell *c(cells[eT->getNumCellAssociee()]);
//   //            cellInterface->setB(BG2P, c);
//   //          }
//   //        }
//   //      } //fin sin*sin <0
//   //    }  //fin boucle voisins
//   //  } //fin if 
//   //} //Fin if cellInterface

//   //  //Determination des ponderations arrieres
//   //Coord MB1; if (cellInterface->getB(BG1P) != 0) { MB1.setFromSubtractedVectors(face->getPos(), cellInterface->getB(BG1P)->getPosition()); }
//   //Coord MB2; if (cellInterface->getB(BG2P) != 0) { MB2.setFromSubtractedVectors(face->getPos(), cellInterface->getB(BG1P)->getPosition()); }
//   //Coord MB; MB.setFromSubtractedVectors(face->getPos(), element->getPosition());
//   //double a, b, beta1(1.), beta2(0.);
//   //a = (MB1.vectoriel(MB)).norm() / MB.norm();
//   //b = (MB2.vectoriel(MB)).norm() / MB.norm();
//   //if (std::fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
//   //cellInterface->setBeta(betaG1P, beta1);
//   //cellInterface->setBeta(betaG2P, beta2);
//   ////Calcul de la distance du vertex M (cellInterface de maille) au vertex H (barycenter)
//   //if (cellInterface->getB(BG1P) != 0) {
//   //  cos = Coord::cos(MB, -1.*MB1);
//   //  double d, e, c;
//   //  d = cos*MB1.norm();
//   //  Coord B1B2; B1B2 = MB1 - MB2;
//   //  e = B1B2.norm()*beta1;
//   //  c = sqrt(e*e - a*a);
//   //  d += c;
//   //  cellInterface->setDistanceH(distanceHGP, d);
//   //}
//   //cout << eG->getIndex() << " : " << endl;
//   //if (cellInterfaces[i]->getB(BG1P) != 0) cout << "BG1P = " << cellInterfaces[i]->getB(BG1P)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHGP) << endl;
//   //if (cellInterfaces[i]->getB(BG2P) != 0) cout << "BG2P = " << cellInterfaces[i]->getB(BG2P)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHGP) << endl;
// }

//***********************************************************************
