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

#include "MeshUnStruct.h"
#include "../Errors.h"

using namespace tinyxml2;

//***********************************************************************

MeshUnStruct::MeshUnStruct(const std::string& meshFile, const std::string& meshExtension) :
  Mesh(),
  m_meshFile(meshFile),
  m_nameMesh(meshFile),
  m_numberNodes(0),
  m_numberInnerNodes(0),
  m_numberInnerElements(0),
  m_numberGhostElements(0),
  m_numberCommunicatingElements(0),
  m_numberFacesParallel(0),
  m_numberGhostCells(0),
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
  m_meshFile = m_nameMesh; // meshFile with extension and full path
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

std::string MeshUnStruct::readMeshFileExtension(const std::string& meshFile)
{ 
	if (rankCpu == 0) {
		std::cout << "------------------------------------------------------" << std::endl;
		std::cout << " A) READING MESH FILE EXTENSION : " + meshFile << std::endl;
	}

  return meshFile.substr(meshFile.length() - 3, 3); //three last characters
}

//***********************************************************************

void MeshUnStruct::attributLimites(std::vector<BoundCond*>& boundCond)
{
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

int MeshUnStruct::initializeGeometrie(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces,
  const int& /*restartSimulation*/, bool pretraitementParallele, std::string ordreCalcul)
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
	std::cout << "    MESH INFORMATIONS:" << std::endl;
	std::cout << "  --------------------------" << std::endl;
	std::cout << "    Number of mesh nodes: " << m_numberNodes << std::endl;
	std::cout << "    ~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	if (m_numberElements0D != 0)
	{
		std::cout << std::endl;
    std::cout << "    Number of 0D elements: " << m_numberElements0D << std::endl;
		std::cout << "    ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		if (m_numberPoints != 0) { std::cout << "      - number of vertices: " << m_numberPoints << std::endl; }
	}
	if (m_numberElements1D != 0)
	{
		std::cout << std::endl; 
    std::cout << "    Number of 1D elements: " << m_numberElements1D << std::endl;
		std::cout << "    ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		if (m_numberSegments != 0) { std::cout << "      - number of segments: " << m_numberSegments << std::endl; }
	}
	if (m_numberElements2D != 0)
	{
		std::cout << std::endl;
    std::cout << "    Number of 2D elements: " << m_numberElements2D << std::endl;
		std::cout << "    ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		if (m_numberTriangles != 0) { std::cout << "      - number of triangles: " << m_numberTriangles << std::endl; }
		if (m_numberQuadrangles != 0) { std::cout << "      - number of quadrangles: " << m_numberQuadrangles << std::endl; }
		std::cout << "      Total surface: " << m_totalSurface << " m2" << std::endl;
	}
	if (m_numberElements3D != 0)
	{
		std::cout << std::endl;
    std::cout << "    Number of 3D elements: " << m_numberElements3D << std::endl;
		std::cout << "    ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		if (m_numberTetrahedrons != 0) { std::cout << "      - number of tetrahedrons: " << m_numberTetrahedrons << std::endl; }
		if (m_numberPyramids != 0) { std::cout << "      - number of pyramids: " << m_numberPyramids << std::endl; }
		if (m_numberHexahedrons != 0) { std::cout << "      - number of hexahedrons: " << m_numberHexahedrons << std::endl; }
		std::cout << "      Total volume: " << m_totalVolume << " m3" << std::endl;
	}
  // Read reference length used for CFL criterion
  double min(1.e5), max(0.), lref(0.);
  for (int i = 0; i < m_numberCellsCalcul; i++) {
    lref = m_elements[i + m_numberBoundFaces]->getLCFL();
    min = std::min(min, lref);
    max = std::max(max, lref);
  }
  std::cout << std::endl;
  std::cout << "    Reference length: " << std::endl;
  std::cout << "    ~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "      - min = " << min << std::endl;
  std::cout << "      - max = " << max << std::endl;

  std::cout << std::endl << "... READING MESH FILE COMPLETE " << std::endl;
}

//**************************************************************************
//******************************** ECRITURE ********************************
//**************************************************************************

void MeshUnStruct::ecritHeaderPiece(std::ofstream& fileStream, TypeMeshContainer<Cell*>* /*cellsLvl*/) const
{
  fileStream << "    <Piece NumberOfPoints=\"" << m_numberNodes << "\" NumberOfCells=\"" << m_numberCellsCalcul - m_numberGhostCells << "\">" << std::endl;
}

//****************************************************************************

void MeshUnStruct::recupereNoeuds(std::vector<double>& jeuDonnees, std::vector<Cell*>* /*cellsLvl*/) const
{
  for (int noeud = 0; noeud < m_numberNodes; noeud++)
  {
    jeuDonnees.push_back(m_nodes[noeud].getX());
    jeuDonnees.push_back(m_nodes[noeud].getY());
    jeuDonnees.push_back(m_nodes[noeud].getZ());
  }
}

//****************************************************************************

void MeshUnStruct::recupereConnectivite(std::vector<double>& jeuDonnees, std::vector<Cell*>* /*cellsLvl*/) const
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

void MeshUnStruct::recupereOffsets(std::vector<double>& jeuDonnees, std::vector<Cell*>* /*cellsLvl*/) const
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

void MeshUnStruct::recupereTypeCell(std::vector<double>& jeuDonnees, std::vector<Cell*>* /*cellsLvl*/) const
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

void MeshUnStruct::recupereDonnees(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& jeuDonnees, const int var, int phase) const
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

void MeshUnStruct::setDataSet(std::vector<double>& jeuDonnees, TypeMeshContainer<Cell*>* cellsLvl, const int var, int phase) const
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

void MeshUnStruct::extractAbsVelocityMRF(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& jeuDonnees, Source *sourceMRF) const
{
  jeuDonnees.clear();
  int numCell;
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      // Absolute velocity is built on the specific region rotating or when whole geometry is rotating.
      if (sourceMRF->getPhysicalEntity() == cellsLvl[0][numCell]->getElement()->getAppartenancePhysique() || sourceMRF->getPhysicalEntity() == 0) {
        Coord absoluteVelocity = sourceMRF->computeAbsVelocity(cellsLvl[0][numCell]->getVelocity(), cellsLvl[0][numCell]->getPosition());
        jeuDonnees.push_back(absoluteVelocity.getX());
        jeuDonnees.push_back(absoluteVelocity.getY());
        jeuDonnees.push_back(absoluteVelocity.getZ());
      }
      // If the region is not rotating absolute velocity = relative velocity
      else {
        jeuDonnees.push_back(cellsLvl[0][numCell]->getVelocity().getX());
        jeuDonnees.push_back(cellsLvl[0][numCell]->getVelocity().getY());
        jeuDonnees.push_back(cellsLvl[0][numCell]->getVelocity().getZ());
      }
    }
  }
}

//***********************************************************************

// void MeshUnStruct::rechercheElementsArrieres(ElementNS *element, FaceNS *face, CellInterface* cellInterface, std::vector<ElementNS *> voisins, Cell** cells) const
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
//   //      Cell* c(cells[eT->getNumCellAssociee()]);
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
//   //            Cell* c(cells[eT->getNumCellAssociee()]);
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

// void MeshUnStruct::rechercheElementsAvants(ElementNS *element, FaceNS *face, CellInterface* cellInterface, std::vector<ElementNS *> voisins, Cell** cells) const
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
//   //      Cell* c(cells[eT->getNumCellAssociee()]);
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
//   //            Cell* c(cells[eT->getNumCellAssociee()]);
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
