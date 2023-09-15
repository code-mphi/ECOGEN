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
  m_nodes(nullptr),
  m_numberInnerElements(0),
  m_numberGhostElements(0),
  m_numberCommunicatingElements(0),
  m_elements(nullptr),
  m_faces(nullptr),
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

  //JC//WARNING this message will be displayed even during the reading of the rough mesh
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
  if (m_faces != nullptr) delete[] m_faces;
  if (m_elements != nullptr) delete[] m_elements;
  if (m_nodes != nullptr) delete[] m_nodes;
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

void MeshUnStruct::assignLimits(std::vector<BoundCond*>& boundCond)
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
    return m_problemDimension;
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
//******************************** WRITING *********************************
//**************************************************************************

void MeshUnStruct::writeHeaderPiece(std::ofstream& fileStream, TypeMeshContainer<Cell*>* /*cellsLvl*/) const
{
  fileStream << "    <Piece NumberOfPoints=\"" << m_numberNodes << "\" NumberOfCells=\"" << m_numberCellsCalcul - m_numberGhostCells << "\">" << std::endl;
}

//****************************************************************************

void MeshUnStruct::getNodes(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const
{
  for (int node = 0; node < m_numberNodes; node++)
  {
    dataset.push_back(m_nodes[node].getX());
    dataset.push_back(m_nodes[node].getY());
    dataset.push_back(m_nodes[node].getZ());
  }
}

//****************************************************************************

void MeshUnStruct::getConnectivity(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const
{
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      for (int node = 0; node < m_elements[i]->getNumberNodes(); node++)
      {
        dataset.push_back(m_elements[i]->getNumNode(node));
      }
    }
  }
}

//****************************************************************************

void MeshUnStruct::getOffsets(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const
{
  int offset(0);
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      offset += m_elements[i]->getNumberNodes();
      dataset.push_back(offset);
    }
  }
}

//****************************************************************************

void MeshUnStruct::getTypeCell(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const
{
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      dataset.push_back(m_elements[i]->getTypeVTK());
    }
  }
}

//****************************************************************************

void MeshUnStruct::getData(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, const int var, int phase) const
{
  dataset.clear();
  int numCell;
  double transport(0.);
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      if (var > 0) { //We want to get the scalar data
        if (phase >= 0) { dataset.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnScalar(var)); }      //data de phases
        else if (phase == -1) { dataset.push_back(cellsLvl[0][numCell]->getMixture()->returnScalar(var)); }   //data de mixture
        else if (phase == -2) {
          transport = cellsLvl[0][numCell]->getTransport(var-1).getValue();
          if (transport < 1.e-20) { transport = 0.; }
          dataset.push_back(transport);
        }
        else if (phase == -3) { dataset.push_back(cellsLvl[0][numCell]->getXi()); }
        else if (phase == -4) { dataset.push_back(cellsLvl[0][numCell]->getDensityGradient()); }
        else if (phase == -7) { // Saturation pressure
          dataset.push_back(cellsLvl[0][numCell]->getPsat());
        }
        else { Errors::errorMessage("MeshUnStruct::getData: unknown number of phase: ", phase); }
      }
      else { //We want to get the vector data
        if (phase >= 0) { //Phases data
          dataset.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getX());
          dataset.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getY());
          dataset.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getZ());
        }
        else if(phase == -1){  //Mixture data
          dataset.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getX());
          dataset.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getY());
          dataset.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getZ());
        }
        else { Errors::errorMessage("MeshUnStruct::getData: unknown number of phase: ", phase); }
      } //End vector
    }
  }
}

//****************************************************************************

void MeshUnStruct::setDataSet(std::vector<double>& dataset, TypeMeshContainer<Cell*>* cellsLvl, const int var, int phase) const
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
        if (phase >= 0) { cellsLvl[0][numCell]->getPhase(phase)->setScalar(var, dataset[iterDataSet++]); } //phases data
        else if (phase == -1) { cellsLvl[0][numCell]->getMixture()->setScalar(var, dataset[iterDataSet++]); }  //mixture data
        else if (phase == -2) { cellsLvl[0][numCell]->getTransport(var - 1).setValue(dataset[iterDataSet++]); } //transport data
        else { Errors::errorMessage("MeshUnStruct::setDataSet: unknown phase number: ", phase); }
      }
      else { //We want to get the vector data
        if (phase >= 0) { //Phases data
          vec.setXYZ(dataset[iterDataSet], dataset[iterDataSet+1], dataset[iterDataSet+2]);
          cellsLvl[0][numCell]->getPhase(phase)->setVector(-var, vec);
          iterDataSet += 3;
        }
        else if (phase == -1) {  //Mixture data
          vec.setXYZ(dataset[iterDataSet], dataset[iterDataSet+1], dataset[iterDataSet+2]);
          cellsLvl[0][numCell]->getMixture()->setVector(-var, vec);
          iterDataSet += 3;
        }
        else { Errors::errorMessage("MeshUnStruct::setDataSet: unknown phase number: ", phase); }
      } //End vector
    }
  }
}

//****************************************************************************

void MeshUnStruct::extractAbsVelocityMRF(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, Source *sourceMRF) const
{
  dataset.clear();
  int numCell;
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      // Absolute velocity is built on the specific region rotating or when whole geometry is rotating.
      if (sourceMRF->getPhysicalEntity() == cellsLvl[0][numCell]->getElement()->getAppartenancePhysique() || sourceMRF->getPhysicalEntity() == 0) {
        Coord absoluteVelocity = sourceMRF->computeAbsVelocity(cellsLvl[0][numCell]->getVelocity(), cellsLvl[0][numCell]->getPosition());
        dataset.push_back(absoluteVelocity.getX());
        dataset.push_back(absoluteVelocity.getY());
        dataset.push_back(absoluteVelocity.getZ());
      }
      // If the region is not rotating absolute velocity = relative velocity
      else {
        dataset.push_back(cellsLvl[0][numCell]->getVelocity().getX());
        dataset.push_back(cellsLvl[0][numCell]->getVelocity().getY());
        dataset.push_back(cellsLvl[0][numCell]->getVelocity().getZ());
      }
    }
  }
}

//***********************************************************************

void MeshUnStruct::extractReferenceLength(std::vector<Cell*>* cellsLvl, std::vector<double>& dataset) const
{
  dataset.clear();
  int numCell;
  for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      dataset.push_back(cellsLvl[0][numCell]->getElement()->getLCFL());
    }
  }
}

//***********************************************************************