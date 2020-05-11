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

//! \file      MUSGmshV4.cpp
//! \author    J. Caze
//! \version   1.0
//! \date      November 26 2019

#include "MUSGmshV4.h"
#include <iterator>

//***********************************************************************

MUSGmshV4::MUSGmshV4(const std::string &meshFile, const std::string &meshExtension) : MUSGmsh(meshFile,meshExtension)
{}

//***********************************************************************

MUSGmshV4::~MUSGmshV4()
{
}
//***********************************************************************

void MUSGmshV4::initGeometryMonoCPU(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string computeOrder)
{
	try {
		// 1) Reading nodes and elements
		// -----------------------------
		this->readMeshMonoCPU();  // Fill m_nodes and m_elements

		// Display mesh information
		this->writeMeshInfoData();

		// CAUTION: Ordering of m_elements is important. Faces first, then cells.

		std::cout << "------------------------------------------------------" << std::endl;
		std::cout << " D) BUILDING GEOMETRY ..." << std::endl;

		// 2) Assignment of cells to their geometric element
		// -------------------------------------------------
		// Cells and boundary counting
		if (m_numberElements3D == 0 && m_numberElements2D == 0) // 1D case
		{
			m_numberCellsCalcul = m_numberElements1D;
			m_numberBoundFaces = m_numberElements0D;
			m_geometrie = 1;
		}
		else if (m_numberElements3D == 0) // 2D case
		{
			m_numberCellsCalcul = m_numberElements2D;
			m_numberBoundFaces = m_numberElements1D;
			m_geometrie = 2;
		}
		else // 3D case
		{
			m_numberCellsCalcul = m_numberElements3D;
			m_numberBoundFaces = m_numberElements2D;
			m_geometrie = 3;
		}
		m_numberCellsTotal = m_numberCellsCalcul;

		// Assignment of non-reflecting boundary for missing limits
		unsigned int nbLimits(0);
		for (int i = 0; i < m_numberBoundFaces; i++) {
			unsigned int appPhys(m_elements[i]->getAppartenancePhysique());
			if (appPhys > nbLimits) { nbLimits = appPhys; }
		}
		for (unsigned int i = m_bound.size(); i < nbLimits; i++) {
			m_bound.push_back(new BoundCondNonReflecting);
		}

		// Assignment of cells to elements and counting inner faces
		m_numberInnerFaces = 0;
		for (int i = 0; i < m_numberCellsCalcul; i++)
		{
			if (computeOrder == "FIRSTORDER") { cells.push_back(new Cell); }
			else { cells.push_back(new CellO2); }
			cells[i]->setElement(m_elements[i + m_numberBoundFaces], i);
			m_numberInnerFaces += m_elements[i + m_numberBoundFaces]->getNumberFaces();
		}

		m_numberInnerFaces -= m_numberBoundFaces; // We remove the boundaries
		m_numberInnerFaces /= 2; // The internal faces are all counted twice => we restore the truth
		m_numberFacesTotal = m_numberInnerFaces + m_numberBoundFaces;

		double volTot(0.);
		for (int i = 0; i < m_numberCellsCalcul; i++)
		{
			volTot += cells[i]->getElement()->getVolume();
		}
		std::cout << "Total volume : " << volTot << std::endl;

		// 3) Building connectivity table
		// -------------------------------
		// Sizing faces array
		m_faces = new FaceNS * [m_numberFacesTotal];
		int** facesBuff; int* sumNodesBuff; // We build a temporary array of faces to speed up the search process
		facesBuff = new int* [m_numberFacesTotal + 1];
		sumNodesBuff = new int[m_numberFacesTotal + 1];
		// Determination of the maximum number of nodes for the faces 
		int sizeFace; // Will be initialize at maximale size
		if (m_numberElements3D != 0)
		{
			if (m_numberQuadrangles != 0) { sizeFace = 4; }
			else if (m_numberTriangles != 0) { sizeFace = 3; }
			else { Errors::errorMessage("Issue in initGeometryMonoCPU for initialization of facesBuff array"); }
		}
		else if (m_numberElements2D != 0) { sizeFace = 2; }
		else { sizeFace = 1; }
		for (int i = 0; i < m_numberFacesTotal + 1; i++)
		{
			facesBuff[i] = new int[sizeFace];
		}

		// Inner faces
		int indexMaxFaces(0);
		clock_t tTemp(clock());
		float t1(0.);
		std::cout << "  1/Building faces ..." << std::endl;
		int frequenceImpressure(std::max((m_numberElements - m_numberBoundFaces) / 10, 1));
		for (int i = m_numberBoundFaces; i < m_numberElements; i++)
		{
			if ((i - m_numberBoundFaces) % frequenceImpressure == 0) { std::cout << "    " << (100 * (i - m_numberBoundFaces) / (m_numberElements - m_numberBoundFaces)) << "% ... " << std::endl; }
			m_elements[i]->construitFaces(m_nodes, m_faces, indexMaxFaces, facesBuff, sumNodesBuff);
		}
		for (int i = 0; i < m_numberFacesTotal + 1; i++) { delete facesBuff[i]; }
		delete[] facesBuff; delete[] sumNodesBuff;
		tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
		std::cout << "    OK in " << t1 << " seconds" << std::endl;

		// Boundaries
		std::cout << "  2/Boundary elements attribution to boundary faces ..." << std::endl;
		tTemp = clock();
		frequenceImpressure = std::max(m_numberBoundFaces / 10, 1);
		for (int i = 0; i < m_numberBoundFaces; i++)
		{
			if (i % frequenceImpressure == 0)
			{
				std::cout << "    " << (100 * i / m_numberBoundFaces) << "% ... " << std::endl;
			}
			// Assigning the boundary
			m_elements[i]->attributFaceLimite(m_nodes, m_faces, indexMaxFaces);
		}
		tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
		std::cout << "    OK in " << t1 << " seconds" << std::endl;

		// Link Geometry/cellInterfaces of compute
		std::cout << "  3/Linking Geometries -> Physics ..." << std::endl;
		tTemp = clock();
		int iMailleG, iMailleD;
		for (int i = 0; i < m_numberFacesTotal; i++)
		{
			if (m_faces[i]->getEstLimite()) // Physical boundary faces
			{
				int appPhys(m_faces[i]->getElementDroite()->getAppartenancePhysique() - 1); // (appartenance - 1) for array starting at zero
				if (appPhys >= static_cast<int>(m_bound.size()) || appPhys < 0) { Errors::errorMessage("Number of boundary conditions not suited"); }
				m_bound[appPhys]->creeLimite(cellInterfaces);
				cellInterfaces[i]->setFace(m_faces[i]);
				iMailleG = m_faces[i]->getElementGauche()->getIndex() - m_numberBoundFaces;
				iMailleD = iMailleG;
			}
			// Inner faces of domain
			else
			{
				if (computeOrder == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
				else { cellInterfaces.push_back(new CellInterfaceO2); }
				cellInterfaces[i]->setFace(m_faces[i]);
				iMailleG = m_faces[i]->getElementGauche()->getIndex() - m_numberBoundFaces;
				iMailleD = m_faces[i]->getElementDroite()->getIndex() - m_numberBoundFaces;
			}
			cellInterfaces[i]->initialize(cells[iMailleG], cells[iMailleD]);
			cells[iMailleG]->addCellInterface(cellInterfaces[i]);
			cells[iMailleD]->addCellInterface(cellInterfaces[i]);

		} // End face
		tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
		std::cout << "    OK in " << t1 << " seconds" << std::endl;
		std::cout << "... BUILDING GEOMETRY COMPLETE " << std::endl;
		std::cout << "------------------------------------------------------" << std::endl;
	}
	catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

//void MUSGmshV4::preProcessMeshFileForParallel()
//{
//	ElementNS** elementsGlobal;
//	Coord* nodesGlobal;
//	int frequenceImpressure(0);
//	int numberNodesGlobal(0), numberElementsGlobal(0);
//	int numberElements0D(0), numberElements1D(0), numberElements2D(0), numberElements3D(0);
//
//	try {
//		// Opening mesh file
//		// -----------------
//		std::cout << "------------------------------------------------------" << std::endl;
//		std::cout << " C0) MESH FILE PRETRAITEMENT " + m_meshFile + " IN PROGRESS ..." << std::endl;
//		clock_t totalTime(clock());
//		std::ifstream meshFile(m_meshFile.c_str(), std::ios::in);
//		if (!meshFile) { throw ErrorECOGEN("mesh file not found :" + m_meshFile, __FILE__, __LINE__); }
//		std::string currentLine;
//
//		// Skiping Gmsh version + go to nodes
//		// ----------------------------------
//		getline(meshFile, currentLine);
//		getline(meshFile, currentLine);
//		getline(meshFile, currentLine);
//		getline(meshFile, currentLine);
//
//		// 2) Reading entities
//		// -------------------
//		std::cout << "  1/Mesh entities reading ..." << std::endl;
//
//		int numPts, numCurves, numSurf, numVol;
//		meshFile >> numPts >> numCurves >> numSurf >> numVol;
//
//		// Map of entityTag -> physicalTag for each entityDim (pt, curve...)
//		std::vector< std::map<int, int> > entities(4);
//		int entityTag, numPhysTag, physTag, numBounds;
//		double useless;
//		int parametric;
//
//		for (int i = 0; i < numPts; i++) { // Points
//			physTag = -1;
//			meshFile >> entityTag;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> numPhysTag;
//			std::cout << "nb physTag : " << numPhysTag << std::endl;
//			for (int j = 0; j < numPhysTag; j++) {
//				meshFile >> physTag;
//			}
//			entities[0][entityTag] = physTag;
//		}
//		for (int i = 0; i < numCurves; i++) { // Curves
//			physTag = -1;
//			meshFile >> entityTag;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> numPhysTag;
//			std::cout << "nb physTag : " << numPhysTag << std::endl;
//			for (int j = 0; j < numPhysTag; j++) {
//				meshFile >> physTag;
//			}
//			entities[1][entityTag] = physTag;
//			meshFile >> numBounds;
//			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
//				meshFile >> useless;
//			}
//		}
//		for (int i = 0; i < numSurf; i++) { // Surfaces
//			physTag = -1;
//			meshFile >> entityTag;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> numPhysTag;
//			std::cout << "nb physTag : " << numPhysTag << std::endl;
//			for (int j = 0; j < numPhysTag; j++) {
//				meshFile >> physTag;
//			}
//			entities[2][entityTag] = physTag;
//			meshFile >> numBounds;
//			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
//				meshFile >> useless;
//			}
//		}
//		for (int i = 0; i < numVol; i++) { // Volumes
//			physTag = -1;
//			meshFile >> entityTag;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> numPhysTag;
//			std::cout << "nb physTag : " << numPhysTag << std::endl;
//			for (int j = 0; j < numPhysTag; j++) {
//				meshFile >> physTag;
//			}
//			entities[3][entityTag] = physTag;
//			meshFile >> numBounds;
//			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
//				meshFile >> useless;
//			}
//		}
//
//		// 3) Reading partitioned entities
//		// -------------------------------
//		int numPartitions, numGhostEntities, ghostEntityTag, partition, partitionTag;
//		int numPtsPart, numCurvesPart, numSurfPart, numVolPart, parentDim, parentTag;
//		int entityTagPart;
//		std::vector< std::map<int, int> > entitiesPartition(4);
//		std::vector< std::map<int, std::vector<int> > > entitiesPartitionTag(4);
//
//		getline(meshFile, currentLine);
//		getline(meshFile, currentLine);
//		getline(meshFile, currentLine);
//		
//		meshFile >> numPartitions >> numGhostEntities;
//		for (int i = 0; i < numGhostEntities; ++i) {
//			meshFile >> ghostEntityTag;
//		}
//		meshFile >> numPtsPart >> numCurvesPart >> numSurfPart >> numVolPart;
//		
//		for (int i = 0; i < numPts; i++) { // Points
//			physTag = -1;
//			meshFile >> entityTagPart;
//			meshFile >> parentDim >> parentTag;
//			meshFile >> numPartitions;
//			for (int j = 0; j < numPartitions; ++j) {
//				meshFile >> partitionTag;
//				entitiesPartitionTag[0][entityTagPart].push_back(partitionTag);
//			}
//			meshFile >> useless >> useless >> useless;
//			meshFile >> numPhysTag;
//			for (int j = 0; j < numPhysTag; j++) {
//				meshFile >> physTag;
//			}
//			entitiesPartition[0][entityTagPart] = physTag;
//		}
//
//		for (int i = 0; i < numCurves; i++) { // Curves
//			physTag = -1;
//			meshFile >> entityTagPart;
//			meshFile >> parentDim >> parentTag;
//			meshFile >> numPartitions;
//			for (int j = 0; j < numPartitions; ++j) {
//				meshFile >> partitionTag;
//				entitiesPartitionTag[1][entityTagPart].push_back(partitionTag);
//			}
//			meshFile >> useless >> useless >> useless;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> numPhysTag;
//			for (int j = 0; j < numPhysTag; j++) {
//				meshFile >> physTag;
//			}
//			entitiesPartition[1][entityTagPart] = physTag;
//			meshFile >> numBounds;
//			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
//				meshFile >> useless;
//			}
//		}
//		for (int i = 0; i < numSurf; i++) { // Surfaces
//			physTag = -1;
//			meshFile >> entityTagPart;
//			meshFile >> parentDim >> parentTag;
//			meshFile >> numPartitions;
//			for (int j = 0; j < numPartitions; ++j) {
//				meshFile >> partitionTag;
//				entitiesPartitionTag[2][entityTagPart].push_back(partitionTag);
//			}
//			meshFile >> useless >> useless >> useless;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> numPhysTag;
//			for (int j = 0; j < numPhysTag; j++) {
//				meshFile >> physTag;
//			}
//			entitiesPartition[2][entityTagPart] = physTag;
//			meshFile >> numBounds;
//			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
//				meshFile >> useless;
//			}
//		}
//		for (int i = 0; i < numVol; i++) { // Volumes
//			physTag = -1;
//			meshFile >> entityTagPart;
//			meshFile >> parentDim >> parentTag;
//			meshFile >> numPartitions;
//			for (int j = 0; j < numPartitions; ++j) {
//				meshFile >> partitionTag;
//				entitiesPartitionTag[3][entityTagPart].push_back(partitionTag);
//			}
//			meshFile >> useless >> useless >> useless;
//			meshFile >> useless >> useless >> useless;
//			meshFile >> numPhysTag;
//			for (int j = 0; j < numPhysTag; j++) {
//				meshFile >> physTag;
//			}
//			entitiesPartition[3][entityTagPart] = physTag;
//			meshFile >> numBounds;
//			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
//				meshFile >> useless;
//			}
//		}
//		
//		getline(meshFile, currentLine);
//		getline(meshFile, currentLine);
//		getline(meshFile, currentLine);
//
//		// 4) Reading nodes
//		// ----------------
//		meshFile.ignore(1, '\n');
//		getline(meshFile, currentLine);
//		getline(meshFile, currentLine); // Skip keyword $EndPartitionedEntities
//		getline(meshFile, currentLine); // Skip keyword $Nodes
//
//		std::cout << "  2/Mesh nodes reading ..." << std::endl;
//
//		int numEntityBlocks(0), numNodes(0), minNodeTag(0), maxNodeTag(0);
//		int entityDim;
//		int numNodesInBlock;
//		double x, y, z;
//
//		meshFile >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag; // STOP
//		nodesGlobal = new Coord[numNodes];
//		for (int b = 0; b < numEntityBlocks; ++b) {
//			meshFile >> entityDim >> entityTag >> parametric >> numNodesInBlock;
//			std::vector<int> nodeTag(numNodesInBlock);
//			if (parametric != 0) throw ErrorECOGEN("parametric notation for nodes not compatible with ECOGEN", __FILE__, __LINE__);
//			for (int n = 0; n < numNodesInBlock; ++n) {
//				meshFile >> nodeTag[n];
//			}
//			for (int n = 0; n < numNodesInBlock; ++n) {
//				meshFile >> x >> y >> z;
//				m_nodes[nodeTag[n] - 1].setXYZ(x, y, z);
//			}
//		}
//	}
//	catch (ErrorECOGEN&) { throw; }
//}

//***********************************************************************

void MUSGmshV4::readMeshMonoCPU()
{
	// Re-open mesh file from now on 
	try {
		// 1) Opening mesh file at Gmsh format V4
		// --------------------------------------
		std::cout << "------------------------------------------------------" << std::endl;
		std::cout << " C) READING MESH FILE " + m_meshFile + " IN PROGRESS ..." << std::endl;
		std::ifstream meshFile(m_meshFile.c_str(), std::ios::in);
		if (!meshFile) { throw ErrorECOGEN("mesh file not found : " + m_meshFile, __FILE__, __LINE__); }
		std::string currentLine;

		// Skiping Gmsh version + go to entities
		// -----------------------------------
		getline(meshFile, currentLine); 
		getline(meshFile, currentLine); 
		getline(meshFile, currentLine); 
		getline(meshFile, currentLine);

		// 2) Reading entities
		// -------------------
		std::cout << "  1/Mesh entities reading ..." << std::endl;

		int numPts, numCurves, numSurf, numVol;
		meshFile >> numPts >> numCurves >> numSurf >> numVol;

		// Map of entityTag -> physicalTag for each entityDim (pt, curve...)
		std::vector< std::map<int, int> > entities(4);
		int entityTag, numPhysTag, physTag, numBounds;
		double useless;

		// Points
		int parametric;
		for (int i = 0; i < numPts; i++) {
			physTag = -1;
			meshFile >> entityTag;
			meshFile >> useless >> useless >> useless;
			meshFile >> numPhysTag;
			for (int j = 0; j < numPhysTag; j++) {
				meshFile >> physTag;
			}
			entities[0][entityTag] = physTag;
		}

		// Curves
		for (int i = 0; i < numCurves; i++) {
			physTag = -1;
			meshFile >> entityTag;
			meshFile >> useless >> useless >> useless;
			meshFile >> useless >> useless >> useless;
			meshFile >> numPhysTag;
			for (int j = 0; j < numPhysTag; j++) {
				meshFile >> physTag;
			}
			entities[1][entityTag] = physTag;
			meshFile >> numBounds;
			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
				meshFile >> useless;
			}
		}

		// Surfaces
		for (int i = 0; i < numSurf; i++) {
			physTag = -1; 
			meshFile >> entityTag;
			meshFile >> useless >> useless >> useless;
			meshFile >> useless >> useless >> useless;
			meshFile >> numPhysTag;
			for (int j = 0; j < numPhysTag; j++) {
				meshFile >> physTag;
			}
			entities[2][entityTag] = physTag;
			meshFile >> numBounds;
			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
				meshFile >> useless;
			}
		}

		// Volumes
		for (int i = 0; i < numVol; i++) {
			physTag = -1;
			meshFile >> entityTag;
			meshFile >> useless >> useless >> useless;
			meshFile >> useless >> useless >> useless;
			meshFile >> numPhysTag;
			for (int j = 0; j < numPhysTag; j++) {
				meshFile >> physTag;
			}
			entities[3][entityTag] = physTag;
			meshFile >> numBounds;
			for (int k = 0; k < numBounds; k++) { // Ignore boundary tags
				meshFile >> useless;
			}
		}

		// 2) Reading nodes
		// ----------------
		meshFile.ignore(1, '\n');
		getline(meshFile, currentLine);
		getline(meshFile, currentLine); // Skip keyword $EndEntities
		getline(meshFile, currentLine); // Skip keyword $Nodes

		int entityDim;
		std::cout << "  2/Mesh nodes reading ..." << std::endl;
		
		int numEntityBlocks(0), numNodes(0), minNodeTag(0), maxNodeTag(0);
		int numNodesInBlock;
		double x, y, z;

		meshFile >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;
		
		if (maxNodeTag != numNodes) throw ErrorECOGEN("sparse node tags not compatible of ECOGEN for mesh fle : " + m_meshFile, __FILE__, __LINE__); 
		m_numberNodes = numNodes;
		m_nodes = new Coord[m_numberNodes];
		for (int b = 0; b < numEntityBlocks; ++b) {
			meshFile >> entityDim >> entityTag >> parametric >> numNodesInBlock; 
			std::vector<int> nodeTag(numNodesInBlock);
			if (parametric != 0) throw ErrorECOGEN("parametric notation for nodes not compatible with ECOGEN", __FILE__, __LINE__);
			for (int n = 0; n < numNodesInBlock; ++n) {
				meshFile >> nodeTag[n];
			}
			for (int n = 0; n < numNodesInBlock; ++n) {
					meshFile >> x >> y >> z;
					m_nodes[nodeTag[n]-1].setXYZ(x, y, z);
			}
		}

		// 3) Reading elements
		// -------------------
		meshFile.ignore(1, '\n');
		getline(meshFile, currentLine); // Skip keyword $EndNodes 
		getline(meshFile, currentLine); // Skip keyword $Elements
		std::cout << "  3/0D/1D/2D/3D elements reading ..." << std::endl;

		int numElts(0), minEltTag(0), maxEltTag(0);
		int numEltsInBlock, eltType;
		int idElt(0);

		meshFile >> numEntityBlocks >> numElts >> minEltTag >> maxEltTag;
		if ( maxEltTag != numElts ) throw ErrorECOGEN("sparse element tags not compatible of ECOGEN for mesh fle : " + m_meshFile, __FILE__, __LINE__);
		m_numberElements = numElts;
		// Allocation array of elements
		m_elements = new ElementNS*[m_numberElements];
		// Reading elements and assigning geometric properties 
		m_numberElements1D = 0, m_numberElements2D = 0, m_numberElements3D = 0;
		
		for (int i = 0; i < numEntityBlocks; ++i) {
			meshFile >> entityDim >> entityTag >> eltType >> numEltsInBlock;
			for (int j = 0; j < numEltsInBlock; ++j) {
				this->readElement(m_nodes, meshFile, &m_elements[idElt], entities, entityDim, entityTag, eltType);
				if (m_elements[idElt]->getTypeGmsh() == 15) { m_numberElements0D++; }
				else if (m_elements[idElt]->getTypeGmsh() == 1) { m_numberElements1D++; }
				else if (m_elements[idElt]->getTypeGmsh() <= 3) { m_numberElements2D++; m_totalSurface += m_elements[idElt]->getVolume(); }
				else if (m_elements[idElt]->getTypeGmsh() <= 7) { m_numberElements3D++; m_totalVolume += m_elements[idElt]->getVolume(); }
				else { throw ErrorECOGEN("Element type of .msh not handled by ECOGEN", __FILE__, __LINE__); }
				idElt++;
			}
			m_numberInnerElements = m_numberElements;
		}
	}
	catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MUSGmshV4::readElement(const Coord* nodesTable, std::ifstream& meshFile, ElementNS** element, std::vector< std::map<int,int> >& entities, int entityDim, int entityTag, int eltType)
{
	int numberElement, numberPhysicEntity, numberGeometricEntity;
	meshFile >> numberElement;
	
	// Assigning tags to element thanks to entityDim, entityTag correspondance
	numberGeometricEntity = entityTag;
	std::map<int, int>::iterator it = entities[entityDim].find(entityTag);
	
	if (it != entities[entityDim].end()) 
		numberPhysicEntity = entities[entityDim][entityTag];
	else 
		Errors::errorMessage("Specific element entity dimension of file .msh unknown for ECOGEN");

	// 1) Assignation of the number of vertices according to element type
	// ------------------------------------------------------------------
	switch (eltType)
	{
	case 1: // Segment (two points)
		*element = new ElementSegment;
		m_numberSegments++;
		break;
	case 2: // Triangle (three points)
		*element = new ElementTriangle;
		m_numberTriangles++;
		break;
	case 3: // Quadrangle (four points)
		*element = new ElementQuadrangle;
		m_numberQuadrangles++;
		break;
	case 4: // Tetrahedron (four points)
		*element = new ElementTetrahedron;
		m_numberTetrahedrons++;
		break;
		//case 7: // Quadrangular pyramid (five points) // This element seems to not work with Gmsh, volumes of elements seem to cause this issue 
		//  *element = new ElementPyramid; 
		//  m_numberPyramids++;
		//  break;
	case 15: // Point (a vertex)
		*element = new ElementPoint;
		m_numberPoints++;
		break;
	case 5: // Hexahedron (eight points)
		*element = new ElementHexahedron;
		m_numberHexahedrons++;
		break;
	case 6: // Prism (six points)
		*element = new ElementPrism;
		m_numberHexahedrons++;
		break;
	default:
		Errors::errorMessage("Element type of file .msh inknown for ECOGEN");
		break;
	}

	// 3) Building the element and its properties
	// ------------------------------------------
	int currentNode(0);
	int* numNode = new int[(*element)->getNumberNoeuds()];
	Coord* node = new Coord[(*element)->getNumberNoeuds()];
	for (int i = 0; i < (*element)->getNumberNoeuds(); i++)
	{
		meshFile >> currentNode;
		numNode[i] = currentNode - 1; // Offset because array start at 0
		node[i] = nodesTable[currentNode - 1];
	}
	int indexElement(numberElement - 1);
	(*element)->construitElement(numNode, node, numberPhysicEntity, numberGeometricEntity, indexElement);
	
	delete[] node;
	delete[] numNode;
}