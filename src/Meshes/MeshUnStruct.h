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

#ifndef MESHUNSTRUCT_H
#define MESHUNSTRUCT_H

#include "Mesh.h"
#include "MeshUnStruct/MUSGmsh/HeaderElements.h"
#include "../Order2/CellInterfaceO2NS.h"
#include "../Order2/CellO2NS.h"
#include "../InputOutput/IO.h"

class MeshUnStruct : public Mesh
{
public:
  MeshUnStruct(const std::string& meshFile, const std::string& meshExtension);
  virtual ~MeshUnStruct();

  static std::string readMeshFileExtension(const std::string& meshFile);

  // --- Mesh virtual member functions ---
  virtual void assignLimits(std::vector<BoundCond*>& boundCond);
  virtual int initializeGeometrie(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces,
    const int& /*restartSimulation*/, bool pretraitementParallele = true, std::string ordreCalcul = "FIRSTORDER");
  virtual std::string whoAmI() const { return 0; };

  // --- MeshUnStruct virtual member functions --- 
  //! \brief     initialize the geometry for single CPU computation 
  //! \param     cells            
  //! \param     cellInterfaces   
  //! \param     computeOrder     scheme order (currently only first order working) 
  virtual void initGeometryMonoCPU(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string computeOrder = "FIRSTORDER") = 0;
  //! \brief     initialize the geometry for multi CPUs computation
  //! \param     cells            
  //! \param     cellInterfaces   
  //! \param     computeOrder     scheme order (currently only first order working)
  virtual void initGeometryParallel(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string computeOrder = "FIRSTORDER") = 0;
  //! \brief     split original mesh file for computation on several CPUs
  virtual void preProcessMeshFileForParallel() = 0;
  virtual void initCpuMeshSequential(TypeMeshContainer<Cell*>& cells, std::string &computeOrder) = 0;
  virtual void initCpuMeshParallel(TypeMeshContainer<Cell*>& cells, std::string &computeOrder, int cpu) = 0;

  // Printing / Reading
  //! \brief    write monocpu mesh information
  void writeMeshInfoData() const;
  virtual void writeHeaderPiece(std::ofstream& fileStream, TypeMeshContainer<Cell*>* /*cellsLvl*/) const;
  virtual void getNodes(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const;
  virtual void getConnectivity(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const;
  virtual void getOffsets(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const;
  virtual void getTypeCell(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const;
  virtual void getData(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, const int var, int phase) const;
  virtual void setDataSet(std::vector<double>& dataset, TypeMeshContainer<Cell*>* cellsLvl, const int var, int phase) const;
  virtual void extractAbsVelocityMRF(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, Source* sourceMRF) const;
  virtual void extractReferenceLength(std::vector<Cell*>* cellsLvl, std::vector<double>& dataset) const;

protected:
  std::string m_meshFile;  //!< Name of the mesh file read
  std::string m_nameMesh;  //!< Name of the mesh file without extension

  int m_numberNodes;                  //!< Number of nodes of the geometric domain
  int m_numberInnerNodes;             //!< Number of inner nodes (except from ghosts)
  Coord* m_nodes;                     //!< Array of node coordinates in the geometric domain
  int m_numberInnerElements;          //!< Number of elements of n dimension of internal compute
  int m_numberGhostElements;          //!< Number of ghost elements of dimension n for parallel computation
  int m_numberCommunicatingElements;  //!< Real number of communicating elements
  ElementNS** m_elements;             //!< Array of internal geometric elements
  FaceNS** m_faces;                   //!< Array of geometrical faces
  std::vector<BoundCond*> m_bound;    //!< Array of boundary conditions

  int m_numberInnerFaces;          //!< Number of faces between two cells of compute
  int m_numberBoundFaces;          //!< Number of faces between a compute cell and a boundary
  int m_numberFacesParallel;       //!< Number of faces between a compute cell and a ghost cell

  int m_numberGhostCells;           //!< Number of ghost cells

  int m_numberElements0D;
  int m_numberElements1D;
  int m_numberElements2D;
  int m_numberElements3D;
  int m_numberSegments;
  int m_numberTriangles;
  int m_numberQuadrangles;
  int m_numberTetrahedrons;
  int m_numberPyramids;
  int m_numberPoints;
  int m_numberHexahedrons;

  //statistics
  double m_totalSurface;    //!< Sum of 2D element surfaces
  double m_totalVolume;     //!< Sum of 3D element volumes
};

#endif // MESHUNSTRUCT_H