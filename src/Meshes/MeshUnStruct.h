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
    ~MeshUnStruct() override;

    static std::string readMeshFileExtension(const std::string& meshFile);

    // --- Mesh virtual member functions ---
    void assignLimits(std::vector<BoundCond*>& boundCond) override;
    int initializeGeometrie(TypeMeshContainer<Cell*>& cells,
                            TypeMeshContainer<Cell*>& cellsGhost,
                            TypeMeshContainer<CellInterface*>& cellInterfaces,
                            const int& /*resumeSimulation*/,
                            bool pretraitementParallele = true,
                            SimulationOrder schemeOrder = SimulationOrder::first) override;
    std::string whoAmI() const override { return Errors::defaultString; };

    // --- MeshUnStruct virtual member functions ---
    //! \brief     initialize the geometry for single CPU computation
    //! \param     cells
    //! \param     cellInterfaces
    //! \param     schemeOrder     scheme order (currently only first order working)
    virtual void initGeometryMonoCPU(TypeMeshContainer<Cell*>& cells,
                                     TypeMeshContainer<CellInterface*>& cellInterfaces,
                                     SimulationOrder schemeOrder = SimulationOrder::first) = 0;
    //! \brief     initialize the geometry for multi CPUs computation
    //! \param     cells
    //! \param     cellInterfaces
    //! \param     schemeOrder     scheme order (currently only first order working)
    virtual void initGeometryParallel(TypeMeshContainer<Cell*>& cells,
                                      TypeMeshContainer<Cell*>& cellsGhost,
                                      TypeMeshContainer<CellInterface*>& cellInterfaces,
                                      SimulationOrder schemeOrder = SimulationOrder::first) = 0;
    //! \brief     split original mesh file for computation on several CPUs
    virtual void preProcessMeshFileForParallel()                                                              = 0;
    void initCpuMeshSequential(TypeMeshContainer<Cell*>& cells, SimulationOrder& schemeOrder) override        = 0;
    void initCpuMeshParallel(TypeMeshContainer<Cell*>& cells, SimulationOrder& schemeOrder, int cpu) override = 0;

    // Printing / Reading
    //! \brief    write monocpu mesh information
    void writeMeshInfoData() const;
    void writeHeaderPiece(std::ofstream& fileStream, TypeMeshContainer<Cell*>* /*cellsLvl*/) const override;
    size_t getNodes(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    size_t getConnectivity(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    size_t getConnectivity(std::vector<int64_t>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    size_t getConnectivity(std::vector<int32_t>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getOffsets(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getOffsets(std::vector<int64_t>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getOffsets(std::vector<int32_t>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getTypeCell(std::vector<double>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getTypeCell(std::vector<unsigned char>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getData(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, const int var, int phase) const override;
    void getData(TypeMeshContainer<Cell*>* cellsLvl, std::vector<float>& dataset, const int var, int phase) const override;

    //void getData(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, const int var, int phase) const override;
    void setDataSet(std::vector<double>& dataset, TypeMeshContainer<Cell*>* cellsLvl, const int var, int phase) const override;
    void extractAbsVelocityMRF(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, Source* sourceMRF) const override;
    void extractAbsVelocityMRF(TypeMeshContainer<Cell*>* cellsLvl, std::vector<float>& dataset, Source* sourceMRF) const override;
    void extractReferenceLength(std::vector<Cell*>* cellsLvl, std::vector<double>& dataset) const override;
    void extractReferenceLength(std::vector<Cell*>* cellsLvl, std::vector<float>& dataset) const override;

  protected:
    std::string m_meshFile; //!< Name of the mesh file read
    std::string m_nameMesh; //!< Name of the mesh file without extension

    int m_numberNodes;                 //!< Number of nodes of the geometric domain
    int m_numberInnerNodes;            //!< Number of inner nodes (except from ghosts)
    Coord* m_nodes;                    //!< Array of node coordinates in the geometric domain
    int m_numberInnerElements;         //!< Number of elements of n dimension of internal compute
    int m_numberGhostElements;         //!< Number of ghost elements of dimension n for parallel computation
    int m_numberCommunicatingElements; //!< Real number of communicating elements
    ElementNS** m_elements;            //!< Array of internal geometric elements
    FaceNS** m_faces;                  //!< Array of geometrical faces
    std::vector<BoundCond*> m_bound;   //!< Array of boundary conditions

    int m_numberInnerFaces;    //!< Number of faces between two cells of compute
    int m_numberBoundFaces;    //!< Number of faces between a compute cell and a boundary
    int m_numberFacesParallel; //!< Number of faces between a compute cell and a ghost cell

    int m_numberGhostCells; //!< Number of ghost cells

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
    double m_totalSurface; //!< Sum of 2D element surfaces
    double m_totalVolume;  //!< Sum of 3D element volumes

  private:
    template <typename T> size_t getConnectivity_private(std::vector<T>& dataset) const;
    template <typename T> void getOffsets_private(std::vector<T>& dataset) const;
    template <typename T> void getTypeCell_private(std::vector<T>& dataset) const;
    template <typename T> void getData_private(TypeMeshContainer<Cell*>* cellsLvl, std::vector<T>& dataset, const int var, int phase) const;
    template <typename T> void extractAbsVelocityMRF_private(TypeMeshContainer<Cell*>* cellsLvl, std::vector<T>& dataset, Source* sourceMRF) const;
    template <typename T> void extractReferenceLength_private(std::vector<Cell*>* cellsLvl, std::vector<T>& dataset) const;
};

#endif // MESHUNSTRUCT_H
