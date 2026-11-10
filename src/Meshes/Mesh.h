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

#ifndef MESH_H
#define MESH_H

#include <ctime>
#include <fstream>
#include "../libTierces/tinyxml2.h"
#include "../Order1/Cell.h"
#include "../Order1/CellGhost.h"
#include "../Order1/CellInterface.h"
#include "../Order2/CellO2GhostCartesian.h"
#include "../Order2/CellO2GhostNS.h"
#include "../BoundConds/HeaderBoundCond.h"
#include "../Maths/Coord.h"
#include "../Parallel/Parallel.h"
#include "../AdditionalPhysics/HeaderQuantitiesAddPhys.h"
#include "../Maths/GeometricObject.h"
#include "../Sources/HeaderSources.h"

#include "../InputOutput/IntTypes.h" // Use same integer types than those used to write output files

//! \class     Mesh
//! \brief     Abstract class for a mesh
//! \details   Can not be instanciated, depend on the mesh properties
class Mesh
{
  public:
    Mesh();
    virtual ~Mesh();

    virtual void assignLimits(std::vector<BoundCond*>& boundCond)                         = 0;
    virtual int initializeGeometrie(TypeMeshContainer<Cell*>& cells,
                                    TypeMeshContainer<Cell*>& cellsGhost,
                                    TypeMeshContainer<CellInterface*>& cellInterfaces,
                                    const int& resumeSimulation,
                                    bool pretraitementParallele = true,
                                    SimulationOrder schemeOrder = SimulationOrder::first) = 0; //!< renvoi le number de dimensions (1,2 ou 3)
    virtual std::string whoAmI() const
    {
      Errors::errorMessage("whoAmI not available for considered mesh");
      return std::string("");
    };
    virtual void setImmersedBoundaries(TypeMeshContainer<CellInterface*>* /*cellInterfacesLvl*/, SimulationOrder /*schemeOrder*/) const {};

    //Accessors
    //---------
    const int& getProblemDimension() const { return m_problemDimension; };
    const int& getNumberCells() const { return m_numberCellsCalcul; };
    const int& getNumberCellsTotal() const { return m_numberCellsTotal; };
    const int& getNumberFaces() const { return m_numberFacesTotal; };
    const int& getNumFichier() const { return m_numFile; };
    const TypeM& getType() const { return m_type; };
    virtual int getNumberCellsX() const { return 0; };
    virtual int getNumberCellsY() const { return 0; };
    virtual int getNumberCellsZ() const { return 0; };
    virtual int getNumberCellsXGlobal() const { return 0; };
    virtual int getNumberCellsYGlobal() const { return 0; };
    virtual int getNumberCellsZGlobal() const { return 0; };
    virtual int getOffsetX() const { return 0; };
    virtual int getOffsetY() const { return 0; };
    virtual int getOffsetZ() const { return 0; };
    virtual int getNumberCpuX() const { return 0; };
    virtual int getNumberCpuY() const { return 0; };
    virtual int getNumberCpuZ() const { return 0; };

    virtual int getLvlMax() const { return 0; };

    //Printing
    //--------
    void writeResultsGnuplot(std::vector<Cell*>* cellsLvl, std::ofstream& fileStream, GeometricObject* objet = 0, bool recordPsat = false) const;
    virtual void writeHeaderPiece(std::ofstream& /*fileStream*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("writeHeaderPiece not available for considered mesh");
    };
    virtual std::string getStringExtent(bool /*global*/ = false) const
    {
      Errors::errorMessage("getStringExtent not available for considered mesh");
      return std::string("");
    };
    //! \brief For structured grid, get global grid extent under the form of an array
    virtual std::array<int, 6> getArrayGlobalExtent() const
    {
      Errors::errorMessage("getArrayGlobalExtent not available for considered mesh");
      return {0, 0, 0, 0, 0, 0};
    };

    //! \brief For structured grid, get coordinates along a given axis
    virtual void getCoord(std::vector<double>& /*dataset*/, Axis /*axis*/, bool global = false) const
    {
      (void)global; // Parameter is intentionally unused: suppresses warning
      Errors::errorMessage("getCoord not available for considered mesh");
    };
    //! \brief Store coordinates of mesh nodes in a flat vector
    //! \return Node count if success, 0 in case of overflow
    //! \note Check that the node count doesn't overflow the maximal value representable by a 'size_t'
    virtual size_t getNodes(std::vector<double>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("getNodes not available for considered mesh");
      return 0;
    };
    //! \brief Store connectivities of mesh cells in a flat vector of double
    //! \return Cell count if success, 0 in case of overflow
    //! \note Check that the cell count doesn't overflow the maximal value representable by a 'size_t'.
    //! Don't check that the connectivity values fit in 'double' variables.
    //!
    //! \remark Method has been kept to not break legacy I/Os but avoiding casts is probably better
    virtual size_t getConnectivity(std::vector<double>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("getConnectivity not available for considered mesh");
      return 0;
    };
    //! \brief Store connectivities of mesh cells in a flat vector of 'ECOGEN_VTKHDF_intType' integers.
    //! \return Cell count if success, 0 in case of overflow
    //! \note Check that the cell count doesn't overflow the maximal value representable by a 'size_t'.
    //! Don't check that the connectivity values fit the vector type.
    virtual size_t getConnectivity(std::vector<int64_t>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("getConnectivity not available for considered mesh and storage in 'int' vector");
      return 0;
    };
    //! \brief Store connectivities of mesh cells in a flat vector of 'ECOGEN_VTKHDF_intType' integers.
    //! \return Cell count if success, 0 in case of overflow
    //! \note Check that the cell count doesn't overflow the maximal value representable by a 'size_t'.
    //! Don't check that the connectivity values fit the vector type.
    virtual size_t getConnectivity(std::vector<int32_t>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("getConnectivity not available for considered mesh and storage in 'int' vector");
      return 0;
    };

    //! \brief Store offsets in connectivity array of mesh cells
    //! \note Check that the max offset doesn't overflow the maximal integer exactly representable by a 'double'.
    //!
    //! \remark Method has been kept to not break legacy I/Os but avoiding casts is probably better
    virtual void getOffsets(std::vector<double>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("getOffsets not available for considered mesh");
    };
    //! \brief Store offsets in connectivity array of mesh cells
    //! \note Check that the max offset doesn't overflow the maximal value representable by a int64_t variable.
    virtual void getOffsets(std::vector<int64_t>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("getOffsets not available for considered mesh and storage in 'int64_t' vector");
    };
    //! \brief Store offsets in connectivity array of mesh cells
    //! \note Check that the max offset doesn't overflow the maximal value representable by a int32_t variable.
    virtual void getOffsets(std::vector<int32_t>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("getOffsets not available for considered mesh and storage in 'int32_t' vector");
    };

    virtual void getTypeCell(std::vector<double>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("getTypeCell not available for considered mesh");
    };
    virtual void getTypeCell(std::vector<unsigned char>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/) const
    {
      Errors::errorMessage("'getTypeCell' not available for considered mesh and storage in 'unsigned char' vector");
    };

    //! \brief     Extracting data for printing results
    //! \details   This method enable to extract a set of data for mixture or phase, scalar or vetor
    //! \param     cellsLvl      data structure containing pointer to cells
    //! \param     dataset       double vector containing the extracted data
    //! \param     var           number of requested varaible to extract (>0 for scalar, <0 for vector)
    //! \param     phase         number of requested phase (-1 for mixture, -2 for transport, -3 for xi, -4 for gradient density mixture)
    virtual void getData(std::vector<Cell*>* /*cellsLvl*/, std::vector<double>& /*dataset*/, const int /*var*/, int /*phase*/) const
    {
      Errors::errorMessage("getData not available for considered mesh");
    };
    virtual void getData(std::vector<Cell*>* /*cellsLvl*/, std::vector<float>& /*dataset*/, const int /*var*/, int /*phase*/) const
    {
      Errors::errorMessage("getData not available for considered mesh");
    };

    //! \brief     Extracting data for printing results
    //! \details   This method enable to extract a set of data for mixture or phase, scalar or vetor
    //! \param     cellsLvl      data structure containing pointer to cells
    //! \param     dataset       double vector containing the extracted data
    //! \param     var           number of requested varaible to extract (>0 for scalar, <0 for vector)
    //! \param     phase         number of requested phase (-1 for mixture, -2 for transport, -3 for xi, -4 for gradient density mixture)
    virtual void setDataSet(std::vector<double>& /*dataset*/, std::vector<Cell*>* /*cellsLvl*/, const int /*var*/, int /*phase*/) const
    {
      Errors::errorMessage("setDataSet not available for requested mesh");
    };
    virtual void refineCellAndCellInterfaces(Cell* /*cell*/, const std::vector<AddPhys*>& /*addPhys*/, int& /*nbCellsTotalAMR*/)
    {
      Errors::errorMessage("refineCellAndCellInterfaces not available for requested mesh");
    };
    //! \brief     Extracting absolute velocity for specific Moving Reference Frame computations
    //! \param     cellsLvl         data structure containing pointer to cells
    //! \param     sourceMRF        pointer to the corresponding MRF source
    //! \param     dataset       double vector containing the extracted data
    virtual void extractAbsVelocityMRF(std::vector<Cell*>* /*cellsLvl*/, std::vector<double>& /*dataset*/, Source* /*sourceMRF*/) const
    {
      Errors::errorMessage("extractAbsVeloxityMRF not available for considered mesh");
    };
    virtual void extractAbsVelocityMRF(std::vector<Cell*>* /*cellsLvl*/, std::vector<float>& /*dataset*/, Source* /*sourceMRF*/) const
    {
      Errors::errorMessage("extractAbsVeloxityMRF not available for considered mesh");
    };

    virtual void extractReferenceLength(std::vector<Cell*>* /*cellsLvl*/, std::vector<double>& /*dataset*/) const
    {
      Errors::errorMessage("extractReferenceLength not available for considered mesh");
    };
    virtual void extractReferenceLength(std::vector<Cell*>* /*cellsLvl*/, std::vector<float>& /*dataset*/) const
    {
      Errors::errorMessage("extractReferenceLength not available for considered mesh");
    };

    virtual void printDomainDecomposition(std::ofstream& /*fileStream*/) {};
    virtual void readDomainDecomposition(std::ifstream& /*fileStream*/) {};

    //Specific to AMR method
    //----------------------
    virtual void procedureRaffinementInitialization(std::vector<Cell*>* /*cellsLvl*/,
                                                    TypeMeshContainer<Cell*>* /*cellsLvlGhost*/,
                                                    std::vector<CellInterface*>* /*cellInterfacesLvl*/,
                                                    const std::vector<AddPhys*>& /*addPhys*/,
                                                    int& nbCellsTotalAMR,
                                                    std::vector<GeometricalDomain*>& /*domains*/,
                                                    Eos** /*eos*/,
                                                    const int& /*resumeSimulation*/,
                                                    SimulationOrder /*schemeOrder*/,
                                                    std::vector<GeometricalDomain*>& /*solidDomains*/)
    {
      nbCellsTotalAMR = m_numberCellsCalcul;
    };
    virtual void procedureRaffinement(std::vector<Cell*>* /*cellsLvl*/,
                                      TypeMeshContainer<Cell*>* /*cellsLvlGhost*/,
                                      std::vector<CellInterface*>* /*cellInterfacesLvl*/,
                                      const int& /*lvl*/,
                                      const std::vector<AddPhys*>& /*addPhys*/,
                                      int& /*nbCellsTotalAMR*/,
                                      Eos** /*eos*/) {};

    //Specific for parallel
    //---------------------
    virtual void initializePersistentCommunications(const TypeMeshContainer<Cell*>& cells, SimulationOrder schemeOrder);
    virtual void finalizeParallele(const int& lvlMax);
    virtual void parallelLoadBalancingAMR(std::vector<Cell*>* /*cellsLvl*/,
                                          TypeMeshContainer<Cell*>* /*cellsLvlGhost*/,
                                          std::vector<CellInterface*>* /*cellInterfacesLvl*/,
                                          SimulationOrder /*schemeOrder*/,
                                          const std::vector<AddPhys*>& /*addPhys*/,
                                          Eos** /*eos*/,
                                          int& /*nbCellsTotalAMR*/,
                                          std::vector<GeometricalDomain*>& /*solidDomains*/,
                                          bool /*init*/ = false) {};

    //Specific for mesh mapping restart
    //---------------------------------
    virtual std::string getMeshExtension() const
    {
      Errors::errorMessage("getMeshExtension not available for requested mesh");
      return std::string("");
    };
    //! \brief     Initialize mesh of a single partition for restart with mesh mapping option
    //! \details   This mesh object has only the elements and nodes filled
    virtual void initCpuMeshSequential(TypeMeshContainer<Cell*>& /*cells*/, SimulationOrder& /*schemeOrder*/)
    {
      Errors::errorMessage("initCpuMeshSequential not available for requested mesh");
    }
    //! \brief     Initialize mesh of a single partition of a partionned mesh for restart with mesh mapping option
    //! \details   This mesh object has only the elements and nodes filled
    virtual void initCpuMeshParallel(TypeMeshContainer<Cell*>& /*cells*/, SimulationOrder& /*schemeOrder*/, int /*cpu*/)
    {
      Errors::errorMessage("initCpuMeshParallel not available for requested mesh");
    }

  protected:
    mutable int m_numFile;

    int m_problemDimension; /*indicator 1D/2D/3D*/
    int m_numberElements; /*Total element count (internal computing cells of dimension n + boundary elements of dimension n-1 + ghost cells of dimension n)*/
    int m_numberFacesTotal;  /*Count of faces at cells interfaces or between cells and boundaries */
    int m_numberCellsCalcul; /*Count of computing cells (domain-internal cells) at level 0*/
    int m_numberCellsTotal;  /*Count of internal + ghosts cells (for parallel communications) at level 0*/

    TypeM m_type;

    //! \brief Access data of cell \a cell and append the wanted variable \a cell in \a dataset vector.
    //! \tparam T Type of dataset vector.
    //! \param cellsLvl Cell array.
    //! \param dataset Vector to fill with cell array variable.
    //! \param var parameter to choose to get scalar (if strictly positive) or vector (in lower or ewal to 0) variable.
    //! \param phase identifier of variable to store.
    //! \param isAMR true if called from AMR.
    template <typename T> void getData_body(Cell* cell, std::vector<T>& dataset, const int var, const int phase, const bool isAMR) const;
};
#endif // MESH_H
