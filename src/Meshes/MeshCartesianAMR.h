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

#ifndef MESHCARTESIANAMR_H
#define MESHCARTESIANAMR_H

#include "MeshCartesian.h"

class MeshCartesianAMR : public MeshCartesian
{
  public:
    MeshCartesianAMR(double lX,
                     int numberCellsX,
                     double lY,
                     int numberCellsY,
                     double lZ,
                     int numberCellsZ,
                     std::vector<stretchZone> stretchX,
                     std::vector<stretchZone> stretchY,
                     std::vector<stretchZone> stretchZ,
                     int lvlMax         = 0,
                     double criteriaVar = 1.e10,
                     bool varRho        = false,
                     bool varP          = false,
                     bool varU          = false,
                     bool varAlpha      = false,
                     double xiSplit     = 1.,
                     double xiJoin      = 1.);
    ~MeshCartesianAMR() override;

    int initializeGeometrie(TypeMeshContainer<Cell*>& cells,
                            TypeMeshContainer<Cell*>& cellsGhost,
                            TypeMeshContainer<CellInterface*>& cellInterfaces,
                            const int& resumeSimulation,
                            bool /*pretraitementParallele*/,
                            SimulationOrder schemeOrder) override;
    void initializeGeometrieAMR(TypeMeshContainer<Cell*>& cells,
                                TypeMeshContainer<Cell*>& cellsGhost,
                                TypeMeshContainer<CellInterface*>& cellInterfaces,
                                const int& resumeSimulation,
                                SimulationOrder schemeOrder);
    void assignElementProperties(TypeMeshContainer<Cell*>& cells, std::vector<decomposition::Key<3>>& keys);
    void createCellInterfacesFacesAndGhostCells(TypeMeshContainer<Cell*>& cells,
                                                TypeMeshContainer<Cell*>& cellsGhost,
                                                TypeMeshContainer<CellInterface*>& cellInterfaces,
                                                SimulationOrder schemeOrder);
    void createGhostCell(const TypeMeshContainer<Cell*>& cells,
                         TypeMeshContainer<Cell*>& cellsGhost,
                         const SimulationOrder schemeOrder,
                         const decomposition::Key<3> nKey,
                         const int neighbour,
                         const int cell_idx);
    void procedureRaffinementInitialization(TypeMeshContainer<Cell*>* cellsLvl,
                                            TypeMeshContainer<Cell*>* cellsLvlGhost,
                                            TypeMeshContainer<CellInterface*>* cellInterfacesLvl,
                                            const std::vector<AddPhys*>& addPhys,
                                            int& nbCellsTotalAMR,
                                            std::vector<GeometricalDomain*>& domains,
                                            Eos** eos,
                                            const int& resumeSimulation,
                                            SimulationOrder schemeOrder,
                                            std::vector<GeometricalDomain*>& solidDomains) override;
    void procedureRaffinement(TypeMeshContainer<Cell*>* cellsLvl,
                              TypeMeshContainer<Cell*>* cellsLvlGhost,
                              TypeMeshContainer<CellInterface*>* cellInterfacesLvl,
                              const int& lvl,
                              const std::vector<AddPhys*>& addPhys,
                              int& nbCellsTotalAMR,
                              Eos** eos) override;
    std::string whoAmI() const override;

    //Printing / Reading
    void writeHeaderPiece(std::ofstream& fileStream, TypeMeshContainer<Cell*>* cellsLvl) const override;
    size_t getNodes(std::vector<double>& dataset, std::vector<Cell*>* cellsLvl) const override;
    size_t getConnectivity(std::vector<double>& dataset, std::vector<Cell*>* cellsLvl) const override;
    size_t getConnectivity(std::vector<int64_t>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    size_t getConnectivity(std::vector<int32_t>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getOffsets(std::vector<double>& dataset, std::vector<Cell*>* cellsLvl) const override;
    void getOffsets(std::vector<int64_t>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getOffsets(std::vector<int32_t>& dataset, std::vector<Cell*>* /*cellsLvl*/) const override;
    void getTypeCell(std::vector<double>& dataset, std::vector<Cell*>* cellsLvl) const override;
    void getTypeCell(std::vector<unsigned char>& dataset, std::vector<Cell*>* cellsLvl) const override;
    void getData(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, const int var, int phase) const override;
    void getData(TypeMeshContainer<Cell*>* cellsLvl, std::vector<float>& dataset, const int var, int phase) const override;
    void setDataSet(std::vector<double>& dataset, TypeMeshContainer<Cell*>* cellsLvl, const int var, int phase) const override;
    void refineCellAndCellInterfaces(Cell* cell, const std::vector<AddPhys*>& addPhys, int& nbCellsTotalAMR) override;
    void printDomainDecomposition(std::ofstream& fileStream) override;
    void readDomainDecomposition(std::ifstream& fileStream) override;

    //Accesseurs
    int getLvlMax() const override { return m_lvlMax; };

    //Pour parallele
    void initializePersistentCommunications(const TypeMeshContainer<Cell*>& cells, SimulationOrder schemeOrder) override;
    void finalizeParallele(const int& lvlMax) override;
    void parallelLoadBalancingAMR(TypeMeshContainer<Cell*>* cellsLvl,
                                  TypeMeshContainer<Cell*>* cellsLvlGhost,
                                  TypeMeshContainer<CellInterface*>* cellInterfacesLvl,
                                  SimulationOrder schemeOrder,
                                  const std::vector<AddPhys*>& addPhys,
                                  Eos** eos,
                                  int& nbCellsTotalAMR,
                                  std::vector<GeometricalDomain*>& solidDomains,
                                  bool init = false) override;
    virtual void computePotentialBalancing(TypeMeshContainer<Cell*>* cellsLvl,
                                           bool init,
                                           int lvl,
                                           bool& balance,
                                           std::vector<typename decomposition::Key<3>::value_type>& indicesSendStartGlobal,
                                           std::vector<typename decomposition::Key<3>::value_type>& indicesSendEndGlobal,
                                           std::vector<typename decomposition::Key<3>::value_type>& indicesReceiveStartGlobal,
                                           std::vector<typename decomposition::Key<3>::value_type>& indicesReceiveEndGlobal);
    virtual void balance(TypeMeshContainer<Cell*>* cellsLvl,
                         TypeMeshContainer<Cell*>* cellsLvlGhost,
                         TypeMeshContainer<CellInterface*>* cellInterfacesLvl,
                         SimulationOrder schemeOrder,
                         const std::vector<AddPhys*>& addPhys,
                         Eos** eos,
                         int& nbCellsTotalAMR,
                         std::vector<typename decomposition::Key<3>::value_type>& indicesSendStartGlobal,
                         std::vector<typename decomposition::Key<3>::value_type>& indicesSendEndGlobal,
                         std::vector<typename decomposition::Key<3>::value_type>& indicesReceiveStartGlobal,
                         std::vector<typename decomposition::Key<3>::value_type>& indicesReceiveEndGlobal,
                         std::vector<GeometricalDomain*>& solidDomains);

  private:
    template <typename T> size_t getConnectivity_private(std::vector<T>& dataset, std::vector<Cell*>* cellsLvl) const;
    template <typename T> void getOffsets_private(std::vector<T>& dataset, std::vector<Cell*>* cellsLvl, uint64_t maxVal) const;
    template <typename T> void getTypeCell_private(std::vector<T>& dataset, std::vector<Cell*>* cellsLvl) const;
    template <typename T> void getData_private(TypeMeshContainer<Cell*>* cellsLvl, std::vector<T>& dataset, const int var, int phase) const;

    int m_lvlMax;                              //!<Niveau maximal sur l arbre AMR (si m_lvlMax = 0, pas d AMR)
    double m_criteriaVar;                      //!<Value of criteria to not pass on the variation of a variable for coarsening or refining (put xi=1.)
    bool m_varRho, m_varP, m_varU, m_varAlpha; //!<Choice on which variation we coarsen or refine
    double m_xiSplit, m_xiJoin;                //!<Value of xi to split or join the cells
    decomposition::Decomposition m_decomp;     //!<Parallel domain decomposition based on keys
};

#endif // MESHCARTESIANAMR_H
