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

#ifndef MUSGMSH_H
#define MUSGMSH_H

#include <sstream>
#include "../../MeshUnStruct.h"
#include "HeaderElements.h"

class MUSGmsh : public MeshUnStruct
{
  public:
    MUSGmsh(const std::string& meshFile, const std::string& meshExtension);
    ~MUSGmsh() override;

    static std::string readVersion(const std::string& meshFile);
    void initGeometryMonoCPU_common(TypeMeshContainer<Cell*>& cells);

    // --- MeshUnStruct virtual member functions ---
    void initGeometryMonoCPU(TypeMeshContainer<Cell*>& cells,
                             TypeMeshContainer<CellInterface*>& cellInterfaces,
                             SimulationOrder schemeOrder = SimulationOrder::first) override                   = 0;
    void initGeometryParallel(TypeMeshContainer<Cell*>& cells,
                              TypeMeshContainer<Cell*>& cellsGhost,
                              TypeMeshContainer<CellInterface*>& cellInterfaces,
                              SimulationOrder schemeOrder = SimulationOrder::first) override                  = 0;
    void preProcessMeshFileForParallel() override                                                             = 0;
    void initCpuMeshSequential(TypeMeshContainer<Cell*>& cells, SimulationOrder& schemeOrder) override        = 0;
    void initCpuMeshParallel(TypeMeshContainer<Cell*>& cells, SimulationOrder& schemeOrder, int cpu) override = 0;

    std::string getMeshExtension() const override { return "msh"; };

  private:
};

#endif // MUSGMSH_H
