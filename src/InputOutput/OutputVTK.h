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

#ifndef OUTPUTVTK_H
#define OUTPUTVTK_H

#include "Output.h"

class OutputVTK : public Output
{
  public:
    OutputVTK(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input* entree);
    OutputVTK(std::string nameRun, int fileNumberRestartMeshMapping, Input* input);
    ~OutputVTK() override;

    void initializeSpecificOutput() override;
    void writeResults(Mesh* mesh, std::vector<Cell*>* cellsLvl) override;
    void readResults(Mesh* mesh, std::vector<Cell*>* cellsLvl, int cpu = -1) override;

  protected:
    void ReadPhysicalDataVTK(Mesh* mesh, std::vector<Cell*>* cellsLvl, tinyxml2::XMLElement* nodeCellData, std::string fileName = "Unknown file");

    std::string createFilenameVTK(const char* name, Mesh* mesh = 0, int proc = -1, int numFichier = -1, std::string nameVariable = "defaut");

    void writeResultsVTK(Mesh* mesh, std::vector<Cell*>* cellsLvl);
    void writeCollectionVTK(Mesh* mesh);
    void writePhysicalDataVTK(Mesh* mesh, std::vector<Cell*>* cellsLvl, std::ofstream& fileStream);

    //Dependant du type de mesh
    void writeMeshRectilinearVTK(Mesh* mesh, std::ofstream& fileStream);
    void writeMeshUnstructuredVTK(Mesh* mesh, std::vector<Cell*>* cellsLvl, std::ofstream& fileStream);
    void writeRectilinearVTKFileEnd(std::ofstream& fileStream);
    void writeUnstructuredVTKFileEnd(std::ofstream& fileStream);
};

#endif //OUTPUTVTK_H
