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

#ifndef OUTPUTXML_H
#define OUTPUTXML_H

#include "Output.h"

class OutputXML :  public Output
{
public:
  OutputXML(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input *entree);
  OutputXML(std::string nameRun, int fileNumberRestartMeshMapping, Input *input);
  virtual ~OutputXML();

  virtual void initializeSpecificOutput();
  virtual void writeResults(Mesh *mesh, std::vector<Cell*>* cellsLvl);

  virtual void readResults(Mesh *mesh, std::vector<Cell*>* cellsLvl);
  virtual void readResultsCpu(Mesh *mesh, std::vector<Cell*>* cellsLvl, int cpu);

protected:

  void ReadPhysicalDataXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, tinyxml2::XMLElement* nodeCellData, std::string fileName = "Unknown file");

  std::string createFilenameXML(const char* name, Mesh *mesh=0, int proc=-1, int numFichier=-1, std::string nameVariable ="defaut");

  void writeResultsXML(Mesh *mesh, std::vector<Cell*>* cellsLvl);
  void writeCollectionXML(Mesh *mesh);
  void writePhysicalDataXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, std::ofstream &fileStream, bool parallel = false);

  //Dependant du type de mesh
  void writeMeshRectilinearXML(Mesh *mesh, std::ofstream &fileStream, bool parallel = false);
  void writeMeshUnstructuredXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, std::ofstream &fileStream, bool parallel = false);
  void writeFinFichierRectilinearXML(std::ofstream &fileStream, bool parallel = false);
  void writeFinFichierUnstructuredXML(std::ofstream &fileStream, bool parallel = false);

  //Non used / old
  // void writeFichierParallelXML(Mesh *mesh, std::vector<Cell*>* cellsLvl);
  // void writeFinFichierPolyDataXML(std::ofstream &fileStream, bool parallel = false);
  // void writeMeshPolyDataXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, std::ofstream &fileStream, const int& lvl, bool parallel = false);
};

#endif //OUTPUTXML_H