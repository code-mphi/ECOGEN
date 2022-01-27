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

#include "OutputGlobalGNU.h"
#include "../Config.h"

//***************************************************************

OutputGlobalGNU::OutputGlobalGNU() : m_quantity(0.)
{}

//***************************************************************

OutputGlobalGNU::OutputGlobalGNU(std::string casTest, std::string run, tinyxml2::XMLElement* element, Input *entree, std::string nameQuantity) :
  OutputGNU(element)
{
  try {
    //Attributes settings
    m_writeBinary = false;
    m_simulationName = casTest;
    m_fileNameResults = nameQuantity;
    m_fileNameVisu = "plot_" + m_fileNameResults + ".gnu";
    m_folderOutput = config.getWorkFolder() + "results/" + run + "/globalQuantities/";
    m_folderScriptGnuplot = "";
    m_splitData = 0;
    m_numFichier = 0;
    m_input = entree;
    m_run = m_input->getRun();
    m_quantity = 0.;
  }
  catch (ErrorECOGEN &) { throw; }
}

//***************************************************************

OutputGlobalGNU::~OutputGlobalGNU(){}

//***************************************************************

void OutputGlobalGNU::initializeSpecificOutput()
{
  try {
    // Creating output file
    std::ofstream fileStream;
    std::string file = m_folderOutput + createFilenameGNU(m_fileNameResults.c_str());
    if (m_run->m_restartSimulation > 0) {
      fileStream.open(file.c_str(), std::ios_base::app);
    }
    else {
      fileStream.open(file.c_str());
    }
    if (!fileStream) { throw ErrorECOGEN("Cannot open the file " + file, __FILE__, __LINE__); }
    fileStream.close();

    // Gnuplot script printing for visualization
    ecritScriptGnuplot(m_fileNameResults);
  }
  catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

void OutputGlobalGNU::writeResults(Mesh* /*mesh*/, std::vector<Cell*>* cellsLvl)
{
  try {
    this->extractTotalQuantity(cellsLvl);
    if (rankCpu == 0) {
      std::ofstream fileStream;
      std::string file = m_folderOutput + createFilenameGNU(m_fileNameResults.c_str());
      fileStream.open(file.c_str(), std::ios_base::app);
      if (m_precision != 0) fileStream.precision(m_precision);
      fileStream << m_run->m_physicalTime << " " << m_quantity << std::endl;
      fileStream.close();
    }
  }
  catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

void OutputGlobalGNU::extractTotalQuantity(std::vector<Cell*>* cellsLvl)
{
  m_quantity = 0.;
  if (m_fileNameResults == "mass") {
    for (unsigned int c = 0; c < cellsLvl[0].size(); c++) { cellsLvl[0][c]->computeTotalMass(m_quantity); }
  }
  else if (m_fileNameResults == "totalenergy") {
    for (unsigned int c = 0; c < cellsLvl[0].size(); c++) { cellsLvl[0][c]->computeTotalEnergy(m_quantity); }
  }
  else { m_quantity = Errors::defaultDouble; }
  if (Ncpu > 1) { parallel.computeSum(m_quantity); }
}

//***************************************************************