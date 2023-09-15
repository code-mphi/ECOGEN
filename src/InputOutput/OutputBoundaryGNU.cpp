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

#include "OutputBoundaryGNU.h"
#include "../Config.h"

//***************************************************************

OutputBoundaryGNU::OutputBoundaryGNU(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input* entree) : 
  OutputGNU(element)
{
  try {
    tinyxml2::XMLElement* subElement;
    tinyxml2::XMLError error;

    // Reading boundary number/name
    subElement = element->FirstChildElement("boundaryID");
    if (subElement == NULL) throw ErrorXMLElement("boundaryID", fileName, __FILE__, __LINE__);
    error = subElement->QueryIntAttribute("number", &m_numPhys);
    if (error != tinyxml2::XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName, __FILE__, __LINE__);
    std::string name(subElement->Attribute("name"));
    if (name == "") { throw ErrorXMLAttribut("name", fileName, __FILE__, __LINE__); }

    // Attributes settings
    m_writeBinary = false;
    m_simulationName = casTest;
    m_fileNameResults = name;
    m_fileNameVisu = "plot_" + m_fileNameResults + ".gnu";
    m_folderOutput = config.getWorkFolder() + "results/" + run + "/boundaries/";
    m_folderScriptGnuplot = "";
    m_splitData = false;
    m_numFichier = 0;
    m_input = entree;
    m_run = m_input->getRun();

    // Reading time control
    subElement = element->FirstChildElement("timeControl");
    if (subElement == NULL) throw ErrorXMLElement("timeControl", fileName, __FILE__, __LINE__);
    error = subElement->QueryDoubleAttribute("acqFreq", &m_acqFreq);
    if (error != tinyxml2::XML_NO_ERROR) throw ErrorXMLAttribut("acqFreq", fileName, __FILE__, __LINE__);
  }
  catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

OutputBoundaryGNU::~OutputBoundaryGNU() {}

//***************************************************************

void OutputBoundaryGNU::initializeSpecificOutput(std::vector<CellInterface*>* cellInterfacesLvl)
{
  try {
    // Set time (zero if initial run and restart if restart activated)
    m_nextAcq = m_run->m_physicalTime;

    // cellInterfaceNumbers keeps indexes of cellInterfaces 
    // which are on the recorded boundary (avoid to loop on all
    // cellInterfaces each time flux is computed)
    if (m_run->m_mesh->getType() != AMR) {
      for (unsigned int c = 0; c < cellInterfacesLvl[0].size(); c++) {
        if (cellInterfacesLvl[0][c]->getNumPhys() == m_numPhys) {
          m_cellInterfaceIndexes.push_back(c);
        }
      }
    }
    else {
      throw ErrorECOGEN("Recording of boundary not available with AMR mesh", __FILE__, __LINE__); 
    }

    this->initializeSpecificOutputBound();
  }
  catch (ErrorECOGEN&) { throw; }
}

//***************************************************************