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

#include "OutputBoundaryFluxGNU.h"
#include "../Config.h"

//***************************************************************

OutputBoundaryFluxGNU::OutputBoundaryFluxGNU(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input* entree) : m_flux(0.)
{
  try {
    //Attributes settings
    m_ecritBinaire = false;
    m_simulationName = casTest;
    m_fileNameResults = element->Attribute("name");
    m_fileNameVisu = "visualization_" + m_fileNameResults + ".gnu";
    m_folderOutput = config.getWorkFolder() + "results/" + run + "/boundariesFlux/";
    m_folderScriptGnuplot = "";
    m_donneesSeparees = false;
    m_numFichier = 0;
    m_input = entree;
    m_run = m_input->getRun();

    tinyxml2::XMLElement* sousElement;
    tinyxml2::XMLError error;

    // Reading boundary number
    sousElement = element->FirstChildElement("boundary");
    if (sousElement == NULL) throw ErrorXMLElement("boundary", fileName, __FILE__, __LINE__);
    error = sousElement->QueryIntAttribute("number", &m_numPhys);
    if (error != tinyxml2::XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName, __FILE__, __LINE__);
    error = sousElement->QueryIntAttribute("number", &m_numPhys);

    // Reading flux type (massflow or power flux)
    m_fluxType = sousElement->Attribute("flux");
    if (m_fluxType == "") throw ErrorXMLAttribut("flux", fileName, __FILE__, __LINE__);
    Tools::uppercase(m_fluxType);

    // Reading time control
    sousElement = element->FirstChildElement("timeControl");
    if (sousElement == NULL) throw ErrorXMLElement("timeControl", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("acqFreq", &m_acqFreq);
    if (error != tinyxml2::XML_NO_ERROR) throw ErrorXMLAttribut("acqFreq", fileName, __FILE__, __LINE__);
  }
  catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

OutputBoundaryFluxGNU::~OutputBoundaryFluxGNU()
{
}

//***************************************************************

void OutputBoundaryFluxGNU::prepareSortieSpecifique(std::vector<CellInterface*>* cellInterfacesLvl)
{
  // Set time
  m_nextAcq = 0.;

  // cellInterfaceNumbers keeps indexes of cellInterfaces 
  // which are on the recorded boundary (avoid to loop on all
  // cellInterfaces each time flux is computed)
  //JC//WARNING Currently not working with AMR
  for (unsigned int c = 0; c < cellInterfacesLvl[0].size(); c++) {
    if (cellInterfacesLvl[0][c]->getNumPhys() == m_numPhys) {
      cellInterfaceIndexes.push_back(c);
    }
  }

  try {
    // Create output file
    std::ofstream fileStream;
    std::string file = m_folderOutput + creationNameFichierGNU(m_fileNameResults.c_str());
    fileStream.open(file.c_str());
    if (!fileStream) { throw ErrorECOGEN("Cannot open the file " + file, __FILE__, __LINE__); }
    fileStream.close();

    // Gnuplot script printing for visualization
    ecritScriptGnuplot(m_fileNameResults);
  }
  catch (ErrorECOGEN&) { throw; }

}

//***************************************************************

void OutputBoundaryFluxGNU::ecritSolution(std::vector<CellInterface*>* cellInterfacesLvl)
{
  this->extractFlux(cellInterfacesLvl);

  try {
    if (rankCpu == 0) {
      std::ofstream fileStream;
      std::string file = m_folderOutput + creationNameFichierGNU(m_fileNameResults.c_str());
      fileStream.open(file.c_str(), std::ios_base::app);
      fileStream << m_run->m_physicalTime << " " << m_flux << std::endl;
      fileStream.close();
    }
  }
  catch (ErrorECOGEN&) { throw; }

  m_nextAcq += m_acqFreq;
}

//***************************************************************

void OutputBoundaryFluxGNU::extractFlux(std::vector<CellInterface*>* cellInterfacesLvl)
{
  m_flux = 0.;

  if (m_fluxType == "MASSFLOW") {
    for (unsigned int c = 0; c < cellInterfaceIndexes.size(); c++) {
      m_flux += cellInterfacesLvl[0][cellInterfaceIndexes[c]]->getMassflow();
    }
  }
  else if (m_fluxType == "POWERFLUX") {
    for (unsigned int c = 0; c < cellInterfaceIndexes.size(); c++) {
      m_flux += cellInterfacesLvl[0][cellInterfaceIndexes[c]]->getPowerFlux();
    }
  }
  else {
    m_flux = 0.;
  }
  
  if (Ncpu > 1) {
    parallel.computeSum(m_flux);
  }
}

//***************************************************************