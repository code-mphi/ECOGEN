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

//***************************************************************

OutputBoundaryFluxGNU::OutputBoundaryFluxGNU(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input* entree) : 
  OutputBoundaryGNU(casTest, run, element, fileName, entree), m_flux(0.)
{
  try {
    // Reading flux type
    tinyxml2::XMLElement* subElement;
    subElement = element->FirstChildElement("flux");
    std::string type(subElement->Attribute("type"));
    if (type == "") { throw ErrorXMLAttribut("type", fileName, __FILE__, __LINE__); }
    Tools::uppercase(type);
    if (type == "MASSFLOW") m_fluxType = FluxType::MASSFLOW;
    else if (type == "POWERFLUX") m_fluxType = FluxType::POWERFLUX;
    else { throw ErrorXMLAttribut("type", fileName, __FILE__, __LINE__); }
  }
  catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

OutputBoundaryFluxGNU::~OutputBoundaryFluxGNU()
{
}

//***************************************************************

void OutputBoundaryFluxGNU::initializeSpecificOutputBound()
{
  try {
    // Create output file
    std::ofstream fs;
    std::string file = m_folderOutput + createFilenameGNU(m_fileNameResults.c_str());
    if (m_run->m_restartSimulation > 0) {
      fs.open(file.c_str(), std::ios_base::app);  
    }
    else {
      fs.open(file.c_str());
    }
    if (!fs) { throw ErrorECOGEN("Cannot open the file " + file, __FILE__, __LINE__); }
    fs.close();

    // Gnuplot script printing for visualization
    ecritScriptGnuplot(m_fileNameResults);
  }
  catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

double OutputBoundaryFluxGNU::getFlux(std::vector<CellInterface*>* cellInterfacesLvl)
{
  if (m_fluxType == FluxType::MASSFLOW) {
    return this->extractMassflow(cellInterfacesLvl);
  }
  else if (m_fluxType == FluxType::POWERFLUX) {
    return this->extractEnthalpyFlux(cellInterfacesLvl);
  }
  else {
    return Errors::defaultDouble;
  }
}

//***************************************************************

double OutputBoundaryFluxGNU::extractMassflow(std::vector<CellInterface*>* cellInterfacesLvl)
{
  m_flux = 0.;

  for (unsigned int c = 0; c < m_cellInterfaceIndexes.size(); c++) {   
    m_flux += this->computeMassflowFace(cellInterfacesLvl[0][m_cellInterfaceIndexes[c]]);
  }
  
  if (Ncpu > 1) {
    parallel.computeSum(m_flux);
  }

  return m_flux;
}

//***************************************************************

double OutputBoundaryFluxGNU::extractEnthalpyFlux(std::vector<CellInterface*>* cellInterfacesLvl)
{
  m_flux = 0.;

  // MRF ON
  if (m_run->m_MRF != -1) {
    for (unsigned int c = 0; c < m_cellInterfaceIndexes.size(); c++) {
      // Absolute velocity is built on the specific region rotating or when whole geometry is rotating.
      if (m_run->m_sources[m_run->m_MRF]->getPhysicalEntity() 
        == cellInterfacesLvl[0][m_cellInterfaceIndexes[c]]->getCellGauche()->getElement()->getAppartenancePhysique()
        || m_run->m_sources[m_run->m_MRF]->getPhysicalEntity() == 0)
      {
        m_flux += this->computeTotalEnthalpyFluxFaceMRF(cellInterfacesLvl[0][m_cellInterfaceIndexes[c]]);
      }
      // If the region is not rotating absolute velocity = relative velocity
      else { 
        m_flux += this->computeTotalEnthalpyFluxFace(cellInterfacesLvl[0][m_cellInterfaceIndexes[c]]); 
      }
    }
  }
  // No MRF
  else {
    for (unsigned int c = 0; c < m_cellInterfaceIndexes.size(); c++) {
      m_flux += this->computeTotalEnthalpyFluxFace(cellInterfacesLvl[0][m_cellInterfaceIndexes[c]]);
    }
  }  

  if (Ncpu > 1) {
    parallel.computeSum(m_flux);
  }

  return m_flux;
}

//***************************************************************

double OutputBoundaryFluxGNU::computeMassflowFace(CellInterface *bound)
{
  double buff(0.);
  // Since velocity is extracted in the Riemann solver, velocityU corresponds to the
  // velocity component in the normal direction
  buff = bound->getBoundData(VarBoundary::rho);
  buff *= bound->getBoundData(VarBoundary::velU);
  buff *= bound->getFace()->getSurface();
  return buff;
}

//***************************************************************

double OutputBoundaryFluxGNU::computeTotalEnthalpyFluxFace(CellInterface *bound)
{
  double buff(0.);
  buff = this->computeMassflowFace(bound);
  buff *= TB->eos[0]->computeTotalEnthalpy(bound->getBoundData(VarBoundary::rho), 
        bound->getBoundData(VarBoundary::p), 
        bound->getBoundData(VarBoundary::velU), 
        bound->getBoundData(VarBoundary::velV), 
        bound->getBoundData(VarBoundary::velW));
  return buff;
}

//***************************************************************

double OutputBoundaryFluxGNU::computeTotalEnthalpyFluxFaceMRF(CellInterface *bound)
{
  double buff(0.);
  buff = this->computeMassflowFace(bound);

  // Velocity extracted in the Riemann solver is defined in the face frame
  // and corresponds to the relative velocity in case of MRF.
  // The relative velocity is projected in the global frame to allow the reconstruction
  // of the absolute velocity with the rotationnal velocity.
  Coord absVel(bound->getBoundData(VarBoundary::velU), 
    bound->getBoundData(VarBoundary::velV), 
    bound->getBoundData(VarBoundary::velW));

  absVel.reverseProjection(bound->getFace()->getNormal(), 
    bound->getFace()->getTangent(), 
    bound->getFace()->getBinormal());

  absVel = m_run->m_sources[m_run->m_MRF]->computeAbsVelocity(absVel, bound->getCellGauche()->getPosition());

  buff *= TB->eos[0]->computeTotalEnthalpy(bound->getBoundData(VarBoundary::rho), 
        bound->getBoundData(VarBoundary::p), 
        absVel);
  return buff;
}

//***************************************************************

void OutputBoundaryFluxGNU::writeResults(std::vector<CellInterface*>* cellInterfacesLvl)
{
  double flux = this->getFlux(cellInterfacesLvl);
  
  try {
    if (rankCpu == 0) {
      std::ofstream fs;
      std::string file = m_folderOutput + createFilenameGNU(m_fileNameResults.c_str());
      fs.open(file.c_str(), std::ios_base::app);
      if (m_precision != 0) fs.precision(m_precision);
      fs << m_run->m_physicalTime << " " << flux << std::endl;
      fs.close();
    }
  }
  catch (ErrorECOGEN&) { throw; }

  m_nextAcq += m_acqFreq;
}

//***************************************************************