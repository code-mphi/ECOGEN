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

#include "OutputBoundaryAllGNU.h"

//***************************************************************

OutputBoundaryAllGNU::OutputBoundaryAllGNU(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input* entree) : 
  OutputBoundaryGNU(casTest, run, element, fileName, entree)
{
}

//***************************************************************

OutputBoundaryAllGNU::~OutputBoundaryAllGNU() {}

//***************************************************************

void OutputBoundaryAllGNU::initializeSpecificOutputBound()
{
  try {
    // Create output file
    // Avoid to create empty file for CPUs which have not any boundary faces
    if (m_cellInterfaceIndexes.size() > 0) {
      std::ofstream fs;
      std::string file = m_folderOutput + createFilenameGNU(m_fileNameResults.c_str(), -1, rankCpu, m_numFichier);
      if (m_run->m_restartSimulation > 0) {
        fs.open(file.c_str(), std::ios_base::app);
      }
      else {
        fs.open(file.c_str());
      }
      if (!fs) { throw ErrorECOGEN("Cannot open the file " + file, __FILE__, __LINE__); }
      fs.close();
    }
  }
  catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

void OutputBoundaryAllGNU::computeOutwardNormal(CellInterface* bound, double &nx, double &ny, double &nz)
{
  Coord vFaceToElt(0.);
  vFaceToElt.setFromSubtractedVectors(bound->getFace()->getPos(), bound->getCellGauche()->getElement()->getPosition());
  double scalarProduct(Coord::scalarProduct(bound->getFace()->getNormal(), vFaceToElt));

  // Switch normal if not in outward direction
  if (scalarProduct > 0) { 
    nx = - nx;
    ny = - ny;
    nz = - nz;
  }
}
  
//***************************************************************

void OutputBoundaryAllGNU::writeResults(std::vector<CellInterface*>* cellInterfacesLvl)
{
  double rho, p, S, nx, ny, nz;
  // Velocity extracted in the Riemann solver is defined in the face frame
  // and corresponds to the relative velocity in case of MRF.
  // The relative velocity is projected in the global frame to allow the reconstruction
  // of the absolute velocity with the rotationnal velocity if MRF is on.
  Coord absVel;
  CellInterface *bound;

  // Avoid to create empty file for CPUs which have not any boundary faces
  if (m_cellInterfaceIndexes.size() > 0) {
    std::ofstream fs;
    std::string file = m_folderOutput + createFilenameGNU(m_fileNameResults.c_str(), -1, rankCpu, m_numFichier);
    fs.open(file.c_str(), std::ios_base::app);
    if (m_precision != 0) fs.precision(m_precision);

    for (unsigned int c = 0; c < m_cellInterfaceIndexes.size(); c++) {
      bound = cellInterfacesLvl[0][m_cellInterfaceIndexes[c]];
      rho = bound->getBoundData(VarBoundary::rho);

      absVel.setXYZ(bound->getBoundData(VarBoundary::velU), 
                    bound->getBoundData(VarBoundary::velV),
                    bound->getBoundData(VarBoundary::velW));

      if (m_run->m_MRF != -1) {
        absVel.reverseProjection(bound->getFace()->getNormal(), 
          bound->getFace()->getTangent(), 
          bound->getFace()->getBinormal());
        
        absVel = m_run->m_sources[m_run->m_MRF]->computeAbsVelocity(absVel, bound->getCellGauche()->getPosition());
        
        absVel.localProjection(bound->getFace()->getNormal(), bound->getFace()->getTangent(), bound->getFace()->getBinormal());
      }

      p = bound->getBoundData(VarBoundary::p);
      S = bound->getFace()->getSurface();

      // Normal 
      nx = bound->getFace()->getNormal().getX();
      ny = bound->getFace()->getNormal().getY();
      nz = bound->getFace()->getNormal().getZ();
      this->computeOutwardNormal(bound, nx, ny, nz);

      fs << rho << " " 
        << absVel.getX() << " " << absVel.getY() << " " << absVel.getZ() << " " 
        << p << " " << S << " " 
        << nx << " " << ny << " " << nz << std::endl;
    }

    fs.close();
  }

  m_numFichier++;
  m_nextAcq += m_acqFreq;
}

//***************************************************************