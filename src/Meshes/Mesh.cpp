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

#include "Mesh.h"

//***********************************************************************

Mesh::Mesh() :
  m_numFichier(0)
{
}

//***********************************************************************

Mesh::~Mesh() {}

//***********************************************************************

void Mesh::writeResultsGnuplot(std::vector<Cell*>* cellsLvl, std::ofstream &fileStream, GeometricObject *objet, bool recordPsat) const
{
  for (unsigned int c = 0; c < cellsLvl[0].size(); c++) {
    if (cellsLvl[0][c]->printGnuplotAMR(fileStream, m_problemDimension, objet, recordPsat)) break;
  }
}

//****************************************************************************
//****************************** Parallele ***********************************
//****************************************************************************

void Mesh::initializePersistentCommunications(const TypeMeshContainer<Cell*>& cells, std::string ordreCalcul)
{
	int numberVariablesPhaseToSend(0);
  for (int k = 0; k < numberPhases; k++) {
    numberVariablesPhaseToSend += cells[0]->getPhase(k)->numberOfTransmittedVariables();
  }
  int numberVariablesMixtureToSend = cells[0]->getMixture()->numberOfTransmittedVariables();
  int m_numberPrimitiveVariables = numberVariablesPhaseToSend + numberVariablesMixtureToSend + numberTransports;
  int m_numberSlopeVariables(0);
  if (ordreCalcul == "SECONDORDER") {
    
    int numberSlopesPhaseToSend(0);
    int numberSlopesMixtureToSend(0);
    int numberSlopesTransportToSend(0);

    if (this->getType() != TypeM::UNS) {
      for (int k = 0; k < numberPhases; k++) {
        numberSlopesPhaseToSend += cells[0]->getPhase(k)->numberOfTransmittedSlopes();
      }
      numberSlopesMixtureToSend = cells[0]->getMixture()->numberOfTransmittedSlopes();
      m_numberSlopeVariables = 1 + 1; //+1 for the interface detection + 1 for slope index
      numberSlopesTransportToSend = numberTransports;
    }
    else {
      for (int k = 0; k < numberPhases; k++) {
        numberSlopesPhaseToSend += cells[0]->getGradPhase(k)->numberOfTransmittedGradients();
      }
      numberSlopesMixtureToSend = cells[0]->getGradMixture()->numberOfTransmittedGradients();
      for (int t = 0; t < numberTransports; t++) {
        numberSlopesTransportToSend += cells[0]->getGradTransport(t)->numberOfTransmittedGradients();
      }
    }
    m_numberSlopeVariables += numberSlopesPhaseToSend + numberSlopesMixtureToSend + numberSlopesTransportToSend;
  }
	parallel.initializePersistentCommunications(m_numberPrimitiveVariables, m_numberSlopeVariables, numberTransports, m_problemDimension);
}

//***********************************************************************

void Mesh::finalizeParallele(const int& lvlMax)
{
	parallel.finalize(lvlMax);
}

//***********************************************************************
