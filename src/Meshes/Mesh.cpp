//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

//! \file      Mesh.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include "Mesh.h"

//***********************************************************************

Mesh::Mesh() :
  m_numFichier(0)
{
}

//***********************************************************************

Mesh::~Mesh() {}

//***********************************************************************

void Mesh::ecritSolutionGnuplot(std::vector<Cell *> *cellsLvl, std::ofstream &fileStream, GeometricObject *objet) const
{
  for (unsigned int c = 0; c < cellsLvl[0].size(); c++) {
    if (cellsLvl[0][c]->printGnuplotAMR(fileStream, m_geometrie, objet)) break;
  }
}

//****************************************************************************
//****************************** Parallele ***********************************
//****************************************************************************

void Mesh::initializePersistentCommunications(const int numberPhases, const int numberTransports, const TypeMeshContainer<Cell *> &cells, std::string ordreCalcul)
{
	m_numberPhases = numberPhases;
	m_numberTransports = numberTransports;
	int numberVariablesPhaseATransmettre = cells[0]->getPhase(0)->numberOfTransmittedVariables();
	numberVariablesPhaseATransmettre *= m_numberPhases;
	int numberVariablesMixtureATransmettre = cells[0]->getMixture()->numberOfTransmittedVariables();
	int m_numberPrimitiveVariables = numberVariablesPhaseATransmettre + numberVariablesMixtureATransmettre + m_numberTransports;
  int m_numberSlopeVariables(0);
  if (ordreCalcul == "SECONDORDER") {
    int numberSlopesPhaseATransmettre = cells[0]->getPhase(0)->numberOfTransmittedSlopes();
    numberSlopesPhaseATransmettre *= m_numberPhases;
    int numberSlopesMixtureATransmettre = cells[0]->getMixture()->numberOfTransmittedSlopes();
    m_numberSlopeVariables = numberSlopesPhaseATransmettre + numberSlopesMixtureATransmettre + m_numberTransports + 1 + 1; //+1 for the interface detection + 1 for slope index
  }
	parallel.initializePersistentCommunications(m_numberPrimitiveVariables, m_numberSlopeVariables, m_numberTransports, m_geometrie);
}

//***********************************************************************

void Mesh::finalizeParallele(const int &lvlMax)
{
	parallel.finalize(lvlMax);
}

//***********************************************************************
