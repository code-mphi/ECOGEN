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

//! \file      OutputGlobalGNU.h
//! \author    J. Caze
//! \version   1.0
//! \date      January 02 2020

#include "OutputGlobalGNU.h"

//***************************************************************

OutputGlobalGNU::OutputGlobalGNU() : m_quantity(0.)
{}

//***************************************************************

OutputGlobalGNU::OutputGlobalGNU(std::string casTest, std::string run, std::string fileName, Input *entree, std::string nameQuantity)
{
	try {
		//Attributes settings
		m_ecritBinaire = false;
		m_simulationName = casTest;
		m_fileNameResults = nameQuantity;
		m_fileNameVisu = "visualization_" + m_fileNameResults + ".gnu";
		m_folderOutput = "./results/" + run + "/globalQuantities/";
		m_folderScriptGnuplot = m_folderOutput;
		m_donneesSeparees = 0;
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

void OutputGlobalGNU::ecritSolution(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
	try {
		this->extractTotalQuantity(cellsLvl);
		if (rankCpu == 0) {
			std::ofstream fileStream;
			std::string file = m_folderOutput + creationNameFichierGNU(m_fileNameResults.c_str());
			fileStream.open(file.c_str(), std::ios_base::app);
			if (!fileStream) { throw ErrorECOGEN("Cannot open the file " + file, __FILE__, __LINE__); }
			fileStream << m_run->m_physicalTime << " " << m_quantity << std::endl;
			fileStream.close();
		}
	}
	catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

void OutputGlobalGNU::extractTotalQuantity(std::vector<Cell*>* cellsLvl)
{
	if (m_fileNameResults == "mass") {
		m_quantity = 0.;
		for (unsigned int c = 0; c < cellsLvl[0].size(); c++) { cellsLvl[0][c]->computeTotalMass(m_quantity); }
		if (Ncpu > 1) { parallel.computeMassTotal(m_quantity); }
	}
	else if (m_fileNameResults == "energy") { m_quantity = 0.; }
}

//***************************************************************

void OutputGlobalGNU::prepareSortieSpecifique()
{
	try {
		this->writeSpecificGnuplotScript();
	}
	catch (ErrorECOGEN&) { throw; }
}

//***************************************************************

void OutputGlobalGNU::writeSpecificGnuplotScript()
{
	try {
		std::ofstream fileStream;
		fileStream.open((m_folderOutput + m_fileNameVisu).c_str());
		if (!fileStream) { throw ErrorECOGEN("Cannot open the file" + m_folderOutput + m_fileNameVisu, __FILE__, __LINE__); }
		
		fileStream << "reset" << std::endl;
		fileStream << "set style data lines" << std::endl;
		fileStream << "set nokey" << std::endl << std::endl;

		fileStream << "set xlabel 'Time (s)'" << std::endl;
		fileStream << "set title '" << m_fileNameResults << "'" << std::endl;
		
		fileStream << "plot '" << creationNameFichierGNU(m_fileNameResults.c_str(), -1, -1, -1) << "'" << " u 1:2" << std::endl;
		fileStream << "pause(-1)" << std::endl;
	}
	catch (ErrorECOGEN&) { throw; }
}