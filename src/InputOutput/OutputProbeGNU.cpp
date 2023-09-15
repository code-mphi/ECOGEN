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

#include "OutputProbeGNU.h"
#include "../Maths/GOVertex.h"
#include "../Config.h"

using namespace tinyxml2;

//***************************************************************

OutputProbeGNU::OutputProbeGNU(std::string casTest, std::string run, XMLElement* element, std::string fileName, Input *entree) : 
  OutputGNU(element)
{
  try {
    //Attributes settings
    m_writeBinary = false;
    m_simulationName = casTest;
    m_fileNameResults = element->Attribute("name");
    m_fileNameVisu = "plot_" + m_fileNameResults + ".gnu";
    m_folderOutput = config.getWorkFolder() + "results/" + run + "/probes/";
	  m_folderScriptGnuplot = "";
    m_splitData = 0;
    m_numFichier = 0;
    m_input = entree;
    m_run = m_input->getRun();
    m_possessesProbe = new bool[Ncpu];

    XMLElement* sousElement;
    XMLError error;

    double donnee;

    //Reading probe data
    Coord vertex;
    sousElement = element->FirstChildElement("vertex");
    if (sousElement == NULL) throw ErrorXMLElement("vertex", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("x", &donnee); vertex.setX(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("y", &donnee); vertex.setY(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("z", &donnee); vertex.setZ(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName, __FILE__, __LINE__);

    m_objet = new GOVertex(vertex);

    sousElement = element->FirstChildElement("timeControl");
    if (sousElement == NULL) throw ErrorXMLElement("timeControl", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("acqFreq", &m_acqFreq);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("acqFreq", fileName, __FILE__, __LINE__);
  }
  catch (ErrorECOGEN &) { throw; }
}

//***************************************************************

OutputProbeGNU::~OutputProbeGNU()
{
  if (m_objet != 0) delete m_objet;
  if(m_possessesProbe !=0) delete m_possessesProbe;
}

//***********************************************************************

void OutputProbeGNU::locateProbeInMesh(const TypeMeshContainer<Cell*>& cells, const int& nbCells, bool localSeeking)
{
  //Locate probe in mesh
  double minimumDistance(1.e12), distance;
  for (int i = 0; i < nbCells; i++) {
    distance = m_objet->distancePoint(cells[i]->getPosition());
    if (distance < minimumDistance && !cells[i]->getWall()) {
      minimumDistance = distance;
      m_cell = cells[i];
      if (distance <= cells[i]->getElement()->getLCFL()) break;
    }
  }
  if(!localSeeking) {
    //Is probe belonging to this CPU ?
    if (Ncpu != 1) {
      double minimumAllCPU(minimumDistance);
      MPI_Allreduce(&minimumDistance, &minimumAllCPU, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      if (std::fabs(minimumAllCPU - minimumDistance) > 1.e-10) { m_possessesProbe[rankCpu] = false; }
      else { m_possessesProbe[rankCpu] = true; }
      MPI_Allgather(&m_possessesProbe[rankCpu], 1, MPI_C_BOOL, m_possessesProbe, 1, MPI_C_BOOL, MPI_COMM_WORLD);
      //Checking for belonging a single CPU
      for(int c=0; c<Ncpu; c++){
        if(m_possessesProbe[c]==true){
          for(int c2=c+1; c2<Ncpu; c2++){
            m_possessesProbe[c2]=false;
          }
        }
      }
    }
    else{
      m_possessesProbe[rankCpu] = true;
    }
  }
}

//***********************************************************************

Cell* OutputProbeGNU::locateProbeInAMRSubMesh(std::vector<Cell*>* cells, const int& nbCells)
{
  int index = 0;

  //Locate probe in AMR Sub-mesh
  double minimumDistance(1.e12), distance;
  for (int i = 0; i < nbCells; i++) {
    distance = m_objet->distancePoint((*cells)[i]->getPosition());
    if (distance < minimumDistance) {
      minimumDistance = distance;
      index = i;
      if (distance <= (*cells)[i]->getElement()->getLCFL()) break;
    }
  }

  if (!(*cells)[index]->getSplit()) { return (*cells)[index]; }
  else {
    return locateProbeInAMRSubMesh((*cells)[index]->getChildVector(), (*cells)[index]->getChildVector()->size());
  }
  return 0;
}

//***********************************************************************

void OutputProbeGNU::initializeSpecificOutput()
{
  //settings
  m_nextAcq = m_run->m_physicalTime;

  //Locate probe in mesh
  locateProbeInMesh(m_run->m_cellsLvl[0], m_run->m_mesh->getNumberCells());

  //Preparing output files
  try {
    if (m_possessesProbe[rankCpu]) {
      //Creating output file
      std::ofstream fileStream;
      std::string file = m_folderOutput + createFilenameGNU(m_fileNameResults.c_str(), -1, -1, -1);
      if (m_run->m_restartSimulation > 0) {
        fileStream.open(file.c_str(), std::ios_base::app);
      }
      else {
        fileStream.open(file.c_str());
      }
      fileStream.close();

      //Gnuplot script printing for visualization
      writeScriptGnuplot(0);
    }
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void OutputProbeGNU::writeResults(Mesh* /*mesh*/, std::vector<Cell*>* /*cellsLvl*/)
{
  std::ofstream fileStream;
  std::string file = m_folderOutput + createFilenameGNU(m_fileNameResults.c_str(), -1, -1, -1);
  fileStream.open(file.c_str(), std::ios_base::app);
  if (m_precision != 0) fileStream.precision(m_precision);
  fileStream << m_run->m_physicalTime << " ";

  //Printing solution with AMR treatement if necessary
  if (!m_cell->getSplit()) {  //if cell is not split
    m_cell->printGnuplotAMR(fileStream, 0, m_objet);
  }
  else { //if cell is split, locate subCells in sub AMR mesh and printing
    locateProbeInAMRSubMesh(m_cell->getChildVector(), m_cell->getChildVector()->size())->printGnuplotAMR(fileStream, 0, m_objet);
  }

  fileStream.close();
  m_nextAcq += m_acqFreq;
}

//***************************************************************
