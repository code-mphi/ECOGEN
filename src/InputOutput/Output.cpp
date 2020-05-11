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

//! \file      Output.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "Output.h"
#include "../Run.h"

using namespace tinyxml2;

//***********************************************************************

Output::Output(){}

//***************************************************************
//Constructeur sortie a partir d une lecture au format XML outputMode
//ex :	<outputMode format="XML" binary="false"/>

Output::Output(std::string casTest, std::string nameRun, XMLElement *element, std::string fileName, Input *entree) :
  m_simulationName(casTest), m_folderOutput(nameRun), m_numFichier(0), m_donneesSeparees(0), m_input(entree)
{
  //Affectation pointeur run
  m_run = m_input->getRun();

  //Names communs
  //------------
  m_infoCalcul = "infoCalcul.out";
  m_infoMesh = "infoMesh";
  m_treeStructure = "treeStructure";
  m_domainDecomposition = "domainDecomposition";
  m_fileNameResults = "result";
  m_fileNameCollectionParaview = "collectionParaview";
  m_fileNameCollectionVisIt = "collectionVisIt";

  m_folderOutput = "./results/" + m_folderOutput + "/";
  m_folderSavesInput = m_folderOutput + "savesInput/";
  m_folderDatasets = m_folderOutput + "datasets/";
  m_folderInfoMesh = m_folderOutput + "infoMesh/";
  m_folderCuts = m_folderOutput + "cuts/";
  m_folderProbes = m_folderOutput + "probes/";
  m_folderRecorderGlobalQuantity = m_folderOutput + "globalQuantities/";

  //XMLElement *elementCut;
  XMLError error;

  //Printing precision (digits number)
  if (element->QueryIntAttribute("precision", &m_precision) != XML_NO_ERROR) m_precision = 0; //default if not specified

  //Recuperation mode Ecriture
  error = element->QueryBoolAttribute("binary", &m_ecritBinaire);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("binary", fileName, __FILE__, __LINE__);

  //Creation du dossier de sortie ou vidange /Macro selon OS Windows ou Linux
  if (rankCpu == 0) {
    //Macro pour les interaction systeme (creation/destruction repertoires)
    #ifdef WIN32
      _mkdir("./results");
      _mkdir(m_folderOutput.c_str());
      _mkdir(m_folderSavesInput.c_str());
      _mkdir(m_folderDatasets.c_str());
      _mkdir(m_folderInfoMesh.c_str());
      _mkdir(m_folderCuts.c_str());
      _mkdir(m_folderProbes.c_str());
	  _mkdir(m_folderRecorderGlobalQuantity.c_str());
    #else
      mkdir("./results", S_IRWXU);
      mkdir(m_folderOutput.c_str(), S_IRWXU);
      mkdir(m_folderSavesInput.c_str(), S_IRWXU);
      mkdir(m_folderDatasets.c_str(), S_IRWXU);
      mkdir(m_folderInfoMesh.c_str(), S_IRWXU);
      mkdir(m_folderCuts.c_str(), S_IRWXU);
      mkdir(m_folderProbes.c_str(), S_IRWXU);
	  mkdir(m_folderRecorderGlobalQuantity.c_str(), S_IRWXU);
    #endif
    try {
      //Sauvegarde des fichiers d entrees
      IO::copieFichier(m_input->getMain(), m_simulationName, m_folderSavesInput);
      IO::copieFichier(m_input->getMesh(), m_simulationName, m_folderSavesInput);
      IO::copieFichier(m_input->getCI(), m_simulationName, m_folderSavesInput);
      IO::copieFichier(m_input->getModel(), m_simulationName, m_folderSavesInput);
    }
    catch (ErrorECOGEN &) { throw; }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //Determination du mode Little / Big Endian
  //-----------------------------------------
  int entierTest = 42; //En binary 0x2a
  char *chaineTest = reinterpret_cast<char*>(&entierTest);
  m_endianMode = "LittleEndian";
  if (chaineTest[0] != 0x2a) { m_endianMode = "BigEndian"; }
}

//***********************************************************************

Output::~Output(){}

//***********************************************************************

void Output::prepareOutput(const Cell &cell)
{
  //Preparation de la cell de reference
  //--------------------------------------
  m_cellRef.allocate(m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_addPhys, m_run->m_model);
  for (int k = 0; k < m_run->m_numberPhases; k++) { m_cellRef.copyPhase(k, cell.getPhase(k)); }
  m_cellRef.copyMixture(cell.getMixture());
  for (int k = 0; k < m_run->m_numberTransports; k++) { m_cellRef.setTransport(cell.getTransport(k).getValue(), k); }

  //Preparation propres au type de sortie
  //-------------------------------------
  try {
    this->prepareSortieSpecifique();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void Output::prepareOutputInfos()
{
  try {
    std::ofstream fileStream;
    //Fichier infosCalcul
    if (rankCpu == 0) fileStream.open((m_folderOutput + m_infoCalcul).c_str(),std::ios::trunc); fileStream.close();
    //Fichiers infosMeshes
    std::string file = m_folderInfoMesh + creationNameFichier(m_infoMesh.c_str(), -1, rankCpu);
    fileStream.open(file.c_str(), std::ios::trunc); fileStream.close();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void Output::printTree(Mesh* mesh, std::vector<Cell *> *cellsLvl, int restartAMRsaveFreq)
{
  if (restartAMRsaveFreq != 0) {
    if ((m_numFichier % restartAMRsaveFreq) == 0) {
      try {
        std::ofstream fileStream;
        std::string file;
        //Print domain decomposition
        if (rankCpu == 0) {
          file = m_folderInfoMesh + creationNameFichier(m_domainDecomposition.c_str(), -1, -1, m_numFichier);
          fileStream.open(file.c_str());
          mesh->printDomainDecomposition(fileStream);
          fileStream.close();
        }
        //Print cell tree
        file = m_folderInfoMesh + creationNameFichier(m_treeStructure.c_str(), -1, rankCpu, m_numFichier);
        fileStream.open(file.c_str());
        for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
          for (unsigned int c = 0; c < cellsLvl[lvl].size(); c++) {
            fileStream << cellsLvl[lvl][c]->getSplit() << " ";
          }
        }
        fileStream.close();
      }
      catch (ErrorECOGEN &) { throw; }
    }
  }
}

//***********************************************************************

void Output::readDomainDecompostion(Mesh* mesh, int restartSimulation)
{
  try {
    std::ifstream fileStream;
    std::string file;
    file = m_folderInfoMesh + creationNameFichier(m_domainDecomposition.c_str(), -1, -1, m_numFichier);
    fileStream.open(file.c_str());
    mesh->readDomainDecomposition(fileStream);
    fileStream.close();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void Output::readTree(Mesh *mesh, TypeMeshContainer<Cell *> *cellsLvl, TypeMeshContainer<Cell *> *cellsLvlGhost, TypeMeshContainer<CellInterface *> *cellInterfacesLvl,
  const std::vector<AddPhys*> &addPhys, Model *model, Eos **eos, int &nbCellsTotalAMR)
{
  try {
    std::ifstream fileStream;
    int splitCell(0);
    std::string chaine;
    std::string file = m_folderInfoMesh + creationNameFichier(m_treeStructure.c_str(), -1, rankCpu, m_numFichier);
    fileStream.open(file.c_str());

    for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
      //Refine cells and cell interfaces
      for (unsigned int c = 0; c < cellsLvl[lvl].size(); c++) {
        fileStream >> splitCell;
        if (splitCell) mesh->refineCellAndCellInterfaces(cellsLvl[lvl][c], addPhys, model, nbCellsTotalAMR);
      }

      if (lvl < mesh->getLvlMax()) {
        if (Ncpu > 1) {
          //Refine ghost cells
          parallel.communicationsSplit(lvl);
          cellsLvlGhost[lvl + 1].clear();
          for (unsigned int i = 0; i < cellsLvlGhost[lvl].size(); i++) { cellsLvlGhost[lvl][i]->chooseRefineDeraffineGhost(mesh->getNumberCellsY(), mesh->getNumberCellsZ(), addPhys, model, cellsLvlGhost); }

          //Update of persistent communications of cells lvl + 1
          parallel.communicationsNumberGhostCells(lvl + 1);
          parallel.updatePersistentCommunicationsLvlAMR(lvl + 1, mesh->getGeometrie());
        }

        //Reconstruction of the arrays of cells and cell interfaces of lvl + 1
        cellsLvl[lvl + 1].clear();
        cellInterfacesLvl[lvl + 1].clear();
        for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->buildLvlCellsAndLvlInternalCellInterfacesArrays(cellsLvl, cellInterfacesLvl); }
        for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->constructionTableauCellInterfacesExternesLvl(cellInterfacesLvl); }
      }
    }
    nbCellsTotalAMR = 0;
    for (int i = 0; i < cellsLvl[0].size(); i++) { cellsLvl[0][i]->updateNbCellsTotalAMR(nbCellsTotalAMR); }
    fileStream.close();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void Output::ecritInfos()
{
  if (m_run->m_iteration > 0) {
    afficheInfoEcriture();
  }
  saveInfos();
  std::cout << "T" << m_run->m_numTest << " | Printing file number: " << m_numFichier << "... ";
}

//***********************************************************************

void Output::saveInfosMailles() const
{
  try {
    std::ofstream fileStream;
    std::string file = m_folderInfoMesh + creationNameFichier(m_infoMesh.c_str(), -1, rankCpu);
    fileStream.open(file.c_str(), std::ios::app);
    fileStream << m_run->m_nbCellsTotalAMR << std::endl;
    fileStream.close();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void Output::ecritJeuDonnees(std::vector<double> jeuDonnees, std::ofstream &fileStream, TypeData typeData)
{
  if (m_precision != 0) fileStream.precision(m_precision);
  if (!m_ecritBinaire) {
    for (unsigned int k = 0; k < jeuDonnees.size(); k++) { fileStream << jeuDonnees[k] << " "; }
  }
  else {
    int donneeInt; float donneeFloat; double donneeDouble; char donneeChar;
    int taille;
    switch (typeData) {
    case DOUBLE:
      taille = jeuDonnees.size()*sizeof(double); break;
    case FLOAT:
      taille = jeuDonnees.size()*sizeof(float); break;
    case INT:
      taille = jeuDonnees.size()*sizeof(int); break;
    case CHAR:
      taille = jeuDonnees.size()*sizeof(char); break;
    }
    IO::writeb64(fileStream, taille);
    char *chaineTampon = new char[taille]; int index = 0;
    switch (typeData) {
    case DOUBLE:
      for (unsigned int k = 0; k < jeuDonnees.size(); k++) {
        donneeDouble = static_cast<double>(jeuDonnees[k]);
        IO::ajouteAlaChaine(chaineTampon, index, donneeDouble);
      }
      break;
    case FLOAT:
      for (unsigned int k = 0; k < jeuDonnees.size(); k++) {
        donneeFloat = static_cast<float>(jeuDonnees[k]);
        IO::ajouteAlaChaine(chaineTampon, index, donneeFloat);
      }
      break;
    case INT:
      for (unsigned int k = 0; k < jeuDonnees.size(); k++) {
        donneeInt = static_cast<int>(std::round(jeuDonnees[k]));
        IO::ajouteAlaChaine(chaineTampon, index, donneeInt);
      }
      break;
    case CHAR:
      for (unsigned int k = 0; k < jeuDonnees.size(); k++) {
        donneeChar = static_cast<char>(jeuDonnees[k]);
        IO::ajouteAlaChaine(chaineTampon, index, donneeChar);
      }
      break;
    }
    IO::writeb64Chaine(fileStream, chaineTampon, taille);
    delete[]chaineTampon;
  }
}

//***********************************************************************

void Output::getJeuDonnees(std::istringstream &data, std::vector<double> &jeuDonnees, TypeData typeData)
{
  if (!m_ecritBinaire) {
    for (unsigned int k = 0; k < jeuDonnees.size(); k++) { data >> jeuDonnees[k]; }
  }
  else {
    Errors::errorMessage("resuming on binary results file not available");
    //int donneeInt; float donneeFloat; double donneeDouble; char donneeChar;
    //int taille;
    //switch (typeData) {
    //case DOUBLE:
    //  taille = jeuDonnees.size() * sizeof(double); break;
    //case FLOAT:
    //  taille = jeuDonnees.size() * sizeof(float); break;
    //case INT:
    //  taille = jeuDonnees.size() * sizeof(int); break;
    //case CHAR:
    //  taille = jeuDonnees.size() * sizeof(char); break;
    //}
    //IO::writeb64(fileStream, taille);
    //char *chaineTampon = new char[taille]; int index = 0;
    //switch (typeData) {
    //case DOUBLE:
    //  for (unsigned int k = 0; k < jeuDonnees.size(); k++) {
    //    donneeDouble = static_cast<double>(jeuDonnees[k]);
    //    IO::ajouteAlaChaine(chaineTampon, index, donneeDouble);
    //  }
    //  break;
    //case FLOAT:
    //  for (unsigned int k = 0; k < jeuDonnees.size(); k++) {
    //    donneeFloat = static_cast<float>(jeuDonnees[k]);
    //    IO::ajouteAlaChaine(chaineTampon, index, donneeFloat);
    //  }
    //  break;
    //case INT:
    //  for (unsigned int k = 0; k < jeuDonnees.size(); k++) {
    //    donneeInt = static_cast<int>(jeuDonnees[k]);
    //    IO::ajouteAlaChaine(chaineTampon, index, donneeInt);
    //  }
    //  break;
    //case CHAR:
    //  for (unsigned int k = 0; k < jeuDonnees.size(); k++) {
    //    donneeChar = static_cast<char>(jeuDonnees[k]);
    //    IO::ajouteAlaChaine(chaineTampon, index, donneeChar);
    //  }
    //  break;
    //}
    //IO::writeb64Chaine(fileStream, chaineTampon, taille);
    //delete[]chaineTampon;
  }
}

//***********************************************************************

void Output::afficheInfoEcriture() const
{
  std::cout << "T" << m_run->m_numTest << " | -------------------------------------------" << std::endl;
  std::cout << "T" << m_run->m_numTest << " | RESULTS FILE NUMBER: " << m_numFichier << ", ITERATION " << m_run->m_iteration << std::endl;
  std::cout << "T" << m_run->m_numTest << " |     Physical time       = " << m_run->m_physicalTime << " s " << std::endl;
  std::cout << "T" << m_run->m_numTest << " |     Last time step      = " << m_run->m_dt << " s " << std::endl;
  m_run->m_stat.printScreenStats(m_run->m_numTest);
}

//***********************************************************************

void Output::saveInfos() const
{
  std::ofstream fileStream;
  if (m_precision != 0) fileStream.precision(m_precision);
  if (rankCpu == 0) {
    fileStream.open((m_folderOutput + m_infoCalcul).c_str(), std::ios::app);
    if(m_numFichier == 0) fileStream << Ncpu << std::endl;
    fileStream << m_numFichier << " " << m_run->m_iteration << " " << m_run->m_physicalTime << " " << m_run->m_dtNext
       << " " << m_run->m_stat.getComputationTime() << " " << m_run->m_stat.getAMRTime() << " " << m_run->m_stat.getCommunicationTime();

    //Additional output with purpose to track the radius of a bubble over time and the maximum pressures.
    //To comment if not needed. Be carefull when using it, integration for bubble radius and maximum pressure at the wall are not generalized.
    //-----
    // if (m_run->m_numberPhases > 1) {
    //  double integration(0.);
    //  for (unsigned int c = 0; c < m_run->m_cellsLvl[0].size(); c++) {
    //    m_run->m_cellsLvl[0][c]->computeIntegration(integration);
    //  }
    //  fileStream << " " << integration;
    // }
    // // fileStream << " " << m_run->m_pMax[0] << " " << m_run->m_pMax[1] << " " << m_run->m_pMax[2] << " " << m_run->m_pMax[3];
    // // fileStream << " " << m_run->m_pMaxWall[0] << " " << m_run->m_pMaxWall[1] << " " << m_run->m_pMaxWall[2] << " " << m_run->m_pMaxWall[3];
    // // fileStream << " " << m_run->m_pMax[0];
    // // fileStream << " " << m_run->m_pMaxWall[0];
    // fileStream << " " << m_run->m_alphaWanted;
    //-----

    fileStream << std::endl;
    fileStream.close();
  }
}

//***********************************************************************

void Output::readInfos()
{
  std::fstream fileStream;
  std::vector<std::stringstream*> chaine(m_run->m_restartSimulation + 2); //1 for CPU number, and 1 for initial conditions
  for (unsigned int i = 0; i < chaine.size(); i++) { chaine[i] = new std::stringstream; }
  std::string chaineTemp;
  clock_t compTime;
  clock_t AMRTime;
  clock_t comTime;
  int numberCPURead;
  int iter(0);
  try {
    fileStream.open((m_folderOutput + m_infoCalcul).c_str(), std::ios::in); //Opening in reading mode
    //Verifying CPU number
    std::getline(fileStream, chaineTemp);
    *(chaine[iter]) << chaineTemp;
    *(chaine[iter]) >> numberCPURead;
    if (numberCPURead != Ncpu) { throw ErrorECOGEN("restart simulation not possible - number of CPU differs from read files"); }
    iter++;
    //Finding corresponding results files
    do {
      std::getline(fileStream, chaineTemp);
      *(chaine[iter]) << chaineTemp;
      *(chaine[iter]) >> m_numFichier >> m_run->m_iteration >> m_run->m_physicalTime >> m_run->m_dt >> compTime >> AMRTime >> comTime;
      iter++;
    } while (m_numFichier != m_run->m_restartSimulation && !fileStream.eof());
    if (fileStream.eof()) { throw ErrorECOGEN("restart simulation not possible - check file 'infosCalcul.out'"); }
  }
  catch (ErrorECOGEN &) { fileStream.close(); throw; }

  //Erasing end of file
  MPI_Barrier(MPI_COMM_WORLD);
  if (rankCpu == 0) {
    fileStream.close();
    fileStream.open((m_folderOutput + m_infoCalcul).c_str(), std::ios::out | std::ios::trunc); //Opening in printing mode with erasing
    for (unsigned int i = 0; i < chaine.size(); i++) {
      fileStream << chaine[i]->str() << std::endl;
      delete chaine[i];
    }
    fileStream.close();
  }
  m_run->m_stat.setCompTime(compTime, AMRTime, comTime);
}

//***********************************************************************

std::string Output::creationNameFichier(const char* name, int lvl, int proc, int numFichier) const
{
  try {
    std::stringstream num;
    num << name;
    //Gestion cpu
    if (proc > -1) num << "_CPU" << proc;
    //Gestion niveau AMR
    if (lvl != -1) num << "_AMR" << lvl;
    //Gestion number de file resultat
    if (numFichier != -1) num << "_TIME" << numFichier;
    //Gestion extension
    num << ".out";
    return num.str();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************