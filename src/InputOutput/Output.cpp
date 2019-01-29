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
//! \version   1.0
//! \date      July 20 2018

#include "Output.h"
#include "../Run.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

Output::Output(){}

//***************************************************************
//Constructeur sortie a partir d une lecture au format XML outputMode
//ex :	<outputMode format="XML" binary="false"/>

Output::Output(string casTest, string nameRun, XMLElement *element, string fileName, Input *entree) :
  m_simulationName(casTest), m_dossierSortie(nameRun), m_numFichier(0), m_donneesSeparees(0), m_input(entree)
{
  //Affectation pointeur run
  m_run = m_input->getRun();

  //Names communs
  //------------
  m_infosCalcul = "infoCalcul.out";
  m_infoMailles = "infosMesh";
  m_treeStructure = "treeStructure";
  m_fileNameResults = "result";
  m_fileNameCollection = "collection";

  m_dossierSortie = "./results/" + m_dossierSortie + "/";
  m_dossierSauvegardesInput = m_dossierSortie + "savesInput/";
  m_dossierSauvegardesInfosMailles = m_dossierSortie + "infosMesh/";
  m_dossierCuts = m_dossierSortie + "cuts/";
  m_folderProbes = m_dossierSortie + "probes/";

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
      _mkdir(m_dossierSortie.c_str());
      _mkdir(m_dossierSauvegardesInput.c_str());
      _mkdir(m_dossierSauvegardesInfosMailles.c_str());
      _mkdir(m_dossierCuts.c_str());
      _mkdir(m_folderProbes.c_str());
    #else
      mkdir("./results", S_IRWXU);
      mkdir(m_dossierSortie.c_str(), S_IRWXU);
      mkdir(m_dossierSauvegardesInput.c_str(), S_IRWXU);
      mkdir(m_dossierSauvegardesInfosMailles.c_str(), S_IRWXU);
      mkdir(m_dossierCuts.c_str(), S_IRWXU);
      mkdir(m_folderProbes.c_str(), S_IRWXU);
    #endif
    try {
      //Sauvegarde des fichiers d entrees
      IO::copieFichier(m_input->getMain(), m_simulationName, m_dossierSauvegardesInput);
      IO::copieFichier(m_input->getMesh(), m_simulationName, m_dossierSauvegardesInput);
      IO::copieFichier(m_input->getCI(), m_simulationName, m_dossierSauvegardesInput);
      IO::copieFichier(m_input->getModel(), m_simulationName, m_dossierSauvegardesInput);
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
    ofstream fileStream;
    //Fichier infosCalcul
    if (rankCpu == 0) fileStream.open((m_dossierSortie + m_infosCalcul).c_str(),ios::trunc); fileStream.close();
    //Fichiers infosMeshes
    string file = m_dossierSauvegardesInfosMailles + creationNameFichier(m_infoMailles.c_str(), -1, rankCpu);
    fileStream.open(file.c_str(), ios::trunc); fileStream.close();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void Output::printTree(Mesh* mesh, vector<Cell *> *cellsLvl)
{
  try {
    ofstream fileStream;
    for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
      string file = m_dossierSauvegardesInfosMailles + creationNameFichier(m_treeStructure.c_str(), lvl, rankCpu, m_numFichier);
      fileStream.open(file.c_str());
      for (unsigned int c = 0; c < cellsLvl[lvl].size(); c++) {
        fileStream << cellsLvl[lvl][c]->getSplit() << " ";
      }
      fileStream.close();
    }
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void Output::readTree(Mesh *mesh, vector<Cell *> *cellsLvl, vector<CellInterface *> *boundariesLvl, const int fileNumber, const vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR)
{
  try {
    ifstream fileStream;
    int splitCell(0);
    string chaine;
    for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
      string file = m_dossierSauvegardesInfosMailles + creationNameFichier(m_treeStructure.c_str(), lvl, rankCpu, fileNumber);
      fileStream.open(file.c_str());
      for (unsigned int c = 0; c < cellsLvl[lvl].size(); c++) {
        fileStream >> splitCell;
        //Raffinement de la cellule
        if (splitCell) mesh->refineCell(cellsLvl[lvl][c],addPhys, model, nbCellsTotalAMR);
      }
      fileStream.close();

      if (Ncpu > 1) Errors::errorMessage("Output::readTree: Resuming with AMR not available");

      //Building cells and interface cells vecors
      //-----------------------------------------
      if (lvl < mesh->getLvlMax()) {
        cellsLvl[lvl+1].clear();
        boundariesLvl[lvl + 1].clear();
        for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->buildLvlCellsAndLvlInternalBoundariesArrays(cellsLvl, boundariesLvl); }
        for (unsigned int i = 0; i < boundariesLvl[lvl].size(); i++) { boundariesLvl[lvl][i]->constructionTableauBordsExternesLvl(boundariesLvl); }
      }
    }
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
  saveInfosMailles();
  cout << "T" << m_run->m_numTest << " | printing file number : " << m_numFichier << "... ";
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
        donneeInt = static_cast<int>(jeuDonnees[k]);
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
  cout << "T" << m_run->m_numTest << " | ------------------------------------------" << endl;
  cout << "T" << m_run->m_numTest << " | RESULTS FILE NUMBER : " << m_numFichier << ",  ITERATION " << m_run->m_iteration << endl;
  cout << "T" << m_run->m_numTest << " |     Physical time       = " << m_run->m_physicalTime << " s " << endl;
  cout << "T" << m_run->m_numTest << " |     Last time step      = " << m_run->m_dt << " s " << endl;
  m_run->m_stat.printScreenStats(m_run->m_numTest);
}

//***********************************************************************

void Output::saveInfos() const
{
  ofstream fileStream;
  if (m_precision!=0) fileStream.precision(m_precision);
  if (rankCpu == 0) {
    fileStream.open((m_dossierSortie + m_infosCalcul).c_str(), ios::app);
    if(m_numFichier==0) fileStream << Ncpu << endl;
    fileStream << m_numFichier << " " << m_run->m_iteration << " " << m_run->m_physicalTime << " " << m_run->m_stat.getComputationTime() << " " << m_run->m_dt;

    //Additional output with purpose to track the radius of a bubble over time and the maximum pressures.
    //To comment if not needed. Be carefull when using it, integration for bubble radius and maximum pressure at the wall are not generalized.
    //-----
    //if (m_run->m_numberPhases > 1) {
    //  double integration(0.);
    //  for (unsigned int c = 0; c < m_run->m_cellsLvl[0].size(); c++) {
    //    m_run->m_cellsLvl[0][c]->computeIntegration(integration);
    //  }
    //  fileStream << " " << integration;
    //}
    ////fileStream << " " << m_run->m_pMax[0] << " " << m_run->m_pMax[1] << " " << m_run->m_pMax[2] << " " << m_run->m_pMax[3];
    ////fileStream << " " << m_run->m_pMaxWall[0] << " " << m_run->m_pMaxWall[1] << " " << m_run->m_pMaxWall[2] << " " << m_run->m_pMaxWall[3];
    //fileStream << " " << m_run->m_pMax[0];
    //fileStream << " " << m_run->m_pMaxWall[0];
    //-----

    fileStream << endl;
    fileStream.close();
  }
}

//***********************************************************************

void Output::readInfos()
{
  fstream fileStream;
  vector<stringstream*> chaine(m_run->m_resumeSimulation + 2); //1 for CPU number, and 1 for initial conditions
  for (unsigned int i = 0; i < chaine.size(); i++) { chaine[i] = new stringstream; }
  string chaineTemp;
  clock_t compTime;
  int numberCPURead;
  int iter(0);
  try {
    fileStream.open((m_dossierSortie + m_infosCalcul).c_str(), ios::in); //Opening in reading mode
                                                                         //Verifying CPU number
    std::getline(fileStream, chaineTemp);
    *(chaine[iter]) << chaineTemp;
    *(chaine[iter]) >> numberCPURead;
    if (numberCPURead != Ncpu) { throw ErrorECOGEN("resume simulation not possible - number of CPU differs from read files"); }
    iter++;
    //Finding corresponding results files
    do {
      std::getline(fileStream, chaineTemp);
      *(chaine[iter]) << chaineTemp;
      *(chaine[iter]) >> m_numFichier >> m_run->m_iteration >> m_run->m_physicalTime >> compTime >> m_run->m_dt;
      iter++;
    } while (m_numFichier != m_run->m_resumeSimulation && !fileStream.eof());
    if (fileStream.eof()) { throw ErrorECOGEN("resume simulation not possible - check file 'infosCalcul.out'"); }
  }
  catch (ErrorECOGEN &) { fileStream.close(); throw; }

  //Erasing end of file
  if (rankCpu == 0) {
    fileStream.close();
    fileStream.open((m_dossierSortie + m_infosCalcul).c_str(), ios::out | ios::trunc); //Opening in printing mode with erasing
    for (unsigned int i = 0; i < chaine.size(); i++) {
      fileStream << chaine[i]->str() << endl;
      delete chaine[i];
    }
    fileStream.close();
  }
  m_run->m_stat.setCompTime(compTime);
}

//***********************************************************************

void Output::saveInfosMailles() const
{
  try {
    ofstream fileStream;
    string file = m_dossierSauvegardesInfosMailles + creationNameFichier(m_infoMailles.c_str(), -1, rankCpu);
    fileStream.open(file.c_str(), ios::app);
    fileStream << m_run->m_nbCellsTotalAMR << endl;
    fileStream.close();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

string Output::creationNameFichier(const char* name, int lvl, int proc, int numFichier) const
{
  try {
    stringstream num;
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