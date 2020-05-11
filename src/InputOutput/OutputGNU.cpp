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

//! \file      OutputGNU.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      May 03 2018

#include "OutputGNU.h"
#include "../Run.h"

using namespace tinyxml2;

//***********************************************************************

OutputGNU::OutputGNU(){}

//***********************************************************************

OutputGNU::OutputGNU(std::string casTest, std::string run, XMLElement *element, std::string fileName, Input *entree) :
  Output(casTest, run, element, fileName, entree)
{
  m_fileNameVisu = "visualization.gnu";
  m_folderScriptGnuplot = "./datasets/";
}

//***********************************************************************

OutputGNU::~OutputGNU(){}

//***********************************************************************

void OutputGNU::ecritSolution(Mesh *mesh, std::vector<Cell *> *cellsLvl)
{
  try {
    std::ofstream fileStream;
    std::string file = m_folderDatasets + creationNameFichierGNU(m_fileNameResults.c_str(), -1, rankCpu, m_numFichier);
    fileStream.open(file.c_str());
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + file, __FILE__, __LINE__); }
    mesh->ecritSolutionGnuplot(cellsLvl, fileStream);
    fileStream << std::endl;
    fileStream.close();

    //Creation du file gnuplot pour visualization des resultats
    if (rankCpu == 0) ecritScriptGnuplot(mesh->getGeometrie());
  }
  catch (ErrorECOGEN &) { throw; }
  m_numFichier++;
}

//*********************************************************************** 

void OutputGNU::ecritScriptGnuplot(const int &dim)
{
  try {
    std::ofstream fileStream;

    fileStream.open((m_folderOutput + m_fileNameVisu).c_str());
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_folderOutput + m_fileNameVisu, __FILE__, __LINE__); }
    fileStream << "reset" << std::endl << "set style data lines" << std::endl;
    fileStream << "set nokey" << std::endl;
    fileStream << std::endl;

    int index;
    
    if (dim == 0) { //Special case probe: first columns : t, etc.
      index = 2;
      fileStream << "set xlabel 't (s)'" << std::endl;
    }
    else{
      index = 1 + dim; //premiere colonne file resultat == premiere donnee apres X,(Y,Z)
      fileStream << "set xlabel 'x (m)'" << std::endl;
    }
    //Gestion 2D/3D
    if (dim == 2) {
      fileStream << "set surface" << std::endl;
      fileStream << "set dgrid3d 50,50" << std::endl;
      fileStream << "set contour base" << std::endl;
      fileStream << "set cntrparam levels 25" << std::endl;
      fileStream << "show contour" << std::endl;
      //fileStream << "set view 180,180,1,1" << endl;
    }
    else if (dim == 3) { throw ErrorECOGEN("OutputGNU::ecritScriptGnuplot : print script gnuplot non prevu en 3D", __FILE__, __LINE__); }

    //1) Variables des phases
    //-----------------------
    for (int phase = 0; phase < m_run->getNumberPhases(); phase++)
    {
      //Variables scalars
      for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberScalars(); var++) {
        fileStream << "set title '" << m_cellRef.getPhase(phase)->returnNameScalar(var) << "_" << m_cellRef.getPhase(phase)->getEos()->getName() << "'" << std::endl;
        printBlocGnuplot(fileStream, index, dim);
      } //Fin var scalar
      //Variables vectorielles (u)
      for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++) {
        fileStream << "set title '" << m_cellRef.getPhase(phase)->returnNameVector(var) << "_" << m_cellRef.getPhase(phase)->getEos()->getName() << "'" << std::endl;
        printBlocGnuplot(fileStream, index, dim);
      } //Fin var vectorielle
    } //Fin phase

   //2) Variables mixture
   //--------------------
   //Variables scalars
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++) {
      fileStream << "set title '" << m_cellRef.getMixture()->returnNameScalar(var) << "'" << std::endl;
      printBlocGnuplot(fileStream, index, dim);
    } //Fin var scalar
    //Variables vectorielle
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++) {
      fileStream << "set title '" << m_cellRef.getMixture()->returnNameVector(var) << "'" << std::endl;
      printBlocGnuplot(fileStream, index, dim);
    } //Fin var vectorielle

    //3) Variables transports
    //-----------------------
    for (int var = 1; var <= m_cellRef.getNumberTransports(); var++) {
      fileStream << "set title 'Transport" << var << "'" << std::endl;
      printBlocGnuplot(fileStream, index, dim);
    } //Fin var scalar

    //4) Ecriture niveaux AMR
    //-----------------------
    fileStream << "set title 'Niveau AMR'" << std::endl;
    printBlocGnuplot(fileStream, index, dim);

    //5) Ecriture variable detection gradients
    //----------------------------------------
    fileStream << "set title 'Xi'" << std::endl;
    printBlocGnuplot(fileStream, index, dim);

    fileStream.close();

  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

//***********************************************************************

void OutputGNU::printBlocGnuplot(std::ofstream &fileStream, int &index, const int &dim)
{
  try {
    if (dim <= 1) { fileStream << "plot"; }
    else { fileStream << "splot"; }
    for (int t = 0; t <= m_numFichier; t++) {
      for (int p = 0; p < Ncpu; p++) {
        fileStream << " \"";
        if(dim==0){ fileStream << m_folderScriptGnuplot + creationNameFichierGNU(m_fileNameResults.c_str(), -1, -1, -1); }
        else { fileStream << m_folderScriptGnuplot + creationNameFichierGNU(m_fileNameResults.c_str(), -1, p, t); }
        fileStream << "\"";
        if (dim <= 1) { fileStream << " u 1:" << index; }
        else { fileStream << " u 1:2:" << index; }
        if (dim == 0) {
          fileStream << std::endl;
          break;
        }
        if (p < Ncpu - 1 || t != m_numFichier) fileStream << ",\\";
        fileStream << std::endl;
      }
    }
    fileStream << "pause(-1)" << std::endl << std::endl;
    index++;
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

std::string OutputGNU::creationNameFichierGNU(const char* name, int lvl, int proc, int numFichier, std::string nameVariable) const
{
  try {
    std::stringstream num;

    if (m_donneesSeparees) {
      throw ErrorECOGEN("OutputGNU::creationNameFichierGNU : donnees Separees non prevu", __FILE__, __LINE__);
      return 0;
    }
    num << name;
    //Gestion nameVariable
    if (nameVariable != "defaut") num << "_" << nameVariable << "_";
    //Gestion binary
    if (m_ecritBinaire)
    {
      //num << "B64";
      throw ErrorECOGEN("OutputGNU::creationNameFichierGNU : donnees binaires non prevu", __FILE__, __LINE__);
      return 0;
    }
    //Gestion cpu
    if (proc > -1) num << "_CPU" << proc;
    //Gestion niveau AMR
    if (lvl != -1)
    {
      //num << "_AMR" << lvl;
      throw ErrorECOGEN("OutputGNU::creationNameFichierGNU : donnees AMR par niveau non prevu", __FILE__, __LINE__);
      return 0;
    }
    //Gestion number de file resultat
    if (numFichier != -1) num << "_TIME" << numFichier;
    //Gestion extension
    num << ".out";
    return num.str();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************