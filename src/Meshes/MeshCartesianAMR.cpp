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

//! \file      MeshCartesianAMR.cpp
//! \author    K. Schmidmayer, F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include <iostream>
#include <algorithm>
#include <sstream>

#include "MeshCartesianAMR.h"

using namespace std;

//***********************************************************************

MeshCartesianAMR::MeshCartesianAMR(double lX, int numberCellsX, double lY, int numberCellsY, double lZ, int numberCellsZ,
  std::vector<stretchZone> stretchX, std::vector<stretchZone> stretchY, std::vector<stretchZone> stretchZ,
	int lvlMax, double criteriaVar, bool varRho, bool varP, bool varU, bool varAlpha, double xiSplit, double xiJoin) :
  MeshCartesian(lX, numberCellsX, lY, numberCellsY, lZ, numberCellsZ, stretchX, stretchY, stretchZ),
  m_lvlMax(lvlMax), m_criteriaVar(criteriaVar), m_varRho(varRho), m_varP(varP), m_varU(varU), m_varAlpha(varAlpha), m_xiSplit(xiSplit), m_xiJoin(xiJoin)
{
  m_type = AMR;
}

//***********************************************************************

MeshCartesianAMR::~MeshCartesianAMR(){
  if (Ncpu > 1) delete[] m_cellsLvlGhost;
}

//***********************************************************************

void MeshCartesianAMR::genereTableauxCellsBordsLvl(Cell **cells, CellInterface **bord, vector<Cell *> **cellsLvl,
	vector<CellInterface *> **boundariesLvl)
{
	(*cellsLvl) = new vector<Cell *>[m_lvlMax + 1];
	for (int i = 0; i < m_numberCellsCalcul; i++) { (*cellsLvl)[0].push_back(cells[i]); }

	(*boundariesLvl) = new vector<CellInterface *>[m_lvlMax + 1];
	for (int i = 0; i < m_numberFacesTotal; i++) { (*boundariesLvl)[0].push_back(bord[i]); }

	m_cellsLvl = cellsLvl;

	if (Ncpu > 1) {
		//Genere les tableaux de cells fantomes par niveau
		m_cellsLvlGhost = new vector<Cell *>[m_lvlMax + 1];
		for (int i = m_numberCellsCalcul; i < m_numberCellsTotal; i++) {
			m_cellsLvlGhost[0].push_back(cells[i]);
		}
	}
}

//***********************************************************************

void MeshCartesianAMR::procedureRaffinementInitialization(vector<Cell *> *cellsLvl, vector<CellInterface *> *boundariesLvl,
  const vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR, vector<GeometricalDomain*> &domains,	Cell **cells, Eos **eos, const int &resumeSimulation)
{
  nbCellsTotalAMR = m_numberCellsCalcul;

  if (resumeSimulation == 0) { //Only for simulation from input files
    for (int iterInit = 0; iterInit < 2; iterInit++) {
      for (int lvl = 0; lvl < m_lvlMax; lvl++) {
        if (Ncpu > 1) { parallel.communicationsPrimitivesAMR(cells, eos, lvl); }
        this->procedureRaffinement(cellsLvl, boundariesLvl, lvl, addPhys, model, nbCellsTotalAMR, cells, eos);
        for (unsigned int i = 0; i < cellsLvl[lvl + 1].size(); i++) {
          cellsLvl[lvl + 1][i]->fill(domains, m_lvlMax);
        }
        for (unsigned int i = 0; i < cellsLvl[lvl + 1].size(); i++) {
          cellsLvl[lvl + 1][i]->completeFulfillState();
        }
        for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
          cellsLvl[lvl][i]->averageChildrenInParent();
        }
      }
    }
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      if (Ncpu > 1) { parallel.communicationsPrimitivesAMR(cells, eos, lvl); }
      for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
        if (!cellsLvl[lvl][i]->getSplit()) { cellsLvl[lvl][i]->completeFulfillState(); }
      }
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::procedureRaffinement(vector<Cell *> *cellsLvl, vector<CellInterface *> *boundariesLvl, const int &lvl,
  const vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR, Cell **cells, Eos **eos)
{
  //1) Calcul de Xi dans chaque cell de niveau lvl
  //-------------------------------------------------
  for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->setToZeroXi(); }
  for (unsigned int i = 0; i < boundariesLvl[lvl].size(); i++) { boundariesLvl[lvl][i]->computeXi(m_criteriaVar, m_varRho, m_varP, m_varU, m_varAlpha); }
  //bool varP2(true);
  //if (lvl >= 5) { varP2 = false; }
  //for (unsigned int i = 0; i < boundariesLvl[lvl].size(); i++) { boundariesLvl[lvl][i]->computeXi(m_criteriaVar, m_varRho, varP2, m_varU, m_varAlpha); }
  //for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
  //  double x(0.), y(0.), z(0.);
  //  x = cellsLvl[lvl][i]->getPosition().getX();
  //  y = cellsLvl[lvl][i]->getPosition().getY();
  //  //z = cellsLvl[lvl][i]->getPosition().getZ();
  //  //if (pow((x*x + y*y + z*z), 0.5) > 500.e-6) {
  //  //if (pow((x*x + y*y), 0.5) > 6.e-4) {
  //  //if ((x > 250e-6) || (y > 200.e-6)) {
  //  //if (x > 15.) {
  //  if (pow((x*x + y * y), 0.5) > 5.) {
  //      cellsLvl[lvl][i]->setToZeroXi();
  //  }
  //}
  if (Ncpu > 1) { parallel.communicationsXi(cells, lvl); }
  
  //2) Smoothing de Xi
  //------------------
  for (int iterDiff = 0; iterDiff < 2; iterDiff++) { //Arbitrary number of iterations
		//Mise a zero cons xi
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->setToZeroConsXi(); }

    //Calcul des "flux"
    for (unsigned int i = 0; i < boundariesLvl[lvl].size(); i++) { boundariesLvl[lvl][i]->computeFluxXi(); }

    //Evolution temporelle
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->timeEvolutionXi(); }
		if (Ncpu > 1) { parallel.communicationsXi(cells, lvl); }
  }

	if (lvl < m_lvlMax) {
    int lvlPlus1 = lvl + 1;
    //3) Raffinement des cells et boundaries
    //------------------------------------
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->chooseRefine(m_xiSplit, m_numberCellsY, m_numberCellsZ, addPhys, model, nbCellsTotalAMR); }

    //4) Deraffinement des cells et boundaries
    //--------------------------------------
    bool deraffine = false;
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->chooseUnrefine(m_xiJoin, nbCellsTotalAMR); }

    if (Ncpu > 1) {
      //5) Raffinement et deraffinement des cells fantomes
      //-----------------------------------------------------
      //Communication split + Raffinement et deraffinement des cells fantomes + Reconstruction du tableau de cells fantomes de niveau lvl + 1
      parallel.communicationsSplit(cells, lvl);
      m_cellsLvlGhost[lvlPlus1].clear();
      for (unsigned int i = 0; i < m_cellsLvlGhost[lvl].size(); i++) { m_cellsLvlGhost[lvl][i]->chooseRefineDeraffineGhost(m_numberCellsY, m_numberCellsZ, addPhys, model, m_cellsLvlGhost); }
      //Communications primitives pour mettre a jour les cells deraffinees
      parallel.communicationsPrimitivesAMR(cells, eos, lvl);

      //6) Mise a jour des communications persistantes au niveau lvl + 1
      //----------------------------------------------------------------
      parallel.communicationsNumberGhostCells(cells, lvlPlus1);	//Communication des numbers d'elements a envoyer et a recevoir de chaque cote de la limite parallele
      parallel.updatePersistentCommunicationsLvl(lvlPlus1, m_geometrie);
    }

    //7) Reconstruction des tableaux de cells et boundaries lvl + 1
    //-----------------------------------------------------------
    cellsLvl[lvlPlus1].clear();
    boundariesLvl[lvlPlus1].clear();
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->buildLvlCellsAndLvlInternalBoundariesArrays(cellsLvl, boundariesLvl); }
    for (unsigned int i = 0; i < boundariesLvl[lvl].size(); i++) { boundariesLvl[lvl][i]->constructionTableauBordsExternesLvl(boundariesLvl); }
  }
}

//***********************************************************************

string MeshCartesianAMR::whoAmI() const
{
  return "CARTESIAN_AMR";
}

//**************************************************************************
//******************************** PRINTING ********************************
//**************************************************************************

void MeshCartesianAMR::ecritHeaderPiece(std::ofstream &fileStream, std::vector<Cell *> *cellsLvl, int lvl) const
{
  int numberCells = 0, numberPointsParMaille = 4;
  for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
    if (!cellsLvl[lvl][i]->getSplit()) { numberCells += 1; }
  }
  if (m_numberCellsZ > 1) { numberPointsParMaille = 8; }

  fileStream << "    <Piece NumberOfPoints=\"" << numberPointsParMaille*numberCells << "\" NumberOfCells=\"" << numberCells << "\">" << endl;
}

//***********************************************************************

void MeshCartesianAMR::recupereNoeuds(std::vector<double> &jeuDonnees, int lvl) const
{
  int dimZ = 0;
  if (m_numberCellsZ > 1) dimZ = 1;

  double dXsur2(0.), dYsur2(0.), dZsur2(0.);
  for (unsigned int i = 0; i < (*m_cellsLvl)[lvl].size(); i++) {
    if (!(*m_cellsLvl)[lvl][i]->getSplit()) {

      dXsur2 = 0.5*(*m_cellsLvl)[lvl][i]->getSizeX();
      dYsur2 = 0.5*(*m_cellsLvl)[lvl][i]->getSizeY();
      dZsur2 = 0.5*(*m_cellsLvl)[lvl][i]->getSizeZ();
      //Point 0
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getX() - dXsur2);
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getY() - dYsur2);
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
      //Point 1
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getX() + dXsur2);
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getY() - dYsur2);
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
      //Point 2
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getX() + dXsur2);
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getY() + dYsur2);
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
      //Point 3
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getX() - dXsur2);
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getY() + dYsur2);
      jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getZ() - dZsur2*dimZ);

      if (dimZ > 0.99) {
        //Point 4
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getX() - dXsur2);
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getY() - dYsur2);
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getZ() + dZsur2);
        //Point 5
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getX() + dXsur2);
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getY() - dYsur2);
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getZ() + dZsur2);
        //Point 6
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getX() + dXsur2);
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getY() + dYsur2);
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getZ() + dZsur2);
        //Point 7
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getX() - dXsur2);
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getY() + dYsur2);
        jeuDonnees.push_back((*m_cellsLvl)[lvl][i]->getPosition().getZ() + dZsur2);
      }
    } //Fin cell non split
  } //Fin Cells
}

//***********************************************************************

void MeshCartesianAMR::recupereConnectivite(std::vector<double> &jeuDonnees, int lvl) const
{
  int dimZ(0);
  if (m_numberCellsZ > 1) dimZ = 1;
  int numberPointsParMaille(4);
  if (m_numberCellsZ > 1) { numberPointsParMaille = 8; }

  if (dimZ < 0.99) {
    int numCell(0);
    for (unsigned int i = 0; i < (*m_cellsLvl)[lvl].size(); i++) {
      if (!(*m_cellsLvl)[lvl][i]->getSplit()) {
        jeuDonnees.push_back(numCell*numberPointsParMaille);
        jeuDonnees.push_back(numCell*numberPointsParMaille+1);
        jeuDonnees.push_back(numCell*numberPointsParMaille+2);
        jeuDonnees.push_back(numCell*numberPointsParMaille+3);
        numCell++;
      }
    }
  }
  else {
    int numCell(0);
    for (unsigned int i = 0; i < (*m_cellsLvl)[lvl].size(); i++) {
      if (!(*m_cellsLvl)[lvl][i]->getSplit()) {
        jeuDonnees.push_back(numCell*numberPointsParMaille);
        jeuDonnees.push_back(numCell*numberPointsParMaille + 1);
        jeuDonnees.push_back(numCell*numberPointsParMaille + 2);
        jeuDonnees.push_back(numCell*numberPointsParMaille + 3);
        jeuDonnees.push_back(numCell*numberPointsParMaille + 4);
        jeuDonnees.push_back(numCell*numberPointsParMaille + 5);
        jeuDonnees.push_back(numCell*numberPointsParMaille + 6);
        jeuDonnees.push_back(numCell*numberPointsParMaille + 7);
        numCell++;
      }
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::recupereOffsets(std::vector<double> &jeuDonnees, int lvl) const
{
  int numberPointsParMaille(4);
  if (m_numberCellsZ > 1) { numberPointsParMaille = 8; }
  int numCell(0);
  for (unsigned int i = 0; i < (*m_cellsLvl)[lvl].size(); i++) {
    if (!(*m_cellsLvl)[lvl][i]->getSplit()) {
      jeuDonnees.push_back((numCell + 1)*numberPointsParMaille);
      numCell++;
    }
  }
}

//****************************************************************************

void MeshCartesianAMR::recupereTypeCell(std::vector<double> &jeuDonnees, int lvl) const
{
  int type(9);
  if (m_numberCellsZ > 1) { type = 12; }
  int numCell(0);
  for (unsigned int i = 0; i < (*m_cellsLvl)[lvl].size(); i++) {
    if (!(*m_cellsLvl)[lvl][i]->getSplit()) {
      jeuDonnees.push_back(type);
      numCell++;
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::recupereDonnees(vector<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl) const
{
  jeuDonnees.clear();
  for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
    if (!cellsLvl[lvl][i]->getSplit()) {
      if (var > 0) { //On veut recuperer les donnees scalars
        if (phase >= 0) { jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnScalar(var)); }      //Donnees de phases
        else if (phase == -1) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnScalar(var)); }   //Donnees de mixture
        else if (phase == -2) { jeuDonnees.push_back(cellsLvl[lvl][i]->getTransport(var - 1).getValue()); }
        else if (phase == -3) { jeuDonnees.push_back(cellsLvl[lvl][i]->getXi()); }
        else if (phase == -4) { jeuDonnees.push_back(cellsLvl[lvl][i]->getGradient()); }
        else { Errors::errorMessage("MeshCartesianAMR::recupereDonnees: unknown number of phase: ", phase); }
      }
      else { //On veut recuperer les donnees vectorielles
        if (phase >= 0) { //Donnees de phases
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getX());
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getY());
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getZ());
        }
        else if(phase == -1){  //Donnees de mixture
          jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getX());
          jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getY());
          jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getZ());
        }
        else { Errors::errorMessage("MeshCartesianAMR::recupereDonnees: unknown number of phase: ", phase); }
      } //Fin vecteur
    } //Fin split
  } //fin lvl
}

//****************************************************************************

void MeshCartesianAMR::setDataSet(std::vector<double> &jeuDonnees, vector<Cell *> *cellsLvl, const int var, int phase, int lvl) const
{
  int iterDataSet(0);
  Coord vec;
  for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
    if (!cellsLvl[lvl][i]->getSplit()) {
      if (var > 0) { //Scalars data are first set
        if (phase >= 0) { cellsLvl[lvl][i]->getPhase(phase)->setScalar(var, jeuDonnees[iterDataSet++]); } //phases data
        else if (phase == -1) { cellsLvl[lvl][i]->getMixture()->setScalar(var, jeuDonnees[iterDataSet++]); }  //mixture data
        else if (phase == -2) { cellsLvl[lvl][i]->getTransport(var - 1).setValue(jeuDonnees[iterDataSet++]); } //transport data
        else if (phase == -3) { cellsLvl[lvl][i]->setXi(jeuDonnees[iterDataSet++]); } //xi indicator
        else { Errors::errorMessage("MeshCartesianAMR::setDataSet: unknown phase number: ", phase); }
      }
      else { //On veut recuperer les donnees vectorielles
        if (phase >= 0) { //Phases data
          vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet + 1], jeuDonnees[iterDataSet + 2]);
          cellsLvl[lvl][i]->getPhase(phase)->setVector(-var, vec);
          iterDataSet += 3;
        }
        else if (phase == -1) {  //Mixture data
          vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet + 1], jeuDonnees[iterDataSet + 2]);
          cellsLvl[lvl][i]->getMixture()->setVector(-var, vec);
          iterDataSet += 3;
        }
        else { Errors::errorMessage("MeshCartesianAMR::setDataSet: unknown phase number: ", phase); }
      } //Fin vecteur
    } // Fin split
  } // Fin lvl
}

//***********************************************************************

void MeshCartesianAMR::refineCell(Cell *cell, const vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR)
{
  cell->refineCellAndBoundaries(m_numberCellsY, m_numberCellsZ, addPhys, model);
  nbCellsTotalAMR += cell->getNumberCellsChildren() - 1;
}

//****************************************************************************
//****************************** Parallele ***********************************
//****************************************************************************

void MeshCartesianAMR::initializePersistentCommunications(const int numberPhases, const int numberTransports, Cell **cells, string ordreCalcul)
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
    m_numberSlopeVariables = numberSlopesPhaseATransmettre + numberSlopesMixtureATransmettre + m_numberTransports + 1; //+1 for the interface detection
  }
	parallel.initializePersistentCommunicationsAMR(m_numberPrimitiveVariables, m_numberSlopeVariables, m_numberTransports, m_geometrie, m_lvlMax);
}

//***********************************************************************

void MeshCartesianAMR::communicationsPrimitives(Cell **cells, Eos **eos, const int &lvl, Prim type)
{
	parallel.communicationsPrimitivesAMR(cells, eos, lvl, type);
}

//***********************************************************************

void MeshCartesianAMR::communicationsSlopes(Cell **cells, const int &lvl)
{
	parallel.communicationsSlopesAMR(cells, lvl);
}

//***********************************************************************

void MeshCartesianAMR::communicationsVector(Cell **cells, string nameVector, const int &dim, const int &lvl, int num, int index)
{
	parallel.communicationsVectorAMR(cells, nameVector, m_geometrie, lvl, num, index);
}

//***********************************************************************

void MeshCartesianAMR::communicationsAddPhys(const vector<AddPhys*> &addPhys, Cell **cells, const int &lvl)
{
	for (unsigned int pa = 0; pa < addPhys.size(); pa++) { addPhys[pa]->communicationsAddPhysAMR(cells, m_geometrie, lvl); }
}

//***********************************************************************

void MeshCartesianAMR::communicationsTransports(Cell **cells, const int &lvl)
{
  parallel.communicationsTransportsAMR(cells, lvl);
}

//***********************************************************************

void MeshCartesianAMR::finalizeParallele(const int &lvlMax)
{
	parallel.finalizeAMR(lvlMax);
}

//***********************************************************************