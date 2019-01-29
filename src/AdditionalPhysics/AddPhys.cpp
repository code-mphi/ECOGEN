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

//! \file      AddPhys.cpp
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "AddPhys.h"

using namespace std;

//***********************************************************************

AddPhys::AddPhys(){}

//***********************************************************************

AddPhys::~AddPhys(){}

//***********************************************************************

void AddPhys::computeFluxAddPhys(CellInterface *cellBound, const int &numberPhases)
{
  this->solveFluxAddPhys(cellBound, numberPhases);

  if (cellBound->getCellGauche()->getLvl() == cellBound->getCellDroite()->getLvl()) {     //CoefAMR = 1 for the two
    this->addFluxAddPhys(cellBound, numberPhases, 1.);                                    //Add flux on the right cell
    this->subtractFluxAddPhys(cellBound, numberPhases, 1.);                               //Subtract flux on the left cell
  }
  else if (cellBound->getCellGauche()->getLvl() > cellBound->getCellDroite()->getLvl()) { //CoefAMR = 1 for the left and 0.5 for the right
    this->addFluxAddPhys(cellBound, numberPhases, 0.5);                                   //Add flux on the right cell
    this->subtractFluxAddPhys(cellBound, numberPhases, 1.);                               //Subtract flux on the left cell
  }
  else {                                                                                  //CoefAMR = 0.5 for the left and 1 for the right
    this->addFluxAddPhys(cellBound, numberPhases, 1.);                                    //Add flux on the right cell
    this->subtractFluxAddPhys(cellBound, numberPhases, 0.5);                              //Subtract flux on the left cell
  }
}

//***********************************************************************

void AddPhys::computeFluxAddPhysBoundary(CellInterface *cellBound, const int &numberPhases)
{
  this->solveFluxAddPhysBoundary(cellBound, numberPhases);
  this->subtractFluxAddPhys(cellBound, numberPhases, 1.); //Subtract flux on the left cell
}

//***********************************************************************

void AddPhys::addFluxAddPhys(CellInterface *cellBound, const int &numberPhases, const double &coefAMR)
{
  double volume(cellBound->getCellDroite()->getElement()->getVolume());
  double surface(cellBound->getFace()->getSurface());
  double coefA(surface / volume); //no "time step"
  coefA = coefA*coefAMR;
  cellBound->getCellDroite()->getCons()->addFlux(coefA, numberPhases);
}

//***********************************************************************

void AddPhys::subtractFluxAddPhys(CellInterface *cellBound, const int &numberPhases, const double &coefAMR)
{
  double volume(cellBound->getCellGauche()->getElement()->getVolume());
  double surface(cellBound->getFace()->getSurface());
  double coefA(surface / volume); //no "time step"
  coefA = coefA*coefAMR;
  cellBound->getCellGauche()->getCons()->subtractFlux(coefA, numberPhases);
}

//***********************************************************************

void AddPhys::addNonConsAddPhys(Cell *cell, const int &numberPhases)
{
  this->addNonCons(cell, numberPhases);
}

//***********************************************************************
