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

//! \file      Source.cpp
//! \author    F. Petitpas, J. Caze
//! \version   1.0
//! \date      October 29 2019

#include "Source.h"

enum srcOrder { K1,K2,K3,K4 };

//***********************************************************************

Source::Source() : m_order(1)
{}

//***********************************************************************

Source::Source(int order) : m_order(order)
{}

//***********************************************************************

Source::~Source()
{}

//***********************************************************************

void Source::integrationEuler(Cell *cell, const int &numberPhases, const double &dt)
{
  sourceCons[K1]->setToZero(numberPhases);
  sourceCons[K1]->addFlux(cell->getCons(), numberPhases);
  this->prepSourceTerms(cell, numberPhases, dt);
  sourceCons[K1]->multiply(dt,numberPhases);
}

//***********************************************************************

void Source::integrationRK2(Cell *cell, const int &numberPhases, const double &dt)
{
  // Construct term K1
  this->integrationEuler(cell, numberPhases, dt);

  sourceCons[K2]->setToZero(numberPhases);
  sourceCons[K2]->addFlux(cell->getCons(),numberPhases);
  sourceCons[K2]->addFlux(sourceCons[K1], numberPhases);

  this->prepSourceTerms(cell, numberPhases, dt, K2);
  sourceCons[K2]->multiply(dt, numberPhases);
 
  //RK2 Coefficients
  for (unsigned int m = 0; m < 2; m++){
    sourceCons[m]->multiply(0.5,numberPhases);
  }
}

//***********************************************************************

void Source::integrationRK4(Cell *cell, const int &numberPhases, const double &dt)
{

  // Construct term K1
  this->integrationEuler(cell, numberPhases, dt);

  // Construct term K2
  sourceCons[K2]->setToZero(numberPhases);
  sourceCons[K2]->addFlux(sourceCons[K1], numberPhases);
  sourceCons[K2]->multiply(0.5, numberPhases);
  sourceCons[K2]->addFlux(cell->getCons(), numberPhases);

  this->prepSourceTerms(cell, numberPhases, dt, K2);
  sourceCons[K2]->multiply(dt, numberPhases);

  // Construct term K3
  sourceCons[K3]->setToZero(numberPhases);
  sourceCons[K3]->addFlux(sourceCons[K2], numberPhases);
  sourceCons[K3]->multiply(0.5, numberPhases);
  sourceCons[K3]->addFlux(cell->getCons(), numberPhases);

  this->prepSourceTerms(cell, numberPhases, dt, K3);
  sourceCons[K3]->multiply(dt, numberPhases);

  // Construct term K4
  sourceCons[K4]->setToZero(numberPhases);
  sourceCons[K4]->addFlux(sourceCons[K3], numberPhases);
  sourceCons[K4]->addFlux(cell->getCons(), numberPhases);

  this->prepSourceTerms(cell, numberPhases, dt, K4);
  sourceCons[K4]->multiply(dt, numberPhases);

  // RK4 coefficients
  sourceCons[K1]->multiply(1. / 6., numberPhases);
  sourceCons[K2]->multiply(1. / 3., numberPhases);
  sourceCons[K3]->multiply(1. / 3., numberPhases);
  sourceCons[K4]->multiply(1. / 6., numberPhases);
}

//***********************************************************************

void Source::integrateSourceTerms(Cell *cell, const int &numberPhases, const double &dt)
{
  cell->buildCons(numberPhases); // Initialize conservative vector U^n
  
  //Deleting old stuff
  for (auto s : sourceCons) { s->setToZero(numberPhases); }

  // Integration order
  if (m_order == 1) { this->integrationEuler(cell, numberPhases, dt); }
  else if (m_order == 2) { this->integrationRK2(cell, numberPhases, dt); }
  else if (m_order == 4) { this->integrationRK4(cell, numberPhases, dt); }
  
  // Source scheme
   for (auto s : sourceCons) {
    cell->getCons()->addFlux(s, numberPhases);    
  }
  cell->buildPrim(numberPhases); 
}
