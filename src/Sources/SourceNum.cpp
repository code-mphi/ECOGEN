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

#include "SourceNum.h"

enum srcOrder { K1,K2,K3,K4 };

//***********************************************************************

SourceNum::SourceNum(int order, int physicalEntity) : Source(physicalEntity), m_order(order)
{}

//***********************************************************************

SourceNum::~SourceNum()
{}

//***********************************************************************

void SourceNum::integrationEuler(Cell* cell, const double& dt)
{
  sourceCons[K1]->setToZero();                //sourceCons initialized at zero
  sourceCons[K1]->addFlux(cell->getCons());   //sourceCons receives conservative variables at time n: Un
  this->prepSourceTerms(cell);                //Compute the source terms into sourceCons
  sourceCons[K1]->multiply(dt);               //sourceCons is multiplied by dt
}

//***********************************************************************

void SourceNum::integrationRK2(Cell* cell, const double& dt)
{
  // Construct term K1
  this->integrationEuler(cell, dt);

  sourceCons[K2]->setToZero();
  sourceCons[K2]->addFlux(cell->getCons());
  sourceCons[K2]->addFlux(sourceCons[K1]);

  this->prepSourceTerms(cell, K2);
  sourceCons[K2]->multiply(dt);
 
  //RK2 Coefficients
  for (unsigned int m = 0; m < 2; m++){
    sourceCons[m]->multiply(0.5);
  }
}

//***********************************************************************

void SourceNum::integrationRK4(Cell* cell, const double& dt)
{

  // Construct term K1
  this->integrationEuler(cell, dt);

  // Construct term K2
  sourceCons[K2]->setToZero();
  sourceCons[K2]->addFlux(sourceCons[K1]);
  sourceCons[K2]->multiply(0.5);
  sourceCons[K2]->addFlux(cell->getCons());

  this->prepSourceTerms(cell, K2);
  sourceCons[K2]->multiply(dt);

  // Construct term K3
  sourceCons[K3]->setToZero();
  sourceCons[K3]->addFlux(sourceCons[K2]);
  sourceCons[K3]->multiply(0.5);
  sourceCons[K3]->addFlux(cell->getCons());

  this->prepSourceTerms(cell, K3);
  sourceCons[K3]->multiply(dt);

  // Construct term K4
  sourceCons[K4]->setToZero();
  sourceCons[K4]->addFlux(sourceCons[K3]);
  sourceCons[K4]->addFlux(cell->getCons());

  this->prepSourceTerms(cell, K4);
  sourceCons[K4]->multiply(dt);

  // RK4 coefficients
  sourceCons[K1]->multiply(1. / 6.);
  sourceCons[K2]->multiply(1. / 3.);
  sourceCons[K3]->multiply(1. / 3.);
  sourceCons[K4]->multiply(1. / 6.);
}

//***********************************************************************

void SourceNum::integrateSourceTerms(Cell* cell, const double& dt)
{
  if (cell->getElement()->getAppartenancePhysique() == m_physicalEntity || m_physicalEntity == 0) {
    // For unstructured mesh if source term is not applied on specific physicalEntity all cells include it. 
    // For Cartesian mesh, there is no physicalEntity (default value is 0) thus source is applied on all cells.
    cell->buildCons(); // Initialize conservative vector Un
  
    //Deleting old stuff
    for (auto s : sourceCons) { s->setToZero(); }

    // Integration order - Compute Delta t times the sources terms (for Euler)
    if (m_order == 1) { this->integrationEuler(cell, dt); }
    else if (m_order == 2) { this->integrationRK2(cell, dt); }
    else if (m_order == 4) { this->integrationRK4(cell, dt); }

    // Source scheme - Add source terms to Un to obtain Un+1
    for (auto s : sourceCons) {
      cell->getCons()->addFlux(s);    
    }
    cell->buildPrim();
  }
}
