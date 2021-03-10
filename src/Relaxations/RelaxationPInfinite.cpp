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

#include "RelaxationPInfinite.h"

//***********************************************************************

RelaxationPInfinite::RelaxationPInfinite(){}

//***********************************************************************

RelaxationPInfinite::~RelaxationPInfinite(){}

//***********************************************************************

void RelaxationPInfinite::relaxation(Cell* cell, const int& numberPhases, const double& /*dt*/, Prim type)
{
  //Is the pressure-relaxation procedure necessary?
  //If alpha = 0 is activated, a test is done to know if the relaxation procedure is necessary or not
  //Else, i.e. alpha = 0 is desactivated (alpha != 0), the relaxation procedure is always done (relax = true)
  bool relax(true);
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    for (int k = 0; k < numberPhases; k++) {
      if (cell->getPhase(k, type)->getAlpha() > (1. - 1.e-5)) relax = false;
    }
  }

  if (relax) {
    Phase* phase(0);
    //Initial state
    double pStar(0.);
    for (int k = 0; k < numberPhases; k++) {
      phase = cell->getPhase(k, type);
      phase->verifyAndCorrectPhase();
      TB->ak[k] = phase->getAlpha();
      TB->pk[k] = phase->getPressure();
      TB->rhok[k] = phase->getDensity();
      pStar += TB->ak[k] * TB->pk[k];
      //phase->verifyPhase();
    }

    //Iterative process for relaxed pressure determination
    int iteration(0);
    NewtonRaphson(numberPhases, pStar, iteration);

    //Apply the relaxation procedure only if it has converged to a solution.
    if (iteration < 100) {
      //Cell update
      for (int k = 0; k < numberPhases; k++) {
        phase = cell->getPhase(k, type);
        phase->setAlpha(TB->akS[k]);
        phase->setDensity(TB->rhokS[k]);
        phase->setPressure(pStar);
      }
      cell->getMixture(type)->setPressure(pStar);
    }
  }
}