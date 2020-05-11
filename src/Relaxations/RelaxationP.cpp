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

//! \file      RelaxationP.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "RelaxationP.h"

//***********************************************************************

RelaxationP::RelaxationP(){}

//***********************************************************************

RelaxationP::~RelaxationP(){}

//***********************************************************************

void RelaxationP::stiffRelaxation(Cell *cell, const int &numberPhases, Prim type) const
{
  //Is the pressure-relaxation procedure necessary?
  //If alpha = 0 is activated, a test is done to know if the relaxation procedure is necessary or not
  //Else, i.e. alpha = 0 is desactivated (alpha != 0), the relaxation procedure is always done (relax = true)
  bool relax(true);
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    for (int k = 0; k < numberPhases; k++) {
      if (cell->getPhase(k, type)->getAlpha() >(1. - 1.e-5)) relax = false;
    }
  }

  if (relax) {
    double drho(0.), dalpha(0.);
    Phase *phase(0);
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
    double f(0.), df(1.);
    do {
      iteration++;
      pStar -= f / df;
      //Physical pressure?
      for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(pStar); }
      f = -1.; df = 0.;
      for (int k = 0; k < numberPhases; k++) {
        TB->rhokS[k] = TB->eos[k]->computeDensityPfinal(TB->pk[k], TB->rhok[k], pStar, &drho);
        TB->akS[k] = TB->ak[k] * TB->rhok[k] / TB->rhokS[k];
        dalpha = TB->ak[k] * TB->rhok[k] * drho / (TB->rhokS[k] * TB->rhokS[k]);
        f += TB->akS[k];
        df -= dalpha;
      }
      if (iteration > 100) {
        std::cout << "pStar=" << pStar << " f=" << f << " df=" << df << std::endl;
        errors.push_back(Errors("Not converged in relaxPressures", __FILE__, __LINE__));
        break;
      }
    } while (std::fabs(f)>1e-10 && iteration < 100);
    //} while (std::fabs(f) > 1e-10);

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