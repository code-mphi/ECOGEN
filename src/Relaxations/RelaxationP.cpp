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

#include "RelaxationP.h"

//***********************************************************************

RelaxationP::RelaxationP(){}

//***********************************************************************

RelaxationP::~RelaxationP(){}

//***********************************************************************

void RelaxationP::NewtonRaphson(double& pStar, int& iteration)
{
  //Iterative process for relaxed pressure determination
  double f(0.), df(1.), drhodp(0.), dalpha(0.);
  do {
    iteration++;
    pStar -= f / df;
    //Physical pressure?
    for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(pStar); }
    f = -1.; df = 0.;
    for (int k = 0; k < numberPhases; k++) {
      TB->rhokS[k] = TB->eos[k]->computeDensityPfinal(TB->pk[k], TB->rhok[k], pStar, &drhodp);
      TB->akS[k] = TB->ak[k] * TB->rhok[k] / TB->rhokS[k];
      dalpha = TB->ak[k] * TB->rhok[k] * drhodp / (TB->rhokS[k] * TB->rhokS[k]);
      f += TB->akS[k];
      df -= dalpha;
    }
  } while (std::fabs(f)>1e-10 && iteration < 100);

  if (iteration == 100 && std::fabs(pStar) > 1.e-7) {
    std::stringstream warningMessage;
    warningMessage << "Not converged in RelaxationP::NewtonRaphson.";
    warningMessage << "Convergence error = " << std::fabs(f);
    warnings.push_back(Errors(warningMessage.str().c_str(), __FILE__, __LINE__));
  }
}

//***********************************************************************

double RelaxationP::computeInterfacePressure(Cell* cell, Prim type)
{
  //pI = sum_k (p_k sum_j Z_j) / (N - 1) sum_k (Z_k) ; where j is different from k and where p_k also considers solid terms
  double pI(0.), sumZk(0.), sumZj(0.);
  for (int k = 0; k < numberPhases; k++) {
    TB->zk[k] = cell->getPhase(k, type)->getDensity() * cell->getPhase(k, type)->getSoundSpeed();
  }
  for (int k = 0; k < numberPhases; k++) {
    if (!TB->alphaNull[k]) {
      sumZk += TB->zk[k];
      sumZj = 0.;
      for (int j = 0; j < numberPhases; j++) {
        if (!TB->alphaNull[j]) {
          if (j != k) sumZj += TB->zk[j];
        }
      }
      pI += sumZj * (cell->getPhase(k, type)->getPressure() - TB->compactionPk[k]);
    }
  }
  return pI /= sumZk * (numberPhases - 1);
}

//***********************************************************************