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

#include <cmath>
#include <algorithm>
#include "FluxPUEq.h"
#include "../Mixture.h"

//***********************************************************************

FluxPUEq::FluxPUEq(const int& numberPhases) : FluxUEq(numberPhases){}

//***********************************************************************

FluxPUEq::~FluxPUEq(){}

//***********************************************************************

void FluxPUEq::correctionEnergy(Cell* cell, const int& numberPhases, Prim type) const
{
  Phase* phase;

  //Usefull data extraction
  for (int k = 0; k < numberPhases; k++){
    phase = cell->getPhase(k, type);
    TB->ak[k] = phase->getAlpha();
    TB->rhok[k] = phase->getDensity();
  }
  //Mixture pressure calculus from mixture EOS
  double rhoe = cell->getMixture(type)->getDensity() * cell->getMixture(type)->getEnergy();
  double p(rhoe), denom(0.), gamPinfSurGamMoinsUn(0.), eRef(0.), unSurGamMoinsUn(0.), covolume(0.);
  for (int k = 0; k < numberPhases; k++) {
    TB->eos[k]->sendSpecialMixtureEos(gamPinfSurGamMoinsUn, eRef, unSurGamMoinsUn, covolume);
    p -= TB->ak[k] * (gamPinfSurGamMoinsUn * (1. - TB->rhok[k]*covolume) + TB->rhok[k] * eRef);
    denom += TB->ak[k] * unSurGamMoinsUn * (1. - TB->rhok[k]*covolume);
  }
  p /= denom;
  
  //Puisqu'on permet au code d'avoir maintenant des pressions negatives dans le liquide, cette procedure n'est plus utile et on garantie ainsi la conservation de l'enegie.
  //En revanche, pour eviter tout probleme durant la simulation, on modifie tout de meme la pression du gas pour etre non-negative si jamais cela est le cas (depend de l'EOS).
  //Cette configuration ne devrait se produire que dans les localisations ou la fraction volumique du gas est tres faible ou nulle.
  //----------------------------------
  //Verifications for mixture pressure
  //double pCorr(p);
  //for (int k = 0; k < numberPhases; k++) {
  //  TB->eos[k]->verifyAndModifyPressure(pCorr, cell->getPhase(k, type)->getAlpha());
  //}
  //if (p < pCorr) {
  //  //Note: when doing this, we are no more locally conservative but it avoids the run to crash in some particular test cases.
  //  //Furthermore, the non-conservative correction has a very negligible impact on the total conservation of energy.
  //  //std::cout << p << " " << pCorr << cell->getMixture(type)->getPressure() << std::endl;
  //  p = cell->getMixture(type)->getPressure();
  //}
  //----------------------------------

  //Cell update
  for (int k = 0; k < numberPhases; k++){
    phase = cell->getPhase(k, type);
    phase->setPressure(p);
    //Same remark as previous one (instead of "phase->verifyPhase("correctionEnergy: ");").
    phase->verifyAndCorrectPhase();
  }
  cell->getMixture(type)->setPressure(p);
}

//***********************************************************************