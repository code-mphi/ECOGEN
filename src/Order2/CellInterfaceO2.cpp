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

#include "CellInterfaceO2.h"

Phase** slopesPhasesLocal1;
Phase** slopesPhasesLocal2;
Mixture* slopesMixtureLocal1;
Mixture* slopesMixtureLocal2;
double* slopesTransportLocal1;
double* slopesTransportLocal2;

//***********************************************************************

CellInterfaceO2::CellInterfaceO2() : CellInterface()
{}

//***********************************************************************

CellInterfaceO2::CellInterfaceO2(int lvl) : CellInterface(lvl)
{}

//***********************************************************************

CellInterfaceO2::~CellInterfaceO2()
{
}

//***********************************************************************

void CellInterfaceO2::allocateSlopes(int& allocateSlopeLocal)
{  
  // Allocate extern variables
  if (allocateSlopeLocal < 1) {
    slopesPhasesLocal1 = new Phase*[numberPhases];
    slopesPhasesLocal2 = new Phase*[numberPhases];
    for (int k = 0; k < numberPhases; k++) {
      m_cellLeft->getPhase(k)->allocateAndCopyPhase(&slopesPhasesLocal1[k]);
      m_cellLeft->getPhase(k)->allocateAndCopyPhase(&slopesPhasesLocal2[k]);
      slopesPhasesLocal1[k]->setToZero();
      slopesPhasesLocal2[k]->setToZero();
    }

    m_cellLeft->getMixture()->allocateAndCopyMixture(&slopesMixtureLocal1);
    m_cellLeft->getMixture()->allocateAndCopyMixture(&slopesMixtureLocal2);
    slopesMixtureLocal1->setToZero();
    slopesMixtureLocal2->setToZero();

		slopesTransportLocal1 = new double[numberTransports];
		slopesTransportLocal2 = new double[numberTransports];
		for (int k = 0; k < numberTransports; k++) {
			slopesTransportLocal1[k] = 0.;
			slopesTransportLocal2[k] = 0.;
		}

    allocateSlopeLocal = 1;
  }
}

//***********************************************************************

void CellInterfaceO2::computeFlux(double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter,
 Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type)
{
  // Quand on fait le premier computeFlux (donc avec vecPhases) on n'incremente pas m_cons pour les cells de niveau different (inferieur) de "lvl".
  // Sinon ca veut dire qu on l ajoute pour les 2 computeFlux sans le remettre a zero entre les deux, donc 2 fois plus de flux que ce que l on veut.
  this->solveRiemann(dtMax, globalLimiter, interfaceLimiter, globalVolumeFractionLimiter, interfaceVolumeFractionLimiter, type);

  switch (type) {
  case vecPhases:
    if (m_cellLeft->getLvl() == m_cellRight->getLvl()) {       //CoefAMR = 1 pour les deux
      this->addFlux(1.);       //Add of flux on right cell
      this->subtractFlux(1.);  //Subtract of flux on left cell
    }
    else if (m_cellLeft->getLvl() > m_cellRight->getLvl()) {   //CoefAMR = 1 pour la gauche et on n'ajoute rien sur la right cell
      this->subtractFlux(1.);  //Subtract of flux on left cell
    }
    else {                                                     //CoefAMR = 1 pour la droite et on ne retire rien sur la left cell
      this->addFlux(1.);       //Add of flux on right cell
    }
    break;

  case vecPhasesO2:
    if (m_cellLeft->getLvl() == m_cellRight->getLvl()) {       //CoefAMR = 1 pour les deux
      this->addFlux(1.);       //Add of flux on right cell
      this->subtractFlux(1.);  //Subtract of flux on left cell
    }
    else if (m_cellLeft->getLvl() > m_cellRight->getLvl()) {   //CoefAMR = 1 pour la gauche et 0.5 pour la droite
      this->addFlux(0.5);      //Add of flux on right cell
      this->subtractFlux(1.);  //Subtract of flux on left cell
    }
    else {                                                     //CoefAMR = 1 pour la droite et 0.5 pour la gauche
      this->addFlux(1.);       //Add of flux on right cell
      this->subtractFlux(0.5); //Subtract of flux on left cell
    }
    break;

  default: break;
  }
}

//***********************************************************************