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

//! \file      CellO2Ghost.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      July 30 2018

#include "CellO2Ghost.h"
#include <iostream>

using namespace std;

//***********************************************************************

CellO2Ghost::CellO2Ghost() : CellO2(), m_vecPhasesSlopesGhost(0), m_mixtureSlopesGhost(0), m_vecTransportsSlopesGhost(0) {}

//***********************************************************************

CellO2Ghost::CellO2Ghost(int lvl) : CellO2(lvl), m_vecPhasesSlopesGhost(0), m_mixtureSlopesGhost(0), m_vecTransportsSlopesGhost(0) {}

//***********************************************************************

CellO2Ghost::~CellO2Ghost()
{
	for (int k = 0; k < m_numberPhases; k++) {
		delete m_vecPhasesSlopesGhost[k];
	}
	delete[] m_vecPhasesSlopesGhost;
	delete m_mixtureSlopesGhost;
	delete[] m_vecTransportsSlopesGhost;
}

//***********************************************************************

void CellO2Ghost::allocate(const int &numberPhases, const int &numberTransports, const std::vector<AddPhys*> &addPhys, Model *model)
{
	m_numberPhases = numberPhases;
	m_numberTransports = numberTransports;
	m_vecPhases = new Phase*[numberPhases];
	m_vecPhasesO2 = new Phase*[numberPhases];
	for (int k = 0; k < numberPhases; k++) {
		model->allocatePhase(&m_vecPhases[k]);
		model->allocatePhase(&m_vecPhasesO2[k]);
	}
	model->allocateMixture(&m_mixture);
	model->allocateMixture(&m_mixtureO2);
	model->allocateCons(&m_cons, numberPhases);
	model->allocateCons(&m_consSauvegarde, numberPhases);
	if (numberTransports > 0) {
		m_vecTransports = new Transport[numberTransports];
		m_consTransports = new Transport[numberTransports];
		m_consTransportsSauvegarde = new Transport[numberTransports];
		m_vecTransportsO2 = new Transport[numberTransports];
	}
	for (unsigned int k = 0; k < addPhys.size(); k++) {
    addPhys[k]->addQuantityAddPhys(this);
	}

	//Allocation des slopes fantomes, specifique aux limites paralleles
	m_vecPhasesSlopesGhost = new Phase*[numberPhases];
	for (int k = 0; k < numberPhases; k++) {
		model->allocatePhase(&m_vecPhasesSlopesGhost[k]);
		m_vecPhasesSlopesGhost[k]->setToZero();
	}
	model->allocateMixture(&m_mixtureSlopesGhost);
	m_mixtureSlopesGhost->setToZero();
	m_vecTransportsSlopesGhost = new double[numberTransports];
	for (int k = 0; k < numberTransports; k++) {
		m_vecTransportsSlopesGhost[k] = 0.;
	}
  m_model = model;
}

//***********************************************************************

void CellO2Ghost::computeLocalSlopes(const int &numberPhases, const int &numberTransports, CellInterface &bordRef, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, double &alphaCellAfterOppositeSide, double &alphaCell, double &alphaCellOtherInterfaceSide, double &epsInterface)
{
	//Mise a zero des slopes locales
	//------------------------------
	double sommeCoeff(0.);
	for (int k = 0; k < numberPhases; k++) {
		slopesPhasesLocal1[k]->setToZero();
	}
	slopesMixtureLocal1->setToZero();
	for (int k = 0; k < numberTransports; k++) {
		slopesTransportLocal1[k] = 0.;
	}

	//Boucle sur les boundaries pour la détermination de la slope du cote de bordRef
	//-------------------------------------------------------------------------
	for (unsigned int b = 0; b < m_boundaries.size(); b++) {
		if (!m_boundaries[b]->getSplit()) {
			for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
			slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
			for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += m_boundaries[b]->getSlopesTransport(k)->getValue(); }
			sommeCoeff += 1.;
		}
	}

	//Normalisation des slopes
	//------------------------
	if (sommeCoeff > 1.e-5) {
		for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->divide(sommeCoeff); }
		slopesMixtureLocal1->divide(sommeCoeff);
		for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] /= sommeCoeff; }
	}
  
	//Limitations des slopes
	//----------------------
  alphaCellAfterOppositeSide = m_alphaCellAfterOppositeSide; //Detection of the interface and THINC method are simplified in parallel
  if ((alphaCell >= epsInterface) && (alphaCell <= 1. - epsInterface) && ((alphaCellOtherInterfaceSide - alphaCell)*(alphaCell - alphaCellAfterOppositeSide) >= 1.e-8)) {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *m_vecPhasesSlopesGhost[k], interfaceLimiter, interfaceVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *m_mixtureSlopesGhost, interfaceLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = interfaceVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], m_vecTransportsSlopesGhost[k]);
    }
  }
  else {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *m_vecPhasesSlopesGhost[k], globalLimiter, globalVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *m_mixtureSlopesGhost, globalLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = globalVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], m_vecTransportsSlopesGhost[k]);
    }
  }
}

//***********************************************************************

void CellO2Ghost::createChildCell(const int &num, const int &lvl)
{
	m_childrenCells.push_back(new CellO2Ghost(lvl + 1));
}

//***********************************************************************

void CellO2Ghost::getBufferSlopes(double *buffer, int &counter)
{
	for (int k = 0; k < m_numberPhases; k++) {
		m_vecPhasesSlopesGhost[k]->getBufferSlopes(buffer, counter);
	}
	m_mixtureSlopesGhost->getBufferSlopes(buffer, counter);
	for (int k = 0; k < m_numberTransports; k++) {
		m_vecTransportsSlopesGhost[k] = buffer[++counter];
	}
  m_alphaCellAfterOppositeSide = buffer[++counter];
}

//***********************************************************************

void CellO2Ghost::getBufferSlopesAMR(double *buffer, int &counter, const int &lvl)
{
	if (m_lvl == lvl) {
		for (int k = 0; k < m_numberPhases; k++) {
			m_vecPhasesSlopesGhost[k]->getBufferSlopes(buffer, counter);
		}
		m_mixtureSlopesGhost->getBufferSlopes(buffer, counter);
		for (int k = 0; k < m_numberTransports; k++) {
			m_vecTransportsSlopesGhost[k] = buffer[++counter];
		}
    m_alphaCellAfterOppositeSide = buffer[++counter];
	}
	else {
		for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
			m_childrenCells[i]->getBufferSlopesAMR(buffer, counter, lvl);
		}
	}
}

//***********************************************************************