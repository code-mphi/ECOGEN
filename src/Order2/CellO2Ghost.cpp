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

#include "CellO2Ghost.h"

//***********************************************************************

CellO2Ghost::CellO2Ghost() : CellO2(), m_rankOfNeighborCPU(0), m_vecPhasesSlopesGhost(0), m_mixtureSlopesGhost(0), m_vecTransportsSlopesGhost(0) {}

//***********************************************************************

CellO2Ghost::CellO2Ghost(int lvl) : CellO2(lvl), m_rankOfNeighborCPU(0), m_vecPhasesSlopesGhost(0), m_mixtureSlopesGhost(0), m_vecTransportsSlopesGhost(0) {}

//***********************************************************************

CellO2Ghost::~CellO2Ghost()
{
	for (unsigned int s = 0; s < m_vecPhasesSlopesGhost.size(); s++) {
		for (int k = 0; k < numberPhases; k++) {
			delete m_vecPhasesSlopesGhost[s][k];
		}
		delete[] m_vecPhasesSlopesGhost[s];
		delete m_mixtureSlopesGhost[s];
		delete[] m_vecTransportsSlopesGhost[s];
	}
}

//***********************************************************************

void CellO2Ghost::pushBackSlope()
{
  m_indexCellInterface.push_back(0);
  m_vecPhasesSlopesGhost.push_back(nullptr);
  m_mixtureSlopesGhost.push_back(nullptr);
  m_vecTransportsSlopesGhost.push_back(nullptr);
  m_alphaCellAfterOppositeSide.push_back(0.);
}

//***********************************************************************

void CellO2Ghost::allocate(const std::vector<AddPhys*>& addPhys)
{
  m_vecPhases = new Phase*[numberPhases];
  m_vecPhasesO2 = new Phase*[numberPhases];
  for (int k = 0; k < numberPhases; k++) {
    model->allocatePhase(&m_vecPhases[k]);
    model->allocatePhase(&m_vecPhasesO2[k]);
  }
  model->allocateMixture(&m_mixture);
  model->allocateMixture(&m_mixtureO2);
  model->allocateCons(&m_cons);
  model->allocateCons(&m_consSauvegarde);
  if (numberTransports > 0) {
    m_vecTransports = new Transport[numberTransports];
    m_consTransports = new Transport[numberTransports];
    m_consTransportsSauvegarde = new Transport[numberTransports];
    m_vecTransportsO2 = new Transport[numberTransports];
  }
  for (unsigned int k = 0; k < addPhys.size(); k++) { addPhys[k]->addQuantityAddPhys(this); }

  //Allocation des slopes fantomes, specifique aux limites paralleles
  for (unsigned int s = 0; s < m_vecPhasesSlopesGhost.size(); s++) {
    m_vecPhasesSlopesGhost[s] = new Phase*[numberPhases];
    for (int k = 0; k < numberPhases; k++) {
      model->allocatePhase(&m_vecPhasesSlopesGhost[s][k]);
      m_vecPhasesSlopesGhost[s][k]->setToZero();
    }
    model->allocateMixture(&m_mixtureSlopesGhost[s]);
    m_mixtureSlopesGhost[s]->setToZero();
    m_vecTransportsSlopesGhost[s] = new double[numberTransports];
    for (int k = 0; k < numberTransports; k++) { m_vecTransportsSlopesGhost[s][k] = 0.; }
  }
}

//***************************************************************************

int CellO2Ghost::getRankOfNeighborCPU() const
{
  return m_rankOfNeighborCPU;
}

//***************************************************************************

void CellO2Ghost::setRankOfNeighborCPU(int rank)
{
  m_rankOfNeighborCPU = rank;
}

//***********************************************************************

void CellO2Ghost::computeLocalSlopes(CellInterface& cellInterfaceRef, Limiter& globalLimiter, Limiter& interfaceLimiter,
	Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter,
	double& alphaCellAfterOppositeSide, double& alphaCell, double& alphaCellOtherInterfaceSide, double& epsInterface)
{
	//Find the corresponding slopes store inside this ghost cell
	//----------------------------------------------------------
	int s(-1);
	int refIndex=-1;
	double epsilon(1.e-08), scalarDiff(0.);
    Coord coordBuffer(0.);
	scalarDiff = m_element->getPosition().scalar(cellInterfaceRef.getFace()->getNormal()) -
	     cellInterfaceRef.getFace()->getPos().scalar(cellInterfaceRef.getFace()->getNormal());
	coordBuffer=scalarDiff*cellInterfaceRef.getFace()->getNormal();
	if      (coordBuffer.getX() < -epsilon) refIndex=0;
	else if (coordBuffer.getX() >  epsilon) refIndex=1;
	else if (coordBuffer.getY() < -epsilon) refIndex=2;
	else if (coordBuffer.getY() >  epsilon) refIndex=3;
	else if (coordBuffer.getZ() < -epsilon) refIndex=4;
	else if (coordBuffer.getZ() >  epsilon) refIndex=5;

	for (unsigned int b = 0; b < m_indexCellInterface.size(); b++) {
		if (m_indexCellInterface[b] == refIndex) {
			s = b;
			break;
		}
	}

	if (s == -1) {
		refIndex = -1;
		for (unsigned int b = 0; b < m_indexCellInterface.size(); b++) {
			if (m_indexCellInterface[b] == refIndex) {
				s = b;
				break;
			}
		}
	}

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

	//Boucle sur les cell interfaces pour la determination de la slope du cote de cellInterfaceRef
	//--------------------------------------------------------------------------------------------
	coordBuffer = 0.;
	for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
		if (!m_cellInterfaces[b]->getSplit()) {
			coordBuffer = m_cellInterfaces[b]->getFace()->getNormal().abs() - cellInterfaceRef.getFace()->getNormal();
            if (coordBuffer.norm() < epsilon) { //Face in the same direction than the reference face
				for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), 1.); }
				slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), 1.);
				for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += m_cellInterfaces[b]->getSlopesTransport(k)->getValue(); }
				sommeCoeff += 1.;
			}
		}
	}

	//Normalisation des slopes
	//------------------------
	if (sommeCoeff > 1.e-8) {
		for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->divide(sommeCoeff); }
		slopesMixtureLocal1->divide(sommeCoeff);
		for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] /= sommeCoeff; }
	}

	//Limitations des slopes
	//----------------------
	if (m_indexCellInterface[s] != -1) { //Slope stores in ghost cell is from a cell interface of type cellInterfaceO2 or BoundCondWallO2
		alphaCellAfterOppositeSide = m_alphaCellAfterOppositeSide[s]; //Detection of the interface and THINC method are simplified in parallel
	}
	else {
		alphaCellAfterOppositeSide = m_vecPhases[0]->getAlpha();
	}
	if ((alphaCell >= epsInterface) && (alphaCell <= 1. - epsInterface) &&
	   ((alphaCellOtherInterfaceSide - alphaCell)*(alphaCell - alphaCellAfterOppositeSide) >= 1.e-8)) {
		for (int k = 0; k < numberPhases; k++) {
			slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *m_vecPhasesSlopesGhost[s][k], interfaceLimiter, interfaceVolumeFractionLimiter);
		}
		slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *m_mixtureSlopesGhost[s], interfaceLimiter);
		for (int k = 0; k < numberTransports; k++) {
			slopesTransportLocal1[k] = interfaceVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], m_vecTransportsSlopesGhost[s][k]);
		}
	}
	else {
		for (int k = 0; k < numberPhases; k++) {
			slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *m_vecPhasesSlopesGhost[s][k], globalLimiter, globalVolumeFractionLimiter);
		}
		slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *m_mixtureSlopesGhost[s], globalLimiter);
		for (int k = 0; k < numberTransports; k++) {
			slopesTransportLocal1[k] = globalVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], m_vecTransportsSlopesGhost[s][k]);
		}
	}
}

//***********************************************************************

void CellO2Ghost::createChildCell(const int& lvl)
{
	m_childrenCells.push_back(new CellO2Ghost(lvl + 1));
	m_childrenCells.back()->setRankOfNeighborCPU(m_rankOfNeighborCPU);
	m_childrenCells.back()->pushBackSlope();
}

//***********************************************************************

void CellO2Ghost::getBufferSlopes(double* buffer, int& counter, const int& lvl)
{
	if (m_lvl == lvl) {
		for (unsigned int s = 0; s < m_vecPhasesSlopesGhost.size(); s++) {
			for (int k = 0; k < numberPhases; k++) {
				m_vecPhasesSlopesGhost[s][k]->getBufferSlopes(buffer, counter);
			}
			m_mixtureSlopesGhost[s]->getBufferSlopes(buffer, counter);
			for (int k = 0; k < numberTransports; k++) {
				m_vecTransportsSlopesGhost[s][k] = buffer[++counter];
			}
	    	m_alphaCellAfterOppositeSide[s] = buffer[++counter];
	    	m_indexCellInterface[s] = static_cast<int>(std::round(buffer[++counter]));
	    }
	}
	else {
		for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
			m_childrenCells[i]->getBufferSlopes(buffer, counter, lvl);
		}
	}
}

//***********************************************************************