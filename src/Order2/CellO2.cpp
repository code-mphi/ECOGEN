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

#include "CellO2.h"
#include "../Models/Phase.h"

//***********************************************************************

CellO2::CellO2() : Cell(), m_vecPhasesO2(0), m_mixtureO2(0), m_vecTransportsO2(0), m_consSauvegarde(0), m_consTransportsSauvegarde(0) {}

//***********************************************************************

CellO2::CellO2(int lvl) : Cell(lvl), m_vecPhasesO2(0), m_mixtureO2(0), m_vecTransportsO2(0), m_consSauvegarde(0), m_consTransportsSauvegarde(0) {}

//***********************************************************************

CellO2::~CellO2()
{
  for (int k = 0; k < m_numberPhases; k++) {
    delete m_vecPhasesO2[k];
  }
  delete[] m_vecPhasesO2;
  if (m_vecTransportsO2 != 0) delete[] m_vecTransportsO2;
  delete m_mixtureO2;
  delete m_consSauvegarde;
	if (m_consTransportsSauvegarde != 0) delete[] m_consTransportsSauvegarde;
}

//***********************************************************************

void CellO2::allocate(const int& numberPhases, const int& numberTransports, const std::vector<AddPhys*>& addPhys, Model* model)
{
  m_numberPhases = numberPhases;
  m_numberTransports = numberTransports;
  m_vecPhases = new Phase*[numberPhases];
  m_vecPhasesO2 = new Phase*[numberPhases];
  for (int k = 0; k < numberPhases; k++){
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
  m_model = model;
}

//***********************************************************************

void CellO2::copyPhase(const int& phaseNumber, Phase* phase)
{
  m_vecPhases[phaseNumber]->copyPhase(*phase);
  m_vecPhasesO2[phaseNumber]->copyPhase(*phase);
}

//***********************************************************************

void CellO2::computeLocalSlopes(const int& numberPhases, const int& numberTransports, CellInterface& cellInterfaceRef,
  Limiter& globalLimiter, Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter,
  double& alphaCellAfterOppositeSide, double& alphaCell, double& alphaCellOtherInterfaceSide, double& epsInterface)
{
	//Mise a zero des slopes locales
	//------------------------------
  double coeff(0.), posCellInterfaceRef(0.);
	double sommeCoeff(0.), sommeCoeff2(0.);
	for (int k = 0; k < numberPhases; k++) {
		slopesPhasesLocal1[k]->setToZero();
		slopesPhasesLocal2[k]->setToZero();
	}
	slopesMixtureLocal1->setToZero();
	slopesMixtureLocal2->setToZero();
	for (int k = 0; k < numberTransports; k++) {
		slopesTransportLocal1[k] = 0.;
		slopesTransportLocal2[k] = 0.;
	}

	//Boucle sur les cell interfaces pour la determination des slopes de chaque cote de la cell
	//-----------------------------------------------------------------------------------------
  int phase0(0);
	for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
    //Calcul de la slope a gauche et a droite de la cell (AMR) en se basant sur celle de reference (cell interface a gauche ou a droite, inconnu)
    if (m_cellInterfaces[b]->getSlopesPhase(0) != 0) { //Cell interface de type CellInterfaceO2 ou BoundCondWallO2
      if (m_cellInterfaces[b] == &cellInterfaceRef) {
				for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), 1.); }
				slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), 1.);
				for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += m_cellInterfaces[b]->getSlopesTransport(k)->getValue(); }
        sommeCoeff += 1.;
      }
      else {
				if (!m_cellInterfaces[b]->getSplit()) {
				//Produit scalar des normals avec celle de reference
					coeff = std::fabs(m_cellInterfaces[b]->getFace()->getNormal().scalar(cellInterfaceRef.getFace()->getNormal()));
					if (coeff > 1.e-6) {
						//Face majoritement selon X
						if (std::fabs(cellInterfaceRef.getFace()->getNormal().getX()) > 0.5) {
							posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getX();
							//Cote cellInterfaceRef
							if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getX()) <= std::fabs(posCellInterfaceRef - m_element->getPosition().getX())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff2 += coeff;
                if (m_cellInterfaces[b]->getCellGauche() == this) {
                  if (m_cellInterfaces[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellDroite()->getPhase(phase0)->getAlpha(); } //Cell interface of type CellInterface/O2 (inner)
                  else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
							}
						}
						//Face majoritement selon Y
						else if (std::fabs(cellInterfaceRef.getFace()->getNormal().getY()) > 0.5) {
							posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getY();
							//Cote cellInterfaceRef
							if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getY()) <= std::fabs(posCellInterfaceRef - m_element->getPosition().getY())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff2 += coeff;
                if (m_cellInterfaces[b]->getCellGauche() == this) {
                  if (m_cellInterfaces[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellDroite()->getPhase(phase0)->getAlpha(); } //Cell interface of type CellInterface/O2 (inner)
                  else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
							}
						}
						//Face majoritement selon Z
						else {
							posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getZ();
							//Cote cellInterfaceRef
							if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getZ()) <= std::fabs(posCellInterfaceRef - m_element->getPosition().getZ())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff2 += coeff;
                if (m_cellInterfaces[b]->getCellGauche() == this) {
                  if (m_cellInterfaces[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellDroite()->getPhase(phase0)->getAlpha(); } //Cell interface of type CellInterface/O2 (inner)
                  else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
							}
						}
					}
				}
      }
    }
	} //fin boucle sur les cell interfaces

	//Normalisation des slopes
	//------------------------
	if (sommeCoeff > 1.e-8) {
		for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->divide(sommeCoeff);	}
		slopesMixtureLocal1->divide(sommeCoeff);
		for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] /= sommeCoeff; }
	}
	if (sommeCoeff2 > 1.e-8) {
		for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->divide(sommeCoeff2);	}
		slopesMixtureLocal2->divide(sommeCoeff2);
		for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] /= sommeCoeff2; }
    alphaCellAfterOppositeSide /= sommeCoeff2;
	}

	//Limitations des slopes
	//----------------------
  //Detection of the interface location
  if ((alphaCell >= epsInterface) && (alphaCell <= 1. - epsInterface) && ((alphaCellOtherInterfaceSide - alphaCell)*(alphaCell - alphaCellAfterOppositeSide) >= 1.e-8)) {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *slopesPhasesLocal2[k], interfaceLimiter, interfaceVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *slopesMixtureLocal2, interfaceLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = interfaceVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], slopesTransportLocal2[k]);
    }
  }
  else {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *slopesPhasesLocal2[k], globalLimiter, globalVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *slopesMixtureLocal2, globalLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = globalVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], slopesTransportLocal2[k]);
    }
  }
}

//***********************************************************************

void CellO2::computeLocalSlopesLimite(const int& numberPhases, const int& numberTransports, CellInterface& cellInterfaceRef,
  Limiter& globalLimiter, Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter,
  double& epsInterface)
{
  //Solution pour multiD cartesian (peut etre une ebauche pour le NS, a voir...)

  //Mise a zero des slopes locales
  //------------------------------
  double coeff(0.), posCellInterfaceRef(0.);
  double sommeCoeff2(0.);
  for (int k = 0; k < numberPhases; k++) {
    slopesPhasesLocal1[k]->setToZero();
    slopesPhasesLocal2[k]->setToZero();
  }
  slopesMixtureLocal1->setToZero();
  slopesMixtureLocal2->setToZero();
  for (int k = 0; k < numberTransports; k++) {
    slopesTransportLocal1[k] = 0.;
    slopesTransportLocal2[k] = 0.;
  }

  //Recupere la slope cote CL
  //-------------------------
  for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*cellInterfaceRef.getSlopesPhase(k), 1.); }
  slopesMixtureLocal1->multiplyAndAdd(*cellInterfaceRef.getSlopesMixture(), 1.);
  for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += cellInterfaceRef.getSlopesTransport(k)->getValue(); }

  //Boucle sur les cell interfaces pour la determination de la slope cote oppose a la CL
  //------------------------------------------------------------------------------------
  for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
    //Calcul de la slope a gauche et a droite de la cell (AMR) en se basant sur celle de reference (cell interface a gauche ou a droite, inconnu)
    if (m_cellInterfaces[b]->getSlopesPhase(0) != 0) { //Cell interface de type CellInterface/O2
      if (m_cellInterfaces[b] != &cellInterfaceRef) {
        if (!m_cellInterfaces[b]->getSplit()) {
          //Produit scalar des normals avec celle de reference
          coeff = std::fabs(m_cellInterfaces[b]->getFace()->getNormal().scalar(cellInterfaceRef.getFace()->getNormal()));
          if (coeff > 1.e-6) {
            //Face majoritement selon X
            if (std::fabs(cellInterfaceRef.getFace()->getNormal().getX()) > 0.5) {
              posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getX();
              //Autre cote
              if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getX()) >= std::fabs(posCellInterfaceRef - m_element->getPosition().getX())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
                sommeCoeff2 += coeff;
              }
            }
            //Face majoritement selon Y
            else if (std::fabs(cellInterfaceRef.getFace()->getNormal().getY()) > 0.5) {
              posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getY();
              //Autre cote
              if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getY()) >= std::fabs(posCellInterfaceRef - m_element->getPosition().getY())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
                sommeCoeff2 += coeff;
              }
            }
            //Face majoritement selon Z
            else {
              posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getZ();
              //Autre cote
              if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getZ()) >= std::fabs(posCellInterfaceRef - m_element->getPosition().getZ())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
                sommeCoeff2 += coeff;
              }
            }
          }
        }
      }
    }
  } //fin boucle sur les cell interfaces

  //Normalisation de la slope
  //-------------------------
  if (sommeCoeff2 > 1.e-8) {
    for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->divide(sommeCoeff2); }
    slopesMixtureLocal2->divide(sommeCoeff2);
    for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] /= sommeCoeff2; }
  }

  //Limitations des slopes
  //----------------------
  //Detection of the interface location
  int phase0(0);
  if ((m_vecPhases[phase0]->getAlpha() >= epsInterface) && (m_vecPhases[phase0]->getAlpha() <= 1. - epsInterface)) {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *slopesPhasesLocal2[k], interfaceLimiter, interfaceVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *slopesMixtureLocal2, interfaceLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = interfaceVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], slopesTransportLocal2[k]);
    }
  }
  else {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *slopesPhasesLocal2[k], globalLimiter, globalVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *slopesMixtureLocal2, globalLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = globalVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], slopesTransportLocal2[k]);
    }
  }
}

//***********************************************************************

void CellO2::saveCons(const int& numberPhases, const int& numberTransports)
{
  m_consSauvegarde->setCons(m_cons, numberPhases);
  for (int k = 0; k < numberTransports; k++) { m_consTransportsSauvegarde[k].setValue(m_consTransports[k].getValue()); }
}

//***********************************************************************

void CellO2::recuperationCons(const int& numberPhases, const int& numberTransports)
{
  m_cons->setCons(m_consSauvegarde, numberPhases);
  for (int k = 0; k < numberTransports; k++) { m_consTransports[k].setValue(m_consTransportsSauvegarde[k].getValue()); }
}

//***********************************************************************

void CellO2::predictionOrdre2(const double& dt, const int& numberPhases, const int& numberTransports, Symmetry* symmetry)
{
  m_cons->setBufferFlux(*this, numberPhases);                   //On determine Un grace au vecteur primitif des phases et de mixture que l on stocke dans fluxTempXXX
  symmetry->addSymmetricTerms(this, numberPhases);              //On ajoute dans m_cons (bilan des flux) les termes des symetries cylindrique ou spherique a partir du vecteur primitif a l'instant n
  m_cons->multiply(0.5*dt, numberPhases);                       //On multiply m_cons (bilan des flux) par dt/2
  m_cons->addFlux(1., numberPhases);                            //On y ajoute fluxTempXXX (Un) -> on obtient Un+1/2 dans m_cons
  m_cons->schemeCorrection(numberPhases);                       //Specific correction for non-conservative models (using Un in fluxTempXXX and Un+1/2 in m_cons)
  m_cons->buildPrim(m_vecPhasesO2, m_mixtureO2, numberPhases);  //On peut reconstruire m_vecPhasesO2 et m_mixtureO2 a partir de m_cons
  
  //Same process for transport (Un construction not needed)
  for (int k = 0; k < numberTransports; k++) {
    m_consTransports[k].multiply(0.5*dt);
    m_vecTransportsO2[k].setValue(m_vecTransports[k].getValue());
    m_vecTransportsO2[k].add(m_consTransports[k].getValue());
  }

  //Relaxations et correction des energies
  //FP//Derniere news, ceci doit etre exclu de cette routine car pas generique selon modele: A mettre dans Run::solveHyperbolicO2 si besoin
  m_model->relaxations(this, numberPhases, dt, vecPhasesO2);
  m_mixtureO2->totalEnergyToInternalEnergy(m_vecQuantitiesAddPhys); //On reconstruit l'energie interne a partir de l energie totale
  m_cons->correctionEnergy(this, numberPhases, vecPhasesO2);
  m_model->fulfillState(m_vecPhasesO2, m_mixtureO2, m_numberPhases);
}

//***********************************************************************

void CellO2::completeFulfillState(Prim type)
{
  //Complete thermodynamical variables
  switch (type) {
  case vecPhases: case restart: //Idem cell ordre 1
    m_model->fulfillState(m_vecPhases, m_mixture, m_numberPhases, type);
    //Extended energies depending on additional physics
    this->prepareAddPhys();
    m_mixture->internalEnergyToTotalEnergy(m_vecQuantitiesAddPhys);
    break;
  case vecPhasesO2: //Not used for now
    m_model->fulfillState(m_vecPhasesO2, m_mixtureO2, m_numberPhases);
    //Extended energies depending on additional physics
    this->prepareAddPhys();
    m_mixtureO2->internalEnergyToTotalEnergy(m_vecQuantitiesAddPhys);
    break;
  default: break;
  }
}

//***********************************************************************

void CellO2::fulfillState(Prim type)
{
  //Complete thermodynamical variables
  switch (type) {
  case vecPhases: case restart: //Identical to cell first order
    m_model->fulfillState(m_vecPhases, m_mixture, m_numberPhases);
    //This routine is used in different configurations and a short note correspond to each one:
    //- Riemann solver: No need to reconstruct the total energy there because it isn't grabbed during the Riemann problem. The total energy is directly reconstruct there.
    //The reason is to avoid calculations on the gradients of additional physics which are not necessary and furthermore wrongly computed.
    //Note that the capillary energy is not required during the Riemann problem because the models are splitted.
    //- Parallel: No need to reconstruct the total energy there because it is already communicated.
    //Note that the gradients of additional physics would also be wrongly computed if done in the ghost cells.
    //- Relaxation or correction: The total energy doesn't have to be updated there.
    break;
  case vecPhasesO2: //Only usefull for the parallel with second order
    m_model->fulfillState(m_vecPhasesO2, m_mixtureO2, m_numberPhases);
    break;
  default: break;
  }
}

//***********************************************************************

Phase* CellO2::getPhase(const int& phaseNumber, Prim type) const
{
  switch (type){
    case vecPhases: return m_vecPhases[phaseNumber]; break;
    case vecPhasesO2: return m_vecPhasesO2[phaseNumber]; break;
    default: return 0; break;
  }
}

//***********************************************************************

Phase** CellO2::getPhases(Prim type) const
{
  switch (type) {
    case vecPhases: return m_vecPhases; break;
    case vecPhasesO2: return m_vecPhasesO2; break;
    default: return 0; break;
  }
}

//***********************************************************************

Mixture* CellO2::getMixture(Prim type) const
{
  switch (type) {
    case vecPhases: return m_mixture; break;
    case vecPhasesO2: return m_mixtureO2; break;
  default: return 0; break;
  }
}

//***********************************************************************

Transport& CellO2::getTransport(const int& numTransport, Prim type) const
{
	switch (type) {
	case vecPhases: return m_vecTransports[numTransport]; break;
	case vecPhasesO2: return m_vecTransportsO2[numTransport]; break;
  default: return m_vecTransports[numTransport]; break; //FP//TODO// trouver un moyen plus intelligent de faire les renvoi par defaut sur objets.
	}
}

//***********************************************************************

Transport* CellO2::getTransports(Prim type) const
{
  switch (type) {
    case vecPhases: return m_vecTransports; break;
    case vecPhasesO2: return m_vecTransportsO2; break;
    default: return 0; break;
  }
}

//***********************************************************************

void CellO2::setTransport(double value, int& numTransport, Prim type)
{
  switch (type) {
  case vecPhases: m_vecTransports[numTransport].setValue(value); break;
  case vecPhasesO2: m_vecTransportsO2[numTransport].setValue(value); break;
  default: break;
  }
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

void CellO2::createChildCell(const int& lvl)
{
  m_childrenCells.push_back(new CellO2(lvl + 1));
}

//****************************************************************************
//********************** Methode Ordre 2 Parallele ***************************
//****************************************************************************

void CellO2::fillBufferSlopes(double* buffer, int& counter, const int& lvl, const int& neighbour) const
{
	if (m_lvl == lvl) {
    std::vector<CellInterface*> cellInterfacesWithNeighboringGhostCell;
    for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
      if (m_cellInterfaces[b]->whoAmI() == 0) { //Cell interface of type CellInterface/O2 (inner)
        if (m_cellInterfaces[b]->getLvl() == m_lvl) {
          if (this == m_cellInterfaces[b]->getCellGauche()) {
            if (m_cellInterfaces[b]->getCellDroite()->isCellGhost()) {
              if (m_cellInterfaces[b]->getCellDroite()->getRankOfNeighborCPU() == neighbour) {
                cellInterfacesWithNeighboringGhostCell.push_back(m_cellInterfaces[b]);
              }
            }
          }
          else {
            if (m_cellInterfaces[b]->getCellGauche()->isCellGhost()) {
              if (m_cellInterfaces[b]->getCellGauche()->getRankOfNeighborCPU() == neighbour) {
                cellInterfacesWithNeighboringGhostCell.push_back(m_cellInterfaces[b]);
              }
            }
          }
        }
      }
    }

    double epsilon(1.e-08), sommeCoeff(0.), scalarDiff(0.), alphaCellAfterOppositeSide(0.);
    Coord coordBuffer(0.);
    int phase0(0);
    for (unsigned int b2 = 0; b2 < cellInterfacesWithNeighboringGhostCell.size(); b2++) {
  		//Reset local slope to send
  		//-------------------------
  		sommeCoeff = 0.;
  		for (int k = 0; k < m_numberPhases; k++) {
  			slopesPhasesLocal1[k]->setToZero();
  		}
  		slopesMixtureLocal1->setToZero();
  		for (int k = 0; k < m_numberTransports; k++) {
  			slopesTransportLocal1[k] = 0.;
  		}

      //Loop over cell interfaces to determine the slope to send
      //--------------------------------------------------------
      alphaCellAfterOppositeSide = 0.;
      int slopeIndex=-1;
      for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
        if (m_cellInterfaces[b] != cellInterfacesWithNeighboringGhostCell[b2]) {
          if (m_cellInterfaces[b]->getSlopesPhase(0) != 0) { //Cell interface de type CellInterfaceO2 ou BoundCondWallO2
            if (!m_cellInterfaces[b]->getSplit()) {
              coordBuffer = m_cellInterfaces[b]->getFace()->getNormal().abs() - cellInterfacesWithNeighboringGhostCell[b2]->getFace()->getNormal();
              if (coordBuffer.norm() < epsilon) { //Face in the same direction than the reference face
                scalarDiff = m_cellInterfaces[b]->getFace()->getPos().scalar(m_cellInterfaces[b]->getFace()->getNormal()) -
                              cellInterfacesWithNeighboringGhostCell[b2]->getFace()->getPos().scalar(cellInterfacesWithNeighboringGhostCell[b2]->getFace()->getNormal());
                if (std::fabs(scalarDiff) > epsilon) { //Face on the opposite side of the cell in regards to the reference face
                  for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), 1.); }
                  slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), 1.);
                  for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
                  sommeCoeff += 1.;
                  if (m_cellInterfaces[b]->getCellGauche() == this) {
                    if (m_cellInterfaces[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellDroite()->getPhase(phase0)->getAlpha(); } //Cell interface of type CellInterface/O2 (inner)
                    else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
                  }
                  else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
                  coordBuffer=scalarDiff*cellInterfacesWithNeighboringGhostCell[b2]->getFace()->getNormal();
                  if      (coordBuffer.getX() < -epsilon) slopeIndex=0;
                  else if (coordBuffer.getX() >  epsilon) slopeIndex=1;
                  else if (coordBuffer.getY() < -epsilon) slopeIndex=2;
                  else if (coordBuffer.getY() >  epsilon) slopeIndex=3;
                  else if (coordBuffer.getZ() < -epsilon) slopeIndex=4;
                  else if (coordBuffer.getZ() >  epsilon) slopeIndex=5;
                }
              }
            }
          }
        }
      }

  		//Normalization of the slope
  		//--------------------------
  		if (sommeCoeff > 1.e-8) {
  			for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->divide(sommeCoeff); }
  			slopesMixtureLocal1->divide(sommeCoeff);
  			for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] /= sommeCoeff; }
        alphaCellAfterOppositeSide /= sommeCoeff;
  		}

  		//Fill buffer to send
  		//-------------------
  		for (int k = 0; k < m_numberPhases; k++) {
  			slopesPhasesLocal1[k]->fillBufferSlopes(buffer, counter);
  		}
  		slopesMixtureLocal1->fillBufferSlopes(buffer, counter);
  		for (int k = 0; k < m_numberTransports; k++) {
  			buffer[++counter] = slopesTransportLocal1[k];
  		}
      buffer[++counter] = alphaCellAfterOppositeSide;
      buffer[++counter] = static_cast<double>(slopeIndex);
    }
	}

	else {
    for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
      if (m_childrenCells[i]->hasNeighboringGhostCellOfCPUneighbour(neighbour)) {
        m_childrenCells[i]->fillBufferSlopes(buffer, counter, lvl, neighbour);
      }
    }
	}
}

//***********************************************************************
