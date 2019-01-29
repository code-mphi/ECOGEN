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

//! \file      CellO2.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      July 30 2018

#include "CellO2.h"
#include <iostream>
#include "../Models/Phase.h"

using namespace std;

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

void CellO2::allocate(const int &numberPhases, const int &numberTransports, const std::vector<AddPhys*> &addPhys, Model *model)
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

void CellO2::copyPhase(const int &phaseNumber, Phase *phase)
{
  m_vecPhases[phaseNumber]->copyPhase(*phase);
  m_vecPhasesO2[phaseNumber]->copyPhase(*phase);
}

//***********************************************************************

void CellO2::computeLocalSlopes(const int &numberPhases, const int &numberTransports, CellInterface &bordRef, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, double &alphaCellAfterOppositeSide, double &alphaCell, double &alphaCellOtherInterfaceSide, double &epsInterface)
{
	//Solution pour multiD cartesian (peut etre une ebauche pour le NS, a voir...)

	//Mise a zero des slopes locales
	//------------------------------
  double coeff(0.), posBordRef(0.);
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

	//Boucle sur les boundaries pour la détermination des slopes de chaque cote de la cell
	//----------------------------------------------------------------------------------
  int phase0(0);
	for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    //Calcul de la slope a gauche et a droite de la cell (AMR) en se basant sur celle de reference (bord a gauche ou a droite, inconnu)
    if (m_boundaries[b]->getSlopesPhase(0) != 0) { //Bord de type CellInterfaceO2, BoundCondWallO2 ou BoundCondSymmetryO2
      if (m_boundaries[b] == &bordRef) {
				for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
				slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
				for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += m_boundaries[b]->getSlopesTransport(k)->getValue(); }
        sommeCoeff += 1.;
      }
      else {
				if (!m_boundaries[b]->getSplit()) {
				//Produit scalar des normals avec celle de reference
					coeff = abs(m_boundaries[b]->getFace()->getNormal().scalar(bordRef.getFace()->getNormal()));
					if (coeff > 1.e-6) {
						//Face majoritement selon X
						if (abs(bordRef.getFace()->getNormal().getX()) > 0.5) {
							posBordRef = bordRef.getFace()->getPos().getX();
							//Cote bordRef
							if (abs(posBordRef - m_boundaries[b]->getFace()->getPos().getX()) <= abs(posBordRef - m_element->getPosition().getX())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff2 += coeff;
                if (m_boundaries[b]->getCellGauche() == this) {
                  if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                  else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
							}
						}
						//Face majoritement selon Y
						else if (abs(bordRef.getFace()->getNormal().getY()) > 0.5) {
							posBordRef = bordRef.getFace()->getPos().getY();
							//Cote bordRef
							if (abs(posBordRef - m_boundaries[b]->getFace()->getPos().getY()) <= abs(posBordRef - m_element->getPosition().getY())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff2 += coeff;
                if (m_boundaries[b]->getCellGauche() == this) {
                  if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                  else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
							}
						}
						//Face majoritement selon Z
						else {
							posBordRef = bordRef.getFace()->getPos().getZ();
							//Cote bordRef
							if (abs(posBordRef - m_boundaries[b]->getFace()->getPos().getZ()) <= abs(posBordRef - m_element->getPosition().getZ())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
								sommeCoeff2 += coeff;
                if (m_boundaries[b]->getCellGauche() == this) {
                  if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                  else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
							}
						}
					}
				}
      }
    }
	} //fin boucle sur les boundaries

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

void CellO2::computeLocalSlopesLimite(const int &numberPhases, const int &numberTransports, CellInterface &bordRef, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter)
{
  //Solution pour multiD cartesian (peut etre une ebauche pour le NS, a voir...)

  //Mise a zero des slopes locales
  //------------------------------
  double coeff(0.), posBordRef(0.);
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
  for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*bordRef.getSlopesPhase(k), 1.); }
  slopesMixtureLocal1->multiplyAndAdd(*bordRef.getSlopesMixture(), 1.);
  for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += bordRef.getSlopesTransport(k)->getValue(); }

  //Boucle sur les boundaries pour la détermination de la slope cote oppose a la CL
  //--------------------------------------------------------------------------
  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    //Calcul de la slope a gauche et a droite de la cell (AMR) en se basant sur celle de reference (bord a gauche ou a droite, inconnu)
    if (m_boundaries[b]->getSlopesPhase(0) != 0) { //Bord de type CellInterface/O2
      if (m_boundaries[b] != &bordRef) {
        if (!m_boundaries[b]->getSplit()) {
          //Produit scalar des normals avec celle de reference
          coeff = abs(m_boundaries[b]->getFace()->getNormal().scalar(bordRef.getFace()->getNormal()));
          if (coeff > 1.e-6) {
            //Face majoritement selon X
            if (abs(bordRef.getFace()->getNormal().getX()) > 0.5) {
              posBordRef = bordRef.getFace()->getPos().getX();
              //Autre cote
              if (abs(posBordRef - m_boundaries[b]->getFace()->getPos().getX()) >= abs(posBordRef - m_element->getPosition().getX())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
                sommeCoeff2 += coeff;
              }
            }
            //Face majoritement selon Y
            else if (abs(bordRef.getFace()->getNormal().getY()) > 0.5) {
              posBordRef = bordRef.getFace()->getPos().getY();
              //Autre cote
              if (abs(posBordRef - m_boundaries[b]->getFace()->getPos().getY()) >= abs(posBordRef - m_element->getPosition().getY())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
                sommeCoeff2 += coeff;
              }
            }
            //Face majoritement selon Z
            else {
              posBordRef = bordRef.getFace()->getPos().getZ();
              //Autre cote
              if (abs(posBordRef - m_boundaries[b]->getFace()->getPos().getZ()) >= abs(posBordRef - m_element->getPosition().getZ())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_boundaries[b]->getSlopesTransport(k)->getValue()); }
                sommeCoeff2 += coeff;
              }
            }
          }
        }
      }
    }
  } //fin boucle sur les boundaries

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
  double epsInterface(1.e-5);
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

void CellO2::saveCons(const int &numberPhases, const int &numberTransports)
{
  m_consSauvegarde->setCons(m_cons, numberPhases);
  for (int k = 0; k < numberTransports; k++) { m_consTransportsSauvegarde[k].setValue(m_consTransports[k].getValue()); }
}

//***********************************************************************

void CellO2::recuperationCons(const int &numberPhases, const int &numberTransports)
{
  m_cons->setCons(m_consSauvegarde, numberPhases);
  for (int k = 0; k < numberTransports; k++) { m_consTransports[k].setValue(m_consTransportsSauvegarde[k].getValue()); }
}

//***********************************************************************

void CellO2::predictionOrdre2(const double &dt, const int &numberPhases, const int &numberTransports, Symmetry *symmetry)
{
  m_cons->setBufferFlux(*this, numberPhases);                   //On determine Un grace au vecteur primitif des phases et de mixture que l on stocke dans fluxTempXXX
  symmetry->addSymmetricTerms(this, numberPhases);              //On ajoute dans m_cons (bilan des flux) les termes des symetries cylindrique ou spherique a partir du vecteur primitif a l'instant n
  m_cons->multiply(0.5*dt, numberPhases);                       //On multiply m_cons (bilan des flux) par dt/2
  m_cons->addFlux(1., numberPhases);                            //On y ajoute fluxTempXXX (Un) -> on obtient Un+1/2 dans m_cons
  m_cons->schemeCorrection(this, numberPhases);
  m_cons->buildPrim(m_vecPhasesO2, m_mixtureO2, numberPhases);  //On peut reconstruire m_vecPhasesO2 et m_mixtureO2 a partir de m_cons
  
  //Same process for transport (Un construction not needed)
  for (int k = 0; k < numberTransports; k++) {
    m_consTransports[k].multiply(0.5*dt);
    m_vecTransportsO2[k].setValue(m_vecTransports[k].getValue());
    m_vecTransportsO2[k].add(m_consTransports[k].getValue());
  }

  //Relaxations et correction des energies
  m_model->relaxations(this, numberPhases, vecPhasesO2);
  m_mixtureO2->totalEnergyToInternalEnergy(m_vecQuantitiesAddPhys); //On reconstruit l'energie interne a partir de l energie totale
  m_cons->correctionEnergy(this, numberPhases, vecPhasesO2);
  m_model->fulfillState(m_vecPhasesO2, m_mixtureO2, m_numberPhases);
}

//***********************************************************************

void CellO2::allocateAndCopyPhase(const int &phaseNumber, Phase *phase)
{
  phase->allocateAndCopyPhase(&m_vecPhases[phaseNumber]);
  phase->allocateAndCopyPhase(&m_vecPhasesO2[phaseNumber]);
}

//***********************************************************************

void CellO2::completeFulfillState(Prim type)
{
  //Complete thermodynamical variables
  switch (type) {
  case vecPhases: case resume: //Idem cell ordre 1
    m_model->fulfillState(m_vecPhases, m_mixture, m_numberPhases, type);
    //Extended energies depending on additional physics
    this->prepareAddPhys();
    m_mixture->internalEnergyToTotalEnergy(m_vecQuantitiesAddPhys);
    break;
  case vecPhasesO2: //Utile seulement pour le parallele ordre 2
    cout << "test: are we passing through here sometime?" << endl; 
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
  case vecPhases: case resume: //Identical to cell first order
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

void CellO2::localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal, const int &numberPhases, Prim type)
{
  switch (type) {
    case vecPhases:
      for (int k = 0; k < numberPhases; k++) {
        m_vecPhases[k]->localProjection(normal, tangent, binormal);
      }
      m_mixture->localProjection(normal, tangent, binormal);
      break;
    case vecPhasesO2:
      for (int k = 0; k < numberPhases; k++) {
        m_vecPhasesO2[k]->localProjection(normal, tangent, binormal);
      }
      m_mixtureO2->localProjection(normal, tangent, binormal);
      break;
    default: break;
  }
}

//***********************************************************************

void CellO2::copyInCell(Cell &cellSource, Prim type) const
{
  cellSource = static_cast<Cell>(*this);
}

//***********************************************************************

Phase* CellO2::getPhase(const int &phaseNumber, Prim type) const
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

Transport& CellO2::getTransport(const int &numTransport, Prim type) const
{
	switch (type) {
	case vecPhases: return m_vecTransports[numTransport]; break;
	case vecPhasesO2: return m_vecTransportsO2[numTransport]; break;
  default: return m_vecTransports[numTransport]; break;
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

void CellO2::setTransport(double value, int &numTransport, Prim type)
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

void CellO2::createChildCell(const int &num, const int &lvl)
{
  m_childrenCells.push_back(new CellO2(lvl + 1));
}

//****************************************************************************
//********************** Methode Ordre 2 Parallele ***************************
//****************************************************************************

void CellO2::fillBufferSlopes(double *buffer, int &counter, string whichCpuAmIForNeighbour) const
{
  //Mise a zero de la slope locale a transmettre
  //--------------------------------------------
  double epsilon(1.e-15), sommeCoeff(0.);
  for (int k = 0; k < m_numberPhases; k++) {
    slopesPhasesLocal1[k]->setToZero();
  }
  slopesMixtureLocal1->setToZero();
  for (int k = 0; k < m_numberTransports; k++) {
    slopesTransportLocal1[k] = 0.;
  }

  //Boucle sur les boundaries pour la détermination de la slope a transmettre
  //--------------------------------------------------------------------
  int phase0(0);
  double alphaCellAfterOppositeSide(0.);
  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    if (!m_boundaries[b]->getSplit()) {
      //Normal of the face in the x-direction
      if (abs(m_boundaries[b]->getFace()->getNormal().getX()) > 0.99) {
        //I am left CPU, I send the average slope of the left side of the right cells
        if (whichCpuAmIForNeighbour == "LEFT") {
          if ((m_boundaries[b]->getFace()->getPos().getX() + epsilon) < m_element->getPosition().getX()) {
            for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
            slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
            for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
            sommeCoeff += 1.;
            if (m_boundaries[b]->getCellGauche() == this) {
              if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
            else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
          }
        }
        //I am right CPU, I send the average slope of the right side of the left cells
        else if (whichCpuAmIForNeighbour == "RIGHT") {
          if ((m_boundaries[b]->getFace()->getPos().getX() - epsilon) > m_element->getPosition().getX()) {
            for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
            slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
            for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
            sommeCoeff += 1.;
            if (m_boundaries[b]->getCellGauche() == this) {
              if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
            else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
          }
        }
      }
      //Normal of the face in the y-direction
      else if (abs(m_boundaries[b]->getFace()->getNormal().getY()) > 0.99) {
        //I am bottom CPU, I send the average slope of the bottom side of the top cells
        if (whichCpuAmIForNeighbour == "BOTTOM") {
          if ((m_boundaries[b]->getFace()->getPos().getY() + epsilon) < m_element->getPosition().getY()) {
            for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
            slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
            for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
            sommeCoeff += 1.;
            if (m_boundaries[b]->getCellGauche() == this) {
              if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
            else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
          }
        }
        //I am top CPU, I send the average slope of the top side of the bottom cells
        else if (whichCpuAmIForNeighbour == "TOP") {
          if ((m_boundaries[b]->getFace()->getPos().getY() - epsilon) > m_element->getPosition().getY()) {
            for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
            slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
            for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
            sommeCoeff += 1.;
            if (m_boundaries[b]->getCellGauche() == this) {
              if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
            else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
          }
        }
      }
      //Normal of the face in the z-direction
      else if (abs(m_boundaries[b]->getFace()->getNormal().getZ()) > 0.99) {
        //I am back CPU, I send the average slope of the back side of the front cells
        if (whichCpuAmIForNeighbour == "BACK") {
          if ((m_boundaries[b]->getFace()->getPos().getZ() + epsilon) < m_element->getPosition().getZ()) {
            for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
            slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
            for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
            sommeCoeff += 1.;
            if (m_boundaries[b]->getCellGauche() == this) {
              if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
            else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
          }
        }
        //I am front CPU, I send the average slope of the front side of the back cells
        else if (whichCpuAmIForNeighbour == "FRONT") {
          if ((m_boundaries[b]->getFace()->getPos().getZ() - epsilon) > m_element->getPosition().getZ()) {
            for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
            slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
            for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
            sommeCoeff += 1.;
            if (m_boundaries[b]->getCellGauche() == this) {
              if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
            else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
          }
        }
      }
    }
  }

  //Normalisation de la slope
  //-------------------------
  if (sommeCoeff > 1.e-5) {
    for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->divide(sommeCoeff); }
    slopesMixtureLocal1->divide(sommeCoeff);
    for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] /= sommeCoeff; }
    alphaCellAfterOppositeSide /= sommeCoeff;
  }

  //Rempli buffer pour envoi
  //------------------------
  for (int k = 0; k < m_numberPhases; k++) {
    slopesPhasesLocal1[k]->fillBufferSlopes(buffer, counter);
  }
  slopesMixtureLocal1->fillBufferSlopes(buffer, counter);
  for (int k = 0; k < m_numberTransports; k++) {
    buffer[++counter] = slopesTransportLocal1[k];
  }
  buffer[++counter] = alphaCellAfterOppositeSide;
}

//***********************************************************************

void CellO2::fillBufferSlopesAMR(double *buffer, int &counter, const int &lvl, string whichCpuAmIForNeighbour) const
{
	if (m_lvl == lvl) {
		//Mise a zero de la slope locale a transmettre
		//--------------------------------------------
		double epsilon(1.e-15), sommeCoeff(0.);
		for (int k = 0; k < m_numberPhases; k++) {
			slopesPhasesLocal1[k]->setToZero();
		}
		slopesMixtureLocal1->setToZero();
		for (int k = 0; k < m_numberTransports; k++) {
			slopesTransportLocal1[k] = 0.;
		}

		//Boucle sur les boundaries pour la détermination de la slope a transmettre
		//--------------------------------------------------------------------
    int phase0(0);
    double alphaCellAfterOppositeSide(0.);
		for (unsigned int b = 0; b < m_boundaries.size(); b++) {
      if (!m_boundaries[b]->getSplit()) {
        //Normal of the face in the x-direction
        if (abs(m_boundaries[b]->getFace()->getNormal().getX()) > 0.99) {
          //I am left CPU, I send the average slope of the left side of the right cells
          if (whichCpuAmIForNeighbour == "LEFT") {
            if ((m_boundaries[b]->getFace()->getPos().getX() + epsilon) < m_element->getPosition().getX()) {
              for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
              slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
              for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
              sommeCoeff += 1.;
              if (m_boundaries[b]->getCellGauche() == this) {
                if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
              }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
          }
          //I am right CPU, I send the average slope of the right side of the left cells
          else if (whichCpuAmIForNeighbour == "RIGHT") {
            if ((m_boundaries[b]->getFace()->getPos().getX() - epsilon) > m_element->getPosition().getX()) {
              for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
              slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
              for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
              sommeCoeff += 1.;
              if (m_boundaries[b]->getCellGauche() == this) {
                if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
              }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
          }
        }
        //Normal of the face in the y-direction
        else if (abs(m_boundaries[b]->getFace()->getNormal().getY()) > 0.99) {
          //I am bottom CPU, I send the average slope of the bottom side of the top cells
          if (whichCpuAmIForNeighbour == "BOTTOM") {
            if ((m_boundaries[b]->getFace()->getPos().getY() + epsilon) < m_element->getPosition().getY()) {
              for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
              slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
              for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
              sommeCoeff += 1.;
              if (m_boundaries[b]->getCellGauche() == this) {
                if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
              }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
          }
          //I am top CPU, I send the average slope of the top side of the bottom cells
          else if (whichCpuAmIForNeighbour == "TOP") {
            if ((m_boundaries[b]->getFace()->getPos().getY() - epsilon) > m_element->getPosition().getY()) {
              for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
              slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
              for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
              sommeCoeff += 1.;
              if (m_boundaries[b]->getCellGauche() == this) {
                if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
              }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
          }
        }
        //Normal of the face in the z-direction
        else if (abs(m_boundaries[b]->getFace()->getNormal().getZ()) > 0.99) {
          //I am back CPU, I send the average slope of the back side of the front cells
          if (whichCpuAmIForNeighbour == "BACK") {
            if ((m_boundaries[b]->getFace()->getPos().getZ() + epsilon) < m_element->getPosition().getZ()) {
              for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
              slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
              for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
              sommeCoeff += 1.;
              if (m_boundaries[b]->getCellGauche() == this) {
                if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
              }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
          }
          //I am front CPU, I send the average slope of the front side of the back cells
          else if (whichCpuAmIForNeighbour == "FRONT") {
            if ((m_boundaries[b]->getFace()->getPos().getZ() - epsilon) > m_element->getPosition().getZ()) {
              for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_boundaries[b]->getSlopesPhase(k), 1.); }
              slopesMixtureLocal1->multiplyAndAdd(*m_boundaries[b]->getSlopesMixture(), 1.);
              for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] += (m_boundaries[b]->getSlopesTransport(k)->getValue()); }
              sommeCoeff += 1.;
              if (m_boundaries[b]->getCellGauche() == this) {
                if (m_boundaries[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_boundaries[b]->getCellDroite()->getPhase(phase0)->getAlpha(); }
                else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
              }
              else { alphaCellAfterOppositeSide += m_boundaries[b]->getCellGauche()->getPhase(phase0)->getAlpha(); }
            }
          }
        }
      }
		}

		//Normalisation de la slope
		//-------------------------
		if (sommeCoeff > 1.e-5) {
			for (int k = 0; k < m_numberPhases; k++) { slopesPhasesLocal1[k]->divide(sommeCoeff); }
			slopesMixtureLocal1->divide(sommeCoeff);
			for (int k = 0; k < m_numberTransports; k++) { slopesTransportLocal1[k] /= sommeCoeff; }
      alphaCellAfterOppositeSide /= sommeCoeff;
		}

		//Rempli buffer pour envoi
		//------------------------
		for (int k = 0; k < m_numberPhases; k++) {
			slopesPhasesLocal1[k]->fillBufferSlopes(buffer, counter);
		}
		slopesMixtureLocal1->fillBufferSlopes(buffer, counter);
		for (int k = 0; k < m_numberTransports; k++) {
			buffer[++counter] = slopesTransportLocal1[k];
		}
    buffer[++counter] = alphaCellAfterOppositeSide;
	}

	else {
    if (whichCpuAmIForNeighbour == "LEFT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am left CPU, I send the average slope of the left side of the right cells
        if ((i % 2) == 1) { m_childrenCells[i]->fillBufferSlopesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "RIGHT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am right CPU, I send the average slope of the right side of the left cells
        if ((i % 2) == 0) { m_childrenCells[i]->fillBufferSlopesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BOTTOM") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am bottom CPU, I send the average slope of the bottom side of the top cells
        if ((i % 4) > 1) { m_childrenCells[i]->fillBufferSlopesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "TOP") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am top CPU, I send the average slope of the top side of the bottom cells
        if ((i % 4) <= 1) { m_childrenCells[i]->fillBufferSlopesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BACK") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am back CPU, I send the average slope of the back side of the front cells
        if (i > 3) { m_childrenCells[i]->fillBufferSlopesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "FRONT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am front CPU, I send the average slope of the front side of the back cells
        if (i <= 3) { m_childrenCells[i]->fillBufferSlopesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
	}
}

//***********************************************************************
