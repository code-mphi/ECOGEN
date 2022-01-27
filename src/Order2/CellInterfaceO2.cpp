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

CellInterfaceO2::CellInterfaceO2() : CellInterface(), m_vecPhasesSlopes(0), m_mixtureSlopes(0), m_vecTransportsSlopes(0)
 //m_BG1M(0), m_BG2M(0), m_BG3M(0), m_BG1P(0), m_BG2P(0), m_BG3P(0), m_BD1M(0),
 //m_BD2M(0), m_BD3M(0), m_BD1P(0), m_BD2P(0), m_BD3P(0),
 //m_betaG1M(1.), m_betaG2M(0.), m_betaG3M(0.), m_betaG1P(1.), m_betaG2P(0.), m_betaG3P(0.),
 //m_betaD1M(1.), m_betaD2M(0.), m_betaD3M(0.), m_betaD1P(1.), m_betaD2P(0.), m_betaD3P(0.),
 //m_distanceHGM(0.), m_distanceHGP(0.), m_distanceHDM(0.), m_distanceHDP(0.)
 {}

//***********************************************************************

CellInterfaceO2::CellInterfaceO2(int lvl) : CellInterface(lvl), m_vecPhasesSlopes(0), m_mixtureSlopes(0), m_vecTransportsSlopes(0)
//m_BG1M(0), m_BG2M(0), m_BG3M(0), m_BG1P(0), m_BG2P(0), m_BG3P(0), m_BD1M(0),
//m_BD2M(0), m_BD3M(0), m_BD1P(0), m_BD2P(0), m_BD3P(0),
//m_betaG1M(1.), m_betaG2M(0.), m_betaG3M(0.), m_betaG1P(1.), m_betaG2P(0.), m_betaG3P(0.),
//m_betaD1M(1.), m_betaD2M(0.), m_betaD3M(0.), m_betaD1P(1.), m_betaD2P(0.), m_betaD3P(0.),
//m_distanceHGM(0.), m_distanceHGP(0.), m_distanceHDM(0.), m_distanceHDP(0.)
{}

//***********************************************************************

CellInterfaceO2::~CellInterfaceO2()
{
  for (int k = 0; k < numberPhases; k++) {
    delete m_vecPhasesSlopes[k];
  }
  delete[] m_vecPhasesSlopes;
  delete m_mixtureSlopes;
  delete[] m_vecTransportsSlopes;
}

//***********************************************************************

void CellInterfaceO2::allocateSlopes(int& allocateSlopeLocal)
{
  //Allocation des slopes des phases
	m_vecPhasesSlopes = new Phase*[numberPhases];

  //On attribut les phases a partir de la cell a gauche (car cell a droite inexistante pour les limites)
  //Necessaire car il faut connaitre le type de phase (ex: PhasePUEq, etc.))
  //Ensuite on met a zero toutes les slopes
  for(int k = 0; k < numberPhases; k++){
    m_cellLeft->getPhase(k)->allocateAndCopyPhase(&m_vecPhasesSlopes[k]);
		m_vecPhasesSlopes[k]->setToZero();
  }
  m_cellLeft->getMixture()->allocateAndCopyMixture(&m_mixtureSlopes);
  m_mixtureSlopes->setToZero();

	//Allocation des slopes sur transports
	m_vecTransportsSlopes = new Transport[numberTransports];
	for (int k = 0; k < numberTransports; k++) {
		m_vecTransportsSlopes[k].setValue(0.);
	}
  
  //Allocation des variables externes
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

void CellInterfaceO2::computeSlopes(Prim type)
{
  if (m_cellInterfacesChildren.size() == 0) {
    //Distance entre les deux mailles en contact
    double distance(m_cellLeft->distance(m_cellRight));
    //Attribution gauche/droite
    //Si la cell gauche ou droite est de niveau inferieur a "lvl", on ne prend pas "type" mais vecPhases (ca evite de prendre vecPhaseO2 alors qu'on ne l'a pas).
    Prim typeGauche = type;
    Prim typeDroite = type;
    if (m_cellLeft->getLvl() < m_lvl) { typeGauche = vecPhases; }
    if (m_cellRight->getLvl() < m_lvl) { typeDroite = vecPhases; }

    for (int k = 0; k < numberPhases; k++) {
			m_vecPhasesSlopes[k]->computeSlopesPhase(*m_cellLeft->getPhase(k, typeGauche), *m_cellRight->getPhase(k, typeDroite), distance);
    }
    m_mixtureSlopes->computeSlopesMixture(*m_cellLeft->getMixture(typeGauche), *m_cellRight->getMixture(typeDroite), distance);
    for (int k = 0; k < numberTransports; k++) {
      m_vecTransportsSlopes[k].computeSlopeTransport(m_cellLeft->getTransport(k, typeGauche).getValue(), m_cellRight->getTransport(k, typeDroite).getValue(), distance);
    }
  }
}

//***********************************************************************

void CellInterfaceO2::computeFlux(double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter,
 Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type)
{
  // Quand on fait le premier computeFlux (donc avec vecPhases) on n'incremente pas m_cons pour les mailles de niveau different (inferieur) de "lvl".
  // Sinon ca veut dire qu on l ajoute pour les 2 computeFlux sans le remettre a zero entre les deux, donc 2 fois plus de flux que ce que l on veut.
  this->solveRiemann(dtMax, globalLimiter, interfaceLimiter, globalVolumeFractionLimiter, interfaceVolumeFractionLimiter, type);

  switch (type) {
  case vecPhases:
    if (m_cellLeft->getLvl() == m_cellRight->getLvl()) {       //CoefAMR = 1 pour les deux
      this->addFlux(1.);       //Ajout du flux sur maille droite
      this->subtractFlux(1.);  //Retrait du flux sur maille gauche
    }
    else if (m_cellLeft->getLvl() > m_cellRight->getLvl()) {   //CoefAMR = 1 pour la gauche et on n'ajoute rien sur la maille droite
      this->subtractFlux(1.);  //Retrait du flux sur maille gauche
    }
    else {                                                     //CoefAMR = 1 pour la droite et on ne retire rien sur la maille gauche
      this->addFlux(1.);       //Ajout du flux sur maille droite
    }
    break;

  case vecPhasesO2:
    if (m_cellLeft->getLvl() == m_cellRight->getLvl()) {       //CoefAMR = 1 pour les deux
      this->addFlux(1.);       //Ajout du flux sur maille droite
      this->subtractFlux(1.);  //Retrait du flux sur maille gauche
    }
    else if (m_cellLeft->getLvl() > m_cellRight->getLvl()) {   //CoefAMR = 1 pour la gauche et 0.5 pour la droite
      this->addFlux(0.5);      //Ajout du flux sur maille droite
      this->subtractFlux(1.);  //Retrait du flux sur maille gauche
    }
    else {                                                     //CoefAMR = 1 pour la droite et 0.5 pour la gauche
      this->addFlux(1.);       //Ajout du flux sur maille droite
      this->subtractFlux(0.5); //Retrait du flux sur maille gauche
    }
    break;

  default: break;
  }
}

//***********************************************************************

void CellInterfaceO2::solveRiemann(double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type)
{
  //Si la cell gauche ou droite est de niveau inferieur a "lvl", on ne prend pas "type" mais vecPhases (ca evite de prendre vecPhaseO2 alors qu'on ne l'a pas).
  if (m_cellLeft->getLvl() == m_lvl) { bufferCellLeft->copyVec(m_cellLeft->getPhases(type), m_cellLeft->getMixture(type), m_cellLeft->getTransports(type)); }
  else { bufferCellLeft->copyVec(m_cellLeft->getPhases(vecPhases), m_cellLeft->getMixture(vecPhases), m_cellLeft->getTransports(vecPhases)); }
  
  if (m_cellRight->getLvl() == m_lvl) { bufferCellRight->copyVec(m_cellRight->getPhases(type), m_cellRight->getMixture(type), m_cellRight->getTransports(type)); }
  else { bufferCellRight->copyVec(m_cellRight->getPhases(vecPhases), m_cellRight->getMixture(vecPhases), m_cellRight->getTransports(vecPhases)); }

  //Calcul des distances cell interface <-> cells pour l extrapolation
  double distanceGauche(this->distance(m_cellLeft));
  double distanceDroite(this->distance(m_cellRight));

  //Initialization of variables to detect the interface for limiters and THINC method
  int phase0(0), phase1(1);
  double alphaCellLeft(0.), alphaCellLeftLeft(0.), alphaCellRight(0.), alphaCellRightRight(0.);
  double beta(1.6), sign(0.), newAlpha(0.), A(0.), B(0.), C(0.), qmin(0.), qmax(0.), epsInterface(1.e-4);
  alphaCellLeft = bufferCellLeft->getPhase(phase0)->getAlpha();
  alphaCellRight = bufferCellRight->getPhase(phase0)->getAlpha();

  //Extrapolation gauche
  m_cellLeft->computeLocalSlopes(*this, globalLimiter, interfaceLimiter, globalVolumeFractionLimiter, interfaceVolumeFractionLimiter, alphaCellLeftLeft, alphaCellLeft, alphaCellRight, epsInterface);
  for (int k = 0; k < numberPhases; k++) {
    bufferCellLeft->getPhase(k)->extrapolate(*slopesPhasesLocal1[k], distanceGauche);
    bufferCellLeft->getPhase(k)->verifyAndCorrectPhase();
    bufferCellLeft->getPhase(k)->verifyAndCorrectDensityMax();
  }
  bufferCellLeft->getMixture()->extrapolate(*slopesMixtureLocal1, distanceGauche);
	for (int k = 0; k < numberTransports; k++) {
		bufferCellLeft->getTransport(k).extrapolate(slopesTransportLocal1[k], distanceGauche);
	}
  //THINC method (for alpha only)
  if (globalVolumeFractionLimiter.AmITHINC() || interfaceVolumeFractionLimiter.AmITHINC()) {
    if ((alphaCellLeft >= epsInterface) && (alphaCellLeft <= 1. - epsInterface) && ((alphaCellRight - alphaCellLeft)*(alphaCellLeft - alphaCellLeftLeft) >= 1.e-8)) {
      if (alphaCellRight - alphaCellLeftLeft > 0.) { sign = 1.; }
      else { sign = -1.; }
      qmin = std::min(alphaCellRight, alphaCellLeftLeft);
      qmax = std::max(alphaCellRight, alphaCellLeftLeft) - qmin;
      C = (alphaCellLeft - qmin + 1.e-20) / (qmax + 1.e-20);
      B = exp(sign*beta*(2.*C - 1.));
      A = (B / cosh(beta) - 1.) / tanh(beta);
      newAlpha = qmin + 0.5*qmax*(1. + sign * (tanh(beta) + A) / (1. + A * tanh(beta)));
      if (newAlpha < epsInterface) { newAlpha = epsInterface; }
      if (newAlpha > 1. - epsInterface) { newAlpha = 1. - epsInterface; }
      bufferCellLeft->getPhase(phase0)->setAlpha(newAlpha);
      bufferCellLeft->getPhase(phase1)->setAlpha(1. - newAlpha);
    }
  }

  //Extrapolation droite
  m_cellRight->computeLocalSlopes(*this, globalLimiter, interfaceLimiter, globalVolumeFractionLimiter, interfaceVolumeFractionLimiter, alphaCellRightRight, alphaCellRight, alphaCellLeft, epsInterface);
  for (int k = 0; k < numberPhases; k++) {
    slopesPhasesLocal1[k]->changeSign(); //On doit soustraire les slopes a droite
    bufferCellRight->getPhase(k)->extrapolate(*slopesPhasesLocal1[k], distanceDroite);
    bufferCellRight->getPhase(k)->verifyAndCorrectPhase();
    bufferCellRight->getPhase(k)->verifyAndCorrectDensityMax();
  }
  slopesMixtureLocal1->changeSign();
  bufferCellRight->getMixture()->extrapolate(*slopesMixtureLocal1, distanceDroite);
	for (int k = 0; k < numberTransports; k++) {
		slopesTransportLocal1[k] = -slopesTransportLocal1[k];
		bufferCellRight->getTransport(k).extrapolate(slopesTransportLocal1[k], distanceDroite);
	}
  //THINC method (for alpha only)
  if (globalVolumeFractionLimiter.AmITHINC() || interfaceVolumeFractionLimiter.AmITHINC()) {
    if ((alphaCellRight >= epsInterface) && (alphaCellRight <= 1. - epsInterface) && ((alphaCellRightRight - alphaCellRight)*(alphaCellRight - alphaCellLeft) >= 1.e-8)) {
      if (alphaCellRightRight - alphaCellLeft > 0.) { sign = 1.; }
      else { sign = -1.; }
      qmin = std::min(alphaCellRightRight, alphaCellLeft);
      qmax = std::max(alphaCellRightRight, alphaCellLeft) - qmin;
      C = (alphaCellRight - qmin + 1.e-20) / (qmax + 1.e-20);
      B = exp(sign*beta*(2.*C - 1.));
      A = (B / cosh(beta) - 1.) / tanh(beta);
      newAlpha = qmin + 0.5*qmax*(1. + sign * A);
      if (newAlpha < epsInterface) { newAlpha = epsInterface; }
      if (newAlpha > 1. - epsInterface) { newAlpha = 1. - epsInterface; }
      bufferCellRight->getPhase(phase0)->setAlpha(newAlpha);
      bufferCellRight->getPhase(phase1)->setAlpha(1. - newAlpha);
    }
  }

  //Projection des velocities sur repere attache a la face
  bufferCellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
  bufferCellRight->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());

  //Calcul des variables etendus (Phases, Mixture, AddPhys)
  bufferCellLeft->fulfillState();
  bufferCellRight->fulfillState();

  //Probleme de Riemann
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  double dxRight(m_cellRight->getElement()->getLCFL());
  dxLeft = dxLeft*std::pow(2., (double)m_lvl);
  dxRight = dxRight*std::pow(2., (double)m_lvl);
  model->solveRiemannIntern(*bufferCellLeft, *bufferCellRight, dxLeft, dxRight, dtMax);
  //Handling of transport functions (m_Sm known: need to be called after Riemann solver)
  if (numberTransports > 0) { model->solveRiemannTransportIntern(*bufferCellLeft, *bufferCellRight); }

  //Projection du flux sur le repere absolu
  model->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
}

//***********************************************************************

Phase* CellInterfaceO2::getSlopesPhase(const int& phaseNumber) const
{
  return m_vecPhasesSlopes[phaseNumber];
}

//***********************************************************************

Mixture* CellInterfaceO2::getSlopesMixture() const
{
  return m_mixtureSlopes;
}

//***********************************************************************

Transport* CellInterfaceO2::getSlopesTransport(const int& numberTransport) const
{
  return &m_vecTransportsSlopes[numberTransport];
}

//***********************************************************************

//Cell*  CellInterfaceO2::getB(BO2 B) const 
//{
//  switch (B){
//  case BG1M: return m_BG1M; break;
//  case BG2M: return m_BG2M; break;
//  case BG1P: return m_BG1P; break;
//  case BG2P: return m_BG2P; break;
//  case BD1M: return m_BD1M; break;
//  case BD2M: return m_BD2M; break;
//  case BD1P: return m_BD1P; break;
//  case BD2P: return m_BD2P; break;
//  default: Errors::errorMessage("probleme enum non connu dans BorDeMailleO2::getB"); return 0; break;
//  }
//}

//***********************************************************************

//double CellInterfaceO2::getBeta(betaO2 beta) const 
//{
//  switch (beta) {
//  case betaG1M: return m_betaG1M; break;
//  case betaG2M: return m_betaG2M; break;
//  case betaG1P: return m_betaG1P; break;
//  case betaG2P: return m_betaG2P; break;
//  case betaD1M: return m_betaD1M; break;
//  case betaD2M: return m_betaD2M; break;
//  case betaD1P: return m_betaD1P; break;
//  case betaD2P: return m_betaD2P; break;
//  default: Errors::errorMessage("probleme enum non connu dans BorDeMailleO2::getBeta"); return 0; break;
//  }
//}

//***********************************************************************

//double CellInterfaceO2::getDistanceH(distanceHO2 dist) const
//{
//  switch (dist) {
//  case distanceHGM: return m_distanceHGM; break;
//  case distanceHGP: return m_distanceHGP; break;
//  case distanceHDM: return m_distanceHDM; break;
//  case distanceHDP: return m_distanceHDP; break;
//  default: Errors::errorMessage("probleme enum non connu dans BorDeMailleO2::getDistanceH"); return 0; break;
//  }
//}

//***********************************************************************

//void CellInterfaceO2::setB(BO2 B, Cell* cell)
//{
//  switch (B) {
//  case BG1M: m_BG1M = cell; break;
//  case BG2M: m_BG2M = cell; break;
//  case BG1P: m_BG1P = cell; break;
//  case BG2P: m_BG2P = cell; break;
//  case BD1M: m_BD1M = cell; break;
//  case BD2M: m_BD2M = cell; break;
//  case BD1P: m_BD1P = cell; break;
//  case BD2P: m_BD2P = cell; break;
//  default: Errors::errorMessage("probleme enum non connu dans BorDeMailleO2::setB"); break;
//  }
//}

//***********************************************************************

//void CellInterfaceO2::setBeta(betaO2 beta, double& value)
//{
//  switch (beta) {
//  case betaG1M: m_betaG1M = value; break;
//  case betaG2M: m_betaG2M = value; break;
//  case betaG1P: m_betaG1P = value; break;
//  case betaG2P: m_betaG2P = value; break;
//  case betaD1M: m_betaD1M = value; break;
//  case betaD2M: m_betaD2M = value; break;
//  case betaD1P: m_betaD1P = value; break;
//  case betaD2P: m_betaD2P = value; break;
//  default: Errors::errorMessage("probleme enum non connu dans BorDeMailleO2::setB"); break;
//  }
//}

//***********************************************************************

//void CellInterfaceO2::setDistanceH(distanceHO2 dist, double& value) 
//{
//  switch (dist) {
//  case distanceHGM: m_distanceHGM = value; break;
//  case distanceHGP: m_distanceHGP = value; break;
//  case distanceHDM: m_distanceHDM = value; break;
//  case distanceHDP: m_distanceHDP = value; break;
//  default: Errors::errorMessage("probleme enum non connu dans BorDeMailleO2::setDistanceH"); break;
//  }
//}

//****************************************************************************
//******************************AMR Method***********************************
//****************************************************************************

void CellInterfaceO2::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new CellInterfaceO2(m_lvl + 1));
}

//***********************************************************************

void CellInterfaceO2::creerCellInterfaceChildInterne(const int& lvl, std::vector<CellInterface*>* childrenInternalCellInterfaces)
{
  (*childrenInternalCellInterfaces).push_back(new CellInterfaceO2(lvl + 1));
}

//***********************************************************************