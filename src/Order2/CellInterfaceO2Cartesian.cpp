#include "CellInterfaceO2Cartesian.h"

//***********************************************************************

CellInterfaceO2Cartesian::CellInterfaceO2Cartesian() : CellInterfaceO2(), m_vecPhasesSlopes(0), m_mixtureSlopes(0), m_vecTransportsSlopes(0)
{}

//***********************************************************************

CellInterfaceO2Cartesian::CellInterfaceO2Cartesian(int lvl) : CellInterfaceO2(lvl), m_vecPhasesSlopes(0), m_mixtureSlopes(0), m_vecTransportsSlopes(0)
{}

//***********************************************************************

CellInterfaceO2Cartesian::~CellInterfaceO2Cartesian()
{
  for (int k = 0; k < numberPhases; k++) {
    delete m_vecPhasesSlopes[k];
  }
  delete[] m_vecPhasesSlopes;
  delete m_mixtureSlopes;
  delete[] m_vecTransportsSlopes;
}

//***********************************************************************

void CellInterfaceO2Cartesian::allocateSlopes(int& allocateSlopeLocal)
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

void CellInterfaceO2Cartesian::computeSlopes(Prim type)
{
  if (m_cellInterfacesChildren.size() == 0) {
    //Distance between the two cells in contact
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

void CellInterfaceO2Cartesian::solveRiemann(double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type)
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

  //Computation of extended variables (Phases, Mixture, AddPhys)
  bufferCellLeft->fulfillState();
  bufferCellRight->fulfillState();

  //Vector and tensor projections on geometric reference frame of the face
  bufferCellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
  bufferCellRight->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());

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

Phase* CellInterfaceO2Cartesian::getSlopesPhase(const int& phaseNumber) const
{
  return m_vecPhasesSlopes[phaseNumber];
}

//***********************************************************************

Mixture* CellInterfaceO2Cartesian::getSlopesMixture() const
{
  return m_mixtureSlopes;
}

//***********************************************************************

Transport* CellInterfaceO2Cartesian::getSlopesTransport(const int& numberTransport) const
{
  return &m_vecTransportsSlopes[numberTransport];
}

//****************************************************************************
//******************************AMR Method***********************************
//****************************************************************************

void CellInterfaceO2Cartesian::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new CellInterfaceO2Cartesian(m_lvl + 1));
}

//***********************************************************************

void CellInterfaceO2Cartesian::creerCellInterfaceChildInterne(const int& lvl, std::vector<CellInterface*>* childrenInternalCellInterfaces)
{
  (*childrenInternalCellInterfaces).push_back(new CellInterfaceO2Cartesian(lvl + 1));
}

//***********************************************************************