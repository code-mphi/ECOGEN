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

#include "BoundCondPeriodicO2Cartesian.h"

//***********************************************************************

BoundCondPeriodicO2Cartesian::BoundCondPeriodicO2Cartesian(const BoundCondPeriodicO2Cartesian& Source, const int& lvl) :
  BoundCondPeriodic(Source, lvl), m_vecPhasesSlopes(0), m_mixtureSlopes(0), m_vecTransportsSlopes(0)
{}

//***********************************************************************

BoundCondPeriodicO2Cartesian::BoundCondPeriodicO2Cartesian(int numPhysique) :
  BoundCondPeriodic(numPhysique), m_vecPhasesSlopes(0), m_mixtureSlopes(0), m_vecTransportsSlopes(0)
{}

//***********************************************************************

BoundCondPeriodicO2Cartesian::~BoundCondPeriodicO2Cartesian()
{
  if (m_vecPhasesSlopes) { // Needed for cases where the vector is not allocated (dealloc of m_limXm for example)
    for (int k = 0; k < numberPhases; k++) {
      delete m_vecPhasesSlopes[k];
    }
    delete[] m_vecPhasesSlopes;
  }
  delete m_mixtureSlopes;
  delete[] m_vecTransportsSlopes;
}

//****************************************************************************

void BoundCondPeriodicO2Cartesian::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondPeriodicO2Cartesian(*(this)));
}

//***********************************************************************

void BoundCondPeriodicO2Cartesian::allocateSlopes(int& /*allocateSlopeLocal*/)
{
  //Allocation of phase slopes
  m_vecPhasesSlopes = new Phase*[numberPhases];
  // We set phases from the left cell (because the right cell doesn't exist for the limits)
  // Necessary because we need to know the phase type (for example: PhasePUEq, etc.)
  // Then all slopes are set to zero
  for (int k = 0; k < numberPhases; k++) {
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
}

//***********************************************************************

void BoundCondPeriodicO2Cartesian::computeSlopes(Prim type)
{
  if (m_cellInterfacesChildren.size() == 0) {
    //Distance between the two cells in contact
    double distance_(this->distance(m_cellLeft) * 2.);

    //Attribution gauche/droite
    //Si la cell gauche ou droite est de niveau inferieur a "lvl", on ne prend pas "type" mais vecPhases (ca evite de prendre vecPhaseO2 alors qu'on ne l'a pas).
    Prim typeGauche = type;
    Prim typeDroite = type;
    if (m_cellLeft->getLvl() < m_lvl) {
      typeGauche = vecPhases;
    }
    if (m_cellRight->getLvl() < m_lvl) {
      typeDroite = vecPhases;
    }

    if (this->getFace()->getNormal().scalar(coordIdentity) > 0.) {
      for (int k = 0; k < numberPhases; k++) {
        m_vecPhasesSlopes[k]->computeSlopesPhase(*m_cellLeft->getPhase(k, typeGauche), *m_cellRight->getPhase(k, typeDroite), distance_);
      }
      m_mixtureSlopes->computeSlopesMixture(*m_cellLeft->getMixture(typeGauche), *m_cellRight->getMixture(typeDroite), distance_);
      for (int k = 0; k < numberTransports; k++) {
        m_vecTransportsSlopes[k].computeSlopeTransport(m_cellLeft->getTransport(k, typeGauche).getValue(),
                                                       m_cellRight->getTransport(k, typeDroite).getValue(),
                                                       distance_);
      }
    }
    else {
      for (int k = 0; k < numberPhases; k++) {
        m_vecPhasesSlopes[k]->computeSlopesPhase(*m_cellRight->getPhase(k, typeDroite), *m_cellLeft->getPhase(k, typeGauche), distance_);
      }
      m_mixtureSlopes->computeSlopesMixture(*m_cellRight->getMixture(typeDroite), *m_cellLeft->getMixture(typeGauche), distance_);
      for (int k = 0; k < numberTransports; k++) {
        m_vecTransportsSlopes[k].computeSlopeTransport(m_cellRight->getTransport(k, typeDroite).getValue(),
                                                       m_cellLeft->getTransport(k, typeGauche).getValue(),
                                                       distance_);
      }
    }
  }
}

//***********************************************************************

void BoundCondPeriodicO2Cartesian::solveRiemann(double& dtMax,
                                                Limiter& globalLimiter,
                                                Limiter& interfaceLimiter,
                                                Limiter& globalVolumeFractionLimiter,
                                                Limiter& interfaceVolumeFractionLimiter,
                                                Prim type)
{
  //Si la cell gauche ou droite est de niveau inferieur a "lvl", on ne prend pas "type" mais vecPhases (ca evite de prendre vecPhaseO2 alors qu'on ne l'a pas).
  if (m_cellLeft->getLvl() == m_lvl) {
    bufferCellLeft->copyVec(m_cellLeft->getPhases(type), m_cellLeft->getMixture(type), m_cellLeft->getTransports(type));
  }
  else {
    bufferCellLeft->copyVec(m_cellLeft->getPhases(vecPhases), m_cellLeft->getMixture(vecPhases), m_cellLeft->getTransports(vecPhases));
  }

  if (m_cellRight->getLvl() == m_lvl) {
    bufferCellRight->copyVec(m_cellRight->getPhases(type), m_cellRight->getMixture(type), m_cellRight->getTransports(type));
  }
  else {
    bufferCellRight->copyVec(m_cellRight->getPhases(vecPhases), m_cellRight->getMixture(vecPhases), m_cellRight->getTransports(vecPhases));
  }

  //Calcul des distances cell interface <-> cells pour l extrapolation
  double distanceGauche(this->distance(m_cellLeft));
  double distanceDroite(distanceGauche);

  //Initialization of variables to detect the interface for limiters and THINC method
  int phase0(0), phase1(1);
  double alphaCellLeft(0.), alphaCellLeftLeft(0.), alphaCellRight(0.), alphaCellRightRight(0.);
  double beta(1.6), sign(0.), newAlpha(0.), A(0.), B(0.), C(0.), qmin(0.), qmax(0.), epsInterface(1.e-4);
  alphaCellLeft  = bufferCellLeft->getPhase(phase0)->getAlpha();
  alphaCellRight = bufferCellRight->getPhase(phase0)->getAlpha();

  //Extrapolation gauche
  m_cellLeft->computeLocalSlopesPeriodic(*this,
                                         globalLimiter,
                                         interfaceLimiter,
                                         globalVolumeFractionLimiter,
                                         interfaceVolumeFractionLimiter,
                                         alphaCellLeftLeft,
                                         alphaCellLeft,
                                         alphaCellRight,
                                         epsInterface);
  for (int k = 0; k < numberPhases; k++) {
    //Face along X
    if (std::fabs(m_face->getNormal().getX()) > 0.5) {
      if (m_cellLeft->getPosition().getX() < m_cellRight->getPosition().getX()) slopesPhasesLocal1[k]->changeSign();
    }
    //Face along Y
    else if (std::fabs(m_face->getNormal().getY()) > 0.5) {
      if (m_cellLeft->getPosition().getY() < m_cellRight->getPosition().getY()) slopesPhasesLocal1[k]->changeSign();
    }
    //Face along Z
    else {
      if (m_cellLeft->getPosition().getZ() < m_cellRight->getPosition().getZ()) slopesPhasesLocal1[k]->changeSign();
    }
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
    if ((alphaCellLeft >= epsInterface) && (alphaCellLeft <= 1. - epsInterface) &&
        ((alphaCellRight - alphaCellLeft) * (alphaCellLeft - alphaCellLeftLeft) >= 1.e-8)) {
      if (alphaCellRight - alphaCellLeftLeft > 0.) {
        sign = 1.;
      }
      else {
        sign = -1.;
      }
      qmin     = std::min(alphaCellRight, alphaCellLeftLeft);
      qmax     = std::max(alphaCellRight, alphaCellLeftLeft) - qmin;
      C        = (alphaCellLeft - qmin + 1.e-20) / (qmax + 1.e-20);
      B        = exp(sign * beta * (2. * C - 1.));
      A        = (B / cosh(beta) - 1.) / tanh(beta);
      newAlpha = qmin + 0.5 * qmax * (1. + sign * (tanh(beta) + A) / (1. + A * tanh(beta)));
      if (newAlpha < epsInterface) {
        newAlpha = epsInterface;
      }
      if (newAlpha > 1. - epsInterface) {
        newAlpha = 1. - epsInterface;
      }
      bufferCellLeft->getPhase(phase0)->setAlpha(newAlpha);
      bufferCellLeft->getPhase(phase1)->setAlpha(1. - newAlpha);
    }
  }

  //Extrapolation droite
  m_cellRight->computeLocalSlopesPeriodic(*this,
                                          globalLimiter,
                                          interfaceLimiter,
                                          globalVolumeFractionLimiter,
                                          interfaceVolumeFractionLimiter,
                                          alphaCellRightRight,
                                          alphaCellRight,
                                          alphaCellLeft,
                                          epsInterface);
  for (int k = 0; k < numberPhases; k++) {
    //---------------------------------------
    //Face along X
    if (std::fabs(m_face->getNormal().getX()) > 0.5) {
      if (m_cellLeft->getPosition().getX() > m_cellRight->getPosition().getX()) slopesPhasesLocal1[k]->changeSign();
    }
    //Face along Y
    else if (std::fabs(m_face->getNormal().getY()) > 0.5) {
      if (m_cellLeft->getPosition().getY() > m_cellRight->getPosition().getY()) slopesPhasesLocal1[k]->changeSign();
    }
    //Face along Z
    else {
      if (m_cellLeft->getPosition().getZ() > m_cellRight->getPosition().getZ()) slopesPhasesLocal1[k]->changeSign();
    }
    bufferCellRight->getPhase(k)->extrapolate(*slopesPhasesLocal1[k], distanceDroite);
    bufferCellRight->getPhase(k)->verifyAndCorrectPhase();
    bufferCellRight->getPhase(k)->verifyAndCorrectDensityMax();
  }
  // slopesMixtureLocal1->changeSign();
  bufferCellRight->getMixture()->extrapolate(*slopesMixtureLocal1, distanceDroite);
  for (int k = 0; k < numberTransports; k++) {
    // slopesTransportLocal1[k] = -slopesTransportLocal1[k];
    bufferCellRight->getTransport(k).extrapolate(slopesTransportLocal1[k], distanceDroite);
  }
  //THINC method (for alpha only)
  if (globalVolumeFractionLimiter.AmITHINC() || interfaceVolumeFractionLimiter.AmITHINC()) {
    if ((alphaCellRight >= epsInterface) && (alphaCellRight <= 1. - epsInterface) &&
        ((alphaCellRightRight - alphaCellRight) * (alphaCellRight - alphaCellLeft) >= 1.e-8)) {
      if (alphaCellRightRight - alphaCellLeft > 0.) {
        sign = 1.;
      }
      else {
        sign = -1.;
      }
      qmin     = std::min(alphaCellRightRight, alphaCellLeft);
      qmax     = std::max(alphaCellRightRight, alphaCellLeft) - qmin;
      C        = (alphaCellRight - qmin + 1.e-20) / (qmax + 1.e-20);
      B        = exp(sign * beta * (2. * C - 1.));
      A        = (B / cosh(beta) - 1.) / tanh(beta);
      newAlpha = qmin + 0.5 * qmax * (1. + sign * A);
      if (newAlpha < epsInterface) {
        newAlpha = epsInterface;
      }
      if (newAlpha > 1. - epsInterface) {
        newAlpha = 1. - epsInterface;
      }
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

  // Probleme de Riemann
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  double dxRight(m_cellRight->getElement()->getLCFL());
  dxLeft  = dxLeft * std::pow(2., (double)m_lvl);
  dxRight = dxRight * std::pow(2., (double)m_lvl);
  model->solveRiemannIntern(*bufferCellLeft, *bufferCellRight, dxLeft, dxRight, dtMax);
  //Handling of transport functions (m_Sm known: need to be called after Riemann solver)
  if (numberTransports > 0) {
    model->solveRiemannTransportIntern(*bufferCellLeft, *bufferCellRight);
  }

  //Projection du flux sur le repere absolu
  model->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
}

//***********************************************************************

Phase* BoundCondPeriodicO2Cartesian::getSlopesPhase(const int& phaseNumber) const { return m_vecPhasesSlopes[phaseNumber]; }

//***********************************************************************

Mixture* BoundCondPeriodicO2Cartesian::getSlopesMixture() const { return m_mixtureSlopes; }

//***********************************************************************

Transport* BoundCondPeriodicO2Cartesian::getSlopesTransport(const int& numberTransport) const { return &m_vecTransportsSlopes[numberTransport]; }

//***************************************************************************
//******************************AMR Method***********************************
//***************************************************************************

void BoundCondPeriodicO2Cartesian::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondPeriodicO2Cartesian(*this, m_lvl + 1));
}

//****************************************************************************
