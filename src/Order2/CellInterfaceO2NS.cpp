#include "CellInterfaceO2NS.h"

//***********************************************************************

CellInterfaceO2NS::CellInterfaceO2NS()
{
}

//***********************************************************************

CellInterfaceO2NS::~CellInterfaceO2NS()
{
}

//***********************************************************************

void CellInterfaceO2NS::solveRiemann(double& dtMax, Limiter& /*globalLimiter*/, Limiter& /*interfaceLimiter*/, Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/, Prim type)
{
  // Initialize buffer cells used in Riemann problem
  bufferCellLeft->copyVec(m_cellLeft->getPhases(type), m_cellLeft->getMixture(type), m_cellLeft->getTransports(type));
  bufferCellRight->copyVec(m_cellRight->getPhases(type), m_cellRight->getMixture(type), m_cellRight->getTransports(type));

  // For NS extrapolation a scalar product between the distance
  // Coord rij and the gradient is already done during CellO2NS::computeLocalSlopes
  // and stored in slopesPhase/MixtureLocal1
  // Therefore, to stay compliant with cartesian 2nd order, distance is equal to 1 here
  double distanceLeft = 1.0; 
  double distanceRight = 1.0;

  // Left side of interface extrapolation
  m_cellLeft->computeLocalSlopes(*this); // Build slopesPhasesLocal1 = theta_i * rij . grad(Wi)
  for (int k = 0; k < numberPhases; k++) {
    bufferCellLeft->getPhase(k)->extrapolate(*slopesPhasesLocal1[k], distanceLeft); // Build Wij_lim = Wi + theta_i * rij . grad(Wi)
    bufferCellLeft->getPhase(k)->verifyAndCorrectPhase();
    bufferCellLeft->getPhase(k)->verifyAndCorrectDensityMax();
  }
  bufferCellLeft->getMixture()->extrapolate(*slopesMixtureLocal1, distanceLeft);
	for (int k = 0; k < numberTransports; k++) {
		bufferCellLeft->getTransport(k).extrapolate(slopesTransportLocal1[k], distanceLeft);
	}

  // Right side of interface extrapolation
  m_cellRight->computeLocalSlopes(*this);
  for (int k = 0; k < numberPhases; k++) {
    bufferCellRight->getPhase(k)->extrapolate(*slopesPhasesLocal1[k], distanceRight);
    bufferCellRight->getPhase(k)->verifyAndCorrectPhase();
    bufferCellRight->getPhase(k)->verifyAndCorrectDensityMax();
  }
  bufferCellRight->getMixture()->extrapolate(*slopesMixtureLocal1, distanceRight);
	for (int k = 0; k < numberTransports; k++) {
		bufferCellRight->getTransport(k).extrapolate(slopesTransportLocal1[k], distanceRight);
	}

  // Compute extended variables since projection change some of them (Phases, Mixture, AddPhys)
  bufferCellLeft->fulfillState();
  bufferCellRight->fulfillState();

  // Vector and tensor projections into reference frame attached to the face
  bufferCellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
  bufferCellRight->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());

  // Riemann problem
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  double dxRight(m_cellRight->getElement()->getLCFL());
  model->solveRiemannIntern(*bufferCellLeft, *bufferCellRight, dxLeft, dxRight, dtMax);
  // Handling of transport functions (m_Sm known: need to be called after Riemann solver)
  if (numberTransports > 0) { model->solveRiemannTransportIntern(*bufferCellLeft, *bufferCellRight); }

  // Flux projection into the global reference frame
  model->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
}

//***********************************************************************