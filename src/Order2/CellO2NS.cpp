#include "CellO2NS.h"

Phase** buffPhasesMin;
Phase** buffPhasesMax;
Mixture* buffMixtureMin;
Mixture* buffMixtureMax;
double* buffTransportMin;
double* buffTransportMax;

//***********************************************************************

CellO2NS::CellO2NS() : CellO2(), m_gradPhase(nullptr), m_gradMixture(nullptr), m_gradTransport(nullptr)
{
}

//***********************************************************************

CellO2NS::~CellO2NS()
{
  for (int k = 0; k < numberPhases; k++) {
    delete m_gradPhase[k];
  }
  delete[] m_gradPhase;
  delete m_gradMixture;
  if (m_gradTransport != nullptr) {
    delete[] m_gradTransport;
  }
}

//***********************************************************************

void CellO2NS::allocateSecondOrderBuffersAndGradientVectors(Phase** phases, Mixture* mixture)
{ 
  buffPhasesMin = new Phase*[numberPhases];
  buffPhasesMax = new Phase*[numberPhases];
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->allocateAndCopyPhase(&buffPhasesMin[k]);
    phases[k]->allocateAndCopyPhase(&buffPhasesMax[k]);
  }
  mixture->allocateAndCopyMixture(&buffMixtureMin);
  mixture->allocateAndCopyMixture(&buffMixtureMax);

  if (numberTransports > 0) {
    buffTransportMin = new double[numberTransports];
    buffTransportMax = new double[numberTransports];
    for (int t = 0; t < numberTransports; t++) {
      buffTransportMin[t] = 0.;
      buffTransportMax[t] = 0.;
    }
  }

  m_gradPhase[0]->initializeGradientVectors();
  m_gradMixture->initializeGradientVectors();
  if (numberTransports > 0) m_gradTransport[0].initializeGradientVectors();
}

//***********************************************************************

void CellO2NS::allocate(const std::vector<AddPhys*>& addPhys)
{
  CellO2::allocate(addPhys);

  m_gradPhase = new GradPhase*[numberPhases];
  for (int k = 0; k < numberSolids; k++) {
    model->allocatePhaseSolidGradient(&m_gradPhase[k]);
  }
  for (int k = numberSolids; k < numberPhases; k++) {
    model->allocatePhaseGradient(&m_gradPhase[k]);
  }
  model->allocateMixtureGradient(&m_gradMixture);

  if (numberTransports > 0) {
    m_gradTransport = new GradTransport[numberTransports];
  }
}

//***********************************************************************

void CellO2NS::computeGradientsO2() 
{
  for (int k = 0; k < numberPhases; k++) {
    m_gradPhase[k]->computeGradients(this, k);
  }
  m_gradMixture->computeGradients(this);
  
  for (int t = 0; t < numberTransports; t++) {
    m_gradTransport[t].computeGradient(this, t);
  }
}

//***********************************************************************

void CellO2NS::limitGradientsO2(Limiter& globalLimiter)
{
  // Reset local slopes
	for (int k = 0; k < numberPhases; k++) {
		slopesPhasesLocal1[k]->setToZero();
		slopesPhasesLocal2[k]->setToZero();
	}
	slopesMixtureLocal1->setToZero();
	slopesMixtureLocal2->setToZero();
	for (int t = 0; t < numberTransports; t++) {
		slopesTransportLocal1[t] = 0.;
		slopesTransportLocal2[t] = 0.;
	}

  // 1) Determine Wmin/Wmax
  buffMixtureMin->copyMixture(*(this->getMixture()));
  buffMixtureMax->copyMixture(*(this->getMixture()));
  for (int k = 0; k < numberPhases; k++) {
    buffPhasesMin[k]->copyPhase(*(this->getPhase(k)));
    buffPhasesMax[k]->copyPhase(*(this->getPhase(k)));
  }
  for (int t = 0; t < numberTransports; t++) {
    buffTransportMin[t] = this->getTransport(t).getValue();
    buffTransportMax[t] = this->getTransport(t).getValue();
  }

  for (unsigned int i = 0; i < m_cellInterfaces.size(); i++) {
    // Inner domain cell interface (not a boundary)
    if (m_cellInterfaces[i]->whoAmI() == 0) {
      for (int k = 0; k < numberPhases; k++) {
        buffPhasesMin[k]->setMin(*(m_cellInterfaces[i]->getCellLeft()->getPhase(k)), *buffPhasesMin[k]);
        buffPhasesMin[k]->setMin(*(m_cellInterfaces[i]->getCellRight()->getPhase(k)), *buffPhasesMin[k]);
      
        buffPhasesMax[k]->setMax(*(m_cellInterfaces[i]->getCellLeft()->getPhase(k)), *buffPhasesMax[k]);
        buffPhasesMax[k]->setMax(*(m_cellInterfaces[i]->getCellRight()->getPhase(k)), *buffPhasesMax[k]);
      }
      buffMixtureMin->setMin(*(m_cellInterfaces[i]->getCellLeft()->getMixture()), *buffMixtureMin);
      buffMixtureMin->setMin(*(m_cellInterfaces[i]->getCellRight()->getMixture()), *buffMixtureMin);

      buffMixtureMax->setMax(*(m_cellInterfaces[i]->getCellLeft()->getMixture()), *buffMixtureMax);
      buffMixtureMax->setMax(*(m_cellInterfaces[i]->getCellRight()->getMixture()), *buffMixtureMax);

      for (int t = 0; t < numberTransports; t++) {
        buffTransportMin[t] = std::min(m_cellInterfaces[i]->getCellLeft()->getTransport(t).getValue(), buffTransportMin[t]);
        buffTransportMin[t] = std::min(m_cellInterfaces[i]->getCellRight()->getTransport(t).getValue(), buffTransportMin[t]);
        buffTransportMax[t] = std::max(m_cellInterfaces[i]->getCellLeft()->getTransport(t).getValue(), buffTransportMax[t]);
        buffTransportMax[t] = std::max(m_cellInterfaces[i]->getCellRight()->getTransport(t).getValue(), buffTransportMax[t]);
      }
    } 
  }

  // 2) Compute gradient limiter of the cell theta_i
  
  // Reset cell gradient limiter theta_i cell to high values for search on minimum cell interfaces theta_ij
  for (int k = 0; k < numberPhases; k++) {
		slopesPhasesLocal2[k]->setToMax();
	}
	slopesMixtureLocal2->setToMax();
  for (int t = 0; t < numberTransports; t++) {
    slopesTransportLocal2[t] = 1.e15;
  }
  
  Coord rij;
  for (unsigned int i = 0; i < m_cellInterfaces.size(); i++) {
    // Compute rij * grad(Wi) => slopesPhases1/MixtureLocal1/slopesTransportLocal1
    rij = m_cellInterfaces[i]->getFace()->getPos() - this->getPosition();

    for (int k = 0; k < numberPhases; k++) {
      m_gradPhase[k]->computeDistanceGradientScalarProduct(rij, slopesPhasesLocal1[k]);
    }
    m_gradMixture->computeDistanceGradientScalarProduct(rij, slopesMixtureLocal1);
    for (int t = 0; t < numberTransports; t++) {
      m_gradTransport[t].computeDistanceGradientScalarProduct(rij, slopesTransportLocal1[t]);
    }

    // Compute interface gradient limiter theta_ij and update cell gradient limiter theta_i of cell with minimum of theta_ij
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal2[k]->computeGradientLimiter(globalLimiter, *(this->getPhase(k)), *buffPhasesMin[k], *buffPhasesMax[k], *slopesPhasesLocal1[k]);
    }
    slopesMixtureLocal2->computeGradientLimiter(globalLimiter, *(this->getMixture()), *buffMixtureMin, *buffMixtureMax, *slopesMixtureLocal1);
    for (int t = 0; t < numberTransports; t++) {
      slopesTransportLocal2[t] = std::min(slopesTransportLocal2[t], globalLimiter.computeGradientLimiter(this->getTransport(t).getValue(), buffTransportMin[t], buffTransportMax[t], slopesTransportLocal1[t]));
    }
  }
  
  // 3) Limit gradient with gradientLimiter * grad(Wi) and stores it in gradient to make it ready for extrapolation
  for (int k = 0; k < numberPhases; k++) {
    m_gradPhase[k]->limitGradients(*slopesPhasesLocal2[k]);
  }
  m_gradMixture->limitGradients(*slopesMixtureLocal2);
  for (int t = 0; t < numberTransports; t++) {
    m_gradTransport[t].limitGradients(slopesTransportLocal2[t]);
  }
}

//***********************************************************************

void CellO2NS::computeLocalSlopes(CellInterface& cellInterfaceRef)
{
  // Compute gradientLimiter * rij . grad(Wi) and stores it in
  // slopesPhases1/MixtureLocal1/slopesTransportLocal1 for extrapolation
  Coord rij = cellInterfaceRef.getFace()->getPos() - this->getPosition();

  for (int k = 0; k < numberPhases; k++) {
    m_gradPhase[k]->computeDistanceGradientScalarProduct(rij, slopesPhasesLocal1[k]);
  }
  m_gradMixture->computeDistanceGradientScalarProduct(rij, slopesMixtureLocal1);

  for (int t = 0; t < numberTransports; t++) {
    m_gradTransport[t].computeDistanceGradientScalarProduct(rij, slopesTransportLocal1[t]);
  }
}

//***********************************************************************

GradPhase* CellO2NS::getGradPhase(const int& phaseNumber) const
{
  return m_gradPhase[phaseNumber];
}

//***********************************************************************

GradMixture* CellO2NS::getGradMixture() const 
{
  return m_gradMixture;
}

//***********************************************************************

GradTransport* CellO2NS::getGradTransport(const int& transportNumber) const 
{
  return &m_gradTransport[transportNumber];
}

//***********************************************************************

void CellO2NS::fillBufferSlopes(double* buffer, int& counter, const int& /*lvl*/, const int& /*neighbour*/) const
{
  //Fill buffer to send
  for (int k = 0; k < numberPhases; k++) {
    m_gradPhase[k]->fillBufferGradients(buffer, counter);
  }
  m_gradMixture->fillBufferGradients(buffer, counter);
  for (int t = 0; t < numberTransports; t++) {
    m_gradTransport[t].fillBufferGradients(buffer, counter);
  }
}

//***********************************************************************