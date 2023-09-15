#include "GradPhasePTUEq.h"

//***************************************************************************

GradPhasePTUEq::GradPhasePTUEq()
{
  m_grads.resize(1);
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i] = 0.;
  }
}

//***************************************************************************

GradPhasePTUEq::~GradPhasePTUEq()
{
}

//***************************************************************************

void GradPhasePTUEq::initializeGradientVectors()
{
  this->initializeGradsVariablesNamesNumerators();

  variableNamesPhases[VarLocal::alpha] = Variable::alpha;
}

//***************************************************************************

void GradPhasePTUEq::computeDistanceGradientScalarProduct(Coord const& distance, Phase* phase) const
{
  static_cast<PhasePTUEq*> (phase)->setAlpha(distance.scalar(m_grads[VarLocal::alpha]));
}

//***************************************************************************

void GradPhasePTUEq::limitGradients(const Phase& gradientLimiter)
{
  m_grads[VarLocal::alpha].setX(m_grads[VarLocal::alpha].getX() * gradientLimiter.getAlpha());
  m_grads[VarLocal::alpha].setY(m_grads[VarLocal::alpha].getY() * gradientLimiter.getAlpha());
  m_grads[VarLocal::alpha].setZ(m_grads[VarLocal::alpha].getZ() * gradientLimiter.getAlpha());
}

//***************************************************************************
//************************** ORDER 2 PARALLEL *******************************
//***************************************************************************

int GradPhasePTUEq::numberOfTransmittedGradients() const
{
  return 3;
}

//***************************************************************************