#include "GradPhaseEulerHomogeneous.h"

//***************************************************************************

GradPhaseEulerHomogeneous::GradPhaseEulerHomogeneous()
{
  m_grads.resize(1);
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i] = 0.;
  }
}

//***************************************************************************

GradPhaseEulerHomogeneous::~GradPhaseEulerHomogeneous()
{
}

//***************************************************************************

void GradPhaseEulerHomogeneous::initializeGradientVectors()
{
  this->initializeGradsVariablesNamesNumerators();

  variableNamesPhases[VarLocal::alpha] = Variable::alpha;
}

//***************************************************************************

void GradPhaseEulerHomogeneous::computeDistanceGradientScalarProduct(Coord const& distance, Phase* phase) const
{
  static_cast<PhaseEulerHomogeneous*> (phase)->setAlpha(distance.scalar(m_grads[VarLocal::alpha]));
}

//***************************************************************************

void GradPhaseEulerHomogeneous::limitGradients(const Phase& gradientLimiter)
{
  m_grads[VarLocal::alpha].setX(m_grads[VarLocal::alpha].getX() * gradientLimiter.getAlpha());
  m_grads[VarLocal::alpha].setY(m_grads[VarLocal::alpha].getY() * gradientLimiter.getAlpha());
  m_grads[VarLocal::alpha].setZ(m_grads[VarLocal::alpha].getZ() * gradientLimiter.getAlpha());
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int GradPhaseEulerHomogeneous::numberOfTransmittedGradients() const
{
  return 3;
}

//****************************************************************************