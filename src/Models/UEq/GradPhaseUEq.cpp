#include "GradPhaseUEq.h"

//***************************************************************************

GradPhaseUEq::GradPhaseUEq()
{
  m_grads.resize(3);
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i] = 0.;
  }
}

//***************************************************************************

GradPhaseUEq::~GradPhaseUEq()
{
}

//***************************************************************************

void GradPhaseUEq::initializeGradientVectors()
{
  this->initializeGradsVariablesNamesNumerators();

  variableNamesPhases[VarLocal::alpha] = Variable::alpha;
  variableNamesPhases[VarLocal::density] = Variable::density;
  variableNamesPhases[VarLocal::pressure] = Variable::pressure;
}

//***************************************************************************

void GradPhaseUEq::computeDistanceGradientScalarProduct(Coord const& distance, Phase* phase) const
{
  static_cast<PhaseUEq*> (phase)->setAlpha(distance.scalar(m_grads[VarLocal::alpha]));
  static_cast<PhaseUEq*> (phase)->setDensity(distance.scalar(m_grads[VarLocal::density]));
  static_cast<PhaseUEq*> (phase)->setPressure(distance.scalar(m_grads[VarLocal::pressure]));
}

//***************************************************************************

void GradPhaseUEq::limitGradients(const Phase& gradientLimiter)
{
  m_grads[VarLocal::alpha].setX(m_grads[VarLocal::alpha].getX() * gradientLimiter.getAlpha());
  m_grads[VarLocal::alpha].setY(m_grads[VarLocal::alpha].getY() * gradientLimiter.getAlpha());
  m_grads[VarLocal::alpha].setZ(m_grads[VarLocal::alpha].getZ() * gradientLimiter.getAlpha());

  m_grads[VarLocal::density].setX(m_grads[VarLocal::density].getX() * gradientLimiter.getDensity());
  m_grads[VarLocal::density].setY(m_grads[VarLocal::density].getY() * gradientLimiter.getDensity());
  m_grads[VarLocal::density].setZ(m_grads[VarLocal::density].getZ() * gradientLimiter.getDensity());

  m_grads[VarLocal::pressure].setX(m_grads[VarLocal::pressure].getX() * gradientLimiter.getPressure());
  m_grads[VarLocal::pressure].setY(m_grads[VarLocal::pressure].getY() * gradientLimiter.getPressure());
  m_grads[VarLocal::pressure].setZ(m_grads[VarLocal::pressure].getZ() * gradientLimiter.getPressure());
}

//***************************************************************************
//************************** ORDER 2 PARALLEL *******************************
//***************************************************************************

int GradPhaseUEq::numberOfTransmittedGradients() const
{
  return 9;
}

//***************************************************************************