#include "GradPhaseEuler.h"

//***************************************************************************

GradPhaseEuler::GradPhaseEuler()
{
  m_grads.resize(5);
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i] = 0.;
  }
}

//***************************************************************************

GradPhaseEuler::~GradPhaseEuler()
{
}

//***************************************************************************

void GradPhaseEuler::initializeGradientVectors()
{
  this->initializeGradsVariablesNamesNumerators();

  variableNamesPhases[VarLocal::density] = Variable::density;
  variableNamesPhases[VarLocal::pressure] = Variable::pressure;
  variableNamesPhases[VarLocal::velocityU] = Variable::velocityU;
  variableNamesPhases[VarLocal::velocityV] = Variable::velocityV;
  variableNamesPhases[VarLocal::velocityW] = Variable::velocityW;
}

//***************************************************************************

void GradPhaseEuler::computeDistanceGradientScalarProduct(Coord const& distance, Phase* phase) const
{
  static_cast<PhaseEuler*> (phase)->setDensity(distance.scalar(m_grads[VarLocal::density]));
  static_cast<PhaseEuler*> (phase)->setPressure(distance.scalar(m_grads[VarLocal::pressure]));
  static_cast<PhaseEuler*> (phase)->setVelocity(distance.scalar(m_grads[VarLocal::velocityU]),
                                                distance.scalar(m_grads[VarLocal::velocityV]),
                                                distance.scalar(m_grads[VarLocal::velocityW]));
}

//****************************************************************************

void GradPhaseEuler::limitGradients(const Phase& gradientLimiter)
{
  m_grads[VarLocal::density].setX(m_grads[VarLocal::density].getX() * gradientLimiter.getDensity());
  m_grads[VarLocal::density].setY(m_grads[VarLocal::density].getY() * gradientLimiter.getDensity());
  m_grads[VarLocal::density].setZ(m_grads[VarLocal::density].getZ() * gradientLimiter.getDensity());

  m_grads[VarLocal::pressure].setX(m_grads[VarLocal::pressure].getX() * gradientLimiter.getPressure());
  m_grads[VarLocal::pressure].setY(m_grads[VarLocal::pressure].getY() * gradientLimiter.getPressure());
  m_grads[VarLocal::pressure].setZ(m_grads[VarLocal::pressure].getZ() * gradientLimiter.getPressure());

  m_grads[VarLocal::velocityU].setX(m_grads[VarLocal::velocityU].getX() * gradientLimiter.getVelocity().getX());
  m_grads[VarLocal::velocityU].setY(m_grads[VarLocal::velocityU].getY() * gradientLimiter.getVelocity().getX());
  m_grads[VarLocal::velocityU].setZ(m_grads[VarLocal::velocityU].getZ() * gradientLimiter.getVelocity().getX());

  m_grads[VarLocal::velocityV].setX(m_grads[VarLocal::velocityV].getX() * gradientLimiter.getVelocity().getY());
  m_grads[VarLocal::velocityV].setY(m_grads[VarLocal::velocityV].getY() * gradientLimiter.getVelocity().getY());
  m_grads[VarLocal::velocityV].setZ(m_grads[VarLocal::velocityV].getZ() * gradientLimiter.getVelocity().getY());

  m_grads[VarLocal::velocityW].setX(m_grads[VarLocal::velocityW].getX() * gradientLimiter.getVelocity().getZ());
  m_grads[VarLocal::velocityW].setY(m_grads[VarLocal::velocityW].getY() * gradientLimiter.getVelocity().getZ());
  m_grads[VarLocal::velocityW].setZ(m_grads[VarLocal::velocityW].getZ() * gradientLimiter.getVelocity().getZ());
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int GradPhaseEuler::numberOfTransmittedGradients() const
{
  return 15;
}

//****************************************************************************