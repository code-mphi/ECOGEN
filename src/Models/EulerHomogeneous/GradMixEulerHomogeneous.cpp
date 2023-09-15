#include "GradMixEulerHomogeneous.h"

//***************************************************************************

GradMixEulerHomogeneous::GradMixEulerHomogeneous()
{
  m_grads.resize(4);
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i] = 0.;
  }
}

//***************************************************************************

GradMixEulerHomogeneous::~GradMixEulerHomogeneous()
{
}

//***************************************************************************

void GradMixEulerHomogeneous::initializeGradientVectors()
{
  this->initializeGradsVariablesNamesNumerators();

  variableNamesMixture[VarLocal::pressure] = Variable::pressure;
  variableNamesMixture[VarLocal::velocityU] = Variable::velocityU;
  variableNamesMixture[VarLocal::velocityV] = Variable::velocityV;
  variableNamesMixture[VarLocal::velocityW] = Variable::velocityW;
}

//***************************************************************************

void GradMixEulerHomogeneous::computeDistanceGradientScalarProduct(Coord const& distance, Mixture* mixture) const
{
  static_cast<MixEulerHomogeneous*> (mixture)->setPressure(distance.scalar(m_grads[VarLocal::pressure]));
  static_cast<MixEulerHomogeneous*> (mixture)->setVelocity(distance.scalar(m_grads[VarLocal::velocityU]),
                                                           distance.scalar(m_grads[VarLocal::velocityV]),
                                                           distance.scalar(m_grads[VarLocal::velocityW]));
}

//***************************************************************************

void GradMixEulerHomogeneous::limitGradients(const Mixture& gradientLimiter)
{
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

int GradMixEulerHomogeneous::numberOfTransmittedGradients() const
{
  return 12;
}

//****************************************************************************