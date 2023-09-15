#include "GradMixUEq.h"

//***************************************************************************

GradMixUEq::GradMixUEq()
{
  m_grads.resize(3);
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i] = 0.;
  }
}

//***************************************************************************

GradMixUEq::~GradMixUEq()
{
}

//***************************************************************************

void GradMixUEq::initializeGradientVectors()
{
  this->initializeGradsVariablesNamesNumerators();

  variableNamesMixture[VarLocal::velocityU] = Variable::velocityU;
  variableNamesMixture[VarLocal::velocityV] = Variable::velocityV;
  variableNamesMixture[VarLocal::velocityW] = Variable::velocityW;
}

//***************************************************************************

void GradMixUEq::computeDistanceGradientScalarProduct(Coord const& distance, Mixture* mixture) const
{
  static_cast<MixUEq*> (mixture)->setVelocity(distance.scalar(m_grads[VarLocal::velocityU]),
                                              distance.scalar(m_grads[VarLocal::velocityV]),
                                              distance.scalar(m_grads[VarLocal::velocityW]));
}

//***************************************************************************

void GradMixUEq::limitGradients(const Mixture& gradientLimiter)
{
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

//***************************************************************************
//************************** ORDER 2 PARALLEL *******************************
//***************************************************************************

int GradMixUEq::numberOfTransmittedGradients() const
{
  return 9;
}

//***************************************************************************