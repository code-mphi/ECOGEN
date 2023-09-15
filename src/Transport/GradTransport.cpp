#include "GradTransport.h"
#include "../Order2/CellO2NS.h"

std::vector<Variable> variableNamesTransports;
std::vector<std::vector<int>> numeratorTransports;

//***************************************************************************

GradTransport::GradTransport() : m_grads(1)
{
  m_grads[0] = 0.;
}

//***************************************************************************

GradTransport::~GradTransport()
{
}

//***************************************************************************

void GradTransport::initializeGradientVectors()
{
  numeratorTransports.resize(numberTransports);
  for (int i = 0; i < numberTransports; ++i) {
    numeratorTransports[i].resize(m_grads.size());
    for (unsigned int j = 0; j < m_grads.size(); ++j) {
      numeratorTransports[i][j] = i;
    }
  }

  variableNamesTransports.resize(m_grads.size());
  variableNamesTransports[0] = Variable::transport;
}

//***************************************************************************

void GradTransport::computeGradient(Cell* cell, const int& numTransport)
{
  cell->computeGradients(m_grads, variableNamesTransports, numeratorTransports[numTransport]);
}

//***************************************************************************

void GradTransport::computeDistanceGradientScalarProduct(Coord const& distance, double &transport) const
{
  transport = distance.scalar(m_grads[0]);
}

//****************************************************************************

void GradTransport::limitGradients(const double& gradientLimiter)
{
  m_grads[0].setX(m_grads[0].getX() * gradientLimiter);
  m_grads[0].setY(m_grads[0].getY() * gradientLimiter);
  m_grads[0].setZ(m_grads[0].getZ() * gradientLimiter);
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int GradTransport::numberOfTransmittedGradients() const
{
  return 3;
}

//***************************************************************************

void GradTransport::getBufferGradients(double * buffer, int& counter)
{
  m_grads[0].setX(buffer[++counter]);
  m_grads[0].setY(buffer[++counter]);
  m_grads[0].setZ(buffer[++counter]);
}

//***************************************************************************

void GradTransport::fillBufferGradients(double * buffer, int& counter)
{
  buffer[++counter] = m_grads[0].getX();
  buffer[++counter] = m_grads[0].getY();
  buffer[++counter] = m_grads[0].getZ();
}

//***************************************************************************