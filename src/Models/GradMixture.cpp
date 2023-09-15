#include "GradMixture.h"
#include "../Order2/CellO2NS.h"

std::vector<Variable> variableNamesMixture;
std::vector<int> numeratorMixture;

//***************************************************************************

GradMixture::GradMixture() {}

//***************************************************************************

GradMixture::~GradMixture() {}

//***************************************************************************

void GradMixture::initializeGradsVariablesNamesNumerators()
{
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i] = 0.;
  }

  numeratorMixture.resize(m_grads.size());
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    numeratorMixture[i] = -1;
  }

  variableNamesMixture.resize(m_grads.size());
}

//***************************************************************************

void GradMixture::computeGradients(Cell* cell)
{
  cell->computeGradients(m_grads, variableNamesMixture, numeratorMixture);
}

//***************************************************************************
//************************** ORDER 2 PARALLEL *******************************
//***************************************************************************

void GradMixture::getBufferGradients(double * buffer, int& counter)
{
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i].setX(buffer[++counter]);
    m_grads[i].setY(buffer[++counter]);
    m_grads[i].setZ(buffer[++counter]);
  }
}

//***************************************************************************

void GradMixture::fillBufferGradients(double * buffer, int& counter)
{
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    buffer[++counter] = m_grads[i].getX();
    buffer[++counter] = m_grads[i].getY();
    buffer[++counter] = m_grads[i].getZ();
  }
}

//***************************************************************************