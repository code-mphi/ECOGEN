#include "GradPhase.h"
#include "../Order2/CellO2NS.h"

std::vector<Variable> variableNamesPhases;
std::vector<std::vector<int>> numeratorPhases;

//***************************************************************************

GradPhase::GradPhase() {}

//***************************************************************************

GradPhase::~GradPhase() {}

//***************************************************************************

void GradPhase::initializeGradsVariablesNamesNumerators()
{
  numeratorPhases.resize(numberPhases);
  for (int i = 0; i < numberPhases; ++i) {
    numeratorPhases[i].resize(m_grads.size());
    for (unsigned int j = 0; j < m_grads.size(); ++j) {
      numeratorPhases[i][j] = i;
    }
  }

  variableNamesPhases.resize(m_grads.size());
}

//***************************************************************************

void GradPhase::computeGradients(Cell* cell, int const& phase)
{
  cell->computeGradients(m_grads, variableNamesPhases, numeratorPhases[phase]);
}

//***************************************************************************
//************************** ORDER 2 PARALLEL *******************************
//***************************************************************************

void GradPhase::getBufferGradients(double * buffer, int& counter)
{
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    m_grads[i].setX(buffer[++counter]);
    m_grads[i].setY(buffer[++counter]);
    m_grads[i].setZ(buffer[++counter]);
  }
}

//***************************************************************************

void GradPhase::fillBufferGradients(double * buffer, int& counter)
{
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    buffer[++counter] = m_grads[i].getX();
    buffer[++counter] = m_grads[i].getY();
    buffer[++counter] = m_grads[i].getZ();
  }
}

//***************************************************************************