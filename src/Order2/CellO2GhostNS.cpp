#include "CellO2GhostNS.h"

//***************************************************************************

CellO2GhostNS::CellO2GhostNS()
{
}

//***************************************************************************

CellO2GhostNS::~CellO2GhostNS()
{
}

//***************************************************************************

int CellO2GhostNS::getRankOfNeighborCPU() const
{
  return m_rankOfNeighborCPU;
}

//***************************************************************************

void CellO2GhostNS::setRankOfNeighborCPU(int rank)
{
  m_rankOfNeighborCPU = rank;
}

//***************************************************************************

void CellO2GhostNS::getBufferSlopes(double* buffer, int& counter, const int& /*lvl*/)
{
  // Fill buffer to send
  for (int k = 0; k < numberPhases; k++) {
    m_gradPhase[k]->getBufferGradients(buffer, counter);
  }
  m_gradMixture->getBufferGradients(buffer, counter);
  for (int t = 0; t < numberTransports; t++) {
    m_gradTransport[t].getBufferGradients(buffer, counter);
  }
}

//***************************************************************************