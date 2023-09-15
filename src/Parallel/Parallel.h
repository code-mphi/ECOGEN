//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//      Official webSite: https://code-mphi.github.io/ECOGEN/
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef PARALLEL_H
#define PARALLEL_H

#include <mpi.h>
#include "../Tools.h"
#include "../Models/Phase.h"
#include "../Order1/Cell.h"

class Parallel
{
public:
  Parallel();
  ~Parallel();

  void initialization();
  void setNeighbour(const int neighbour);
  void addElementToSend(int neighbour, Cell* cell);
  void addElementToReceive(int neighbour, Cell* cell);
  void addSlopesToSend(int neighbour);
  void addSlopesToReceive(int neighbour);
  void deleteSlopesToSend(int neighbour);
  void deleteSlopesToReceive(int neighbour);
  void clearElementsAndSlopesToSendAndReceivePLusNeighbour();
  TypeMeshContainer<Cell*>& getElementsToSend(int neighbour);
  TypeMeshContainer<Cell*>& getElementsToReceive(int neighbour);
  void initializePersistentCommunications(const int& numberPrimitiveVariables, const int& numberSlopeVariables, const int& numberTransportVariables, const int& dim);
  void computeDt(double& dt);
  void computePMax(double& pMax, double& pMaxWall);
  void computeSum(double& var);
  void finalize(const int& lvlMax);
  void stopRun();
  bool verifyStateCPUs();
  
  //Methodes pour toutes les variables primitives
  void initializePersistentCommunicationsPrimitives();
  void finalizePersistentCommunicationsPrimitives(const int& lvlMax);
  void communicationsPrimitives(Eos** eos, int lvl, Prim type = vecPhases);

  //Methodes pour toutes les slopes
  void initializePersistentCommunicationsSlopes();
  void finalizePersistentCommunicationsSlopes(const int& lvlMax);
  void communicationsSlopes(int lvl);

  //Methodes pour une variable scalar
  void initializePersistentCommunicationsScalar();
  void finalizePersistentCommunicationsScalar(const int& lvlMax);

  //Methodes pour une variable vectorielle
  void initializePersistentCommunicationsVector(const int& dim);
  void finalizePersistentCommunicationsVector(const int& lvlMax);
  void communicationsVector(Variable nameVector, const int& dim, int lvl, int num = 0, int index = -1);

  //Methodes pour toutes les variables transports
  void initializePersistentCommunicationsTransports();
  void finalizePersistentCommunicationsTransports(const int& lvlMax);
  void communicationsTransports(int lvl);

  //Methodes pour les variables AMR
  void initializePersistentCommunicationsAMR(const int& numberPrimitiveVariables, const int& numberSlopeVariables, const int& numberTransportVariables, const int& dim, const int& lvlMax);
  void initializePersistentCommunicationsLvlAMR(const int& lvlMax);
  void clearRequestsAndBuffers(int lvl);
  void updatePersistentCommunicationsAMR(const int& dim);
  void updatePersistentCommunicationsLvlAMR(int lvl, const int& dim);
  void finalizeAMR(const int& lvlMax);

  void initializePersistentCommunicationsXi();
  void finalizePersistentCommunicationsXi(const int& lvlMax);
  void communicationsXi(int lvl);

  void initializePersistentCommunicationsSplit();
  void finalizePersistentCommunicationsSplit(const int& lvlMax);
  void communicationsSplit(int lvl);

  void initializePersistentCommunicationsNumberGhostCells();
  void finalizePersistentCommunicationsNumberGhostCells();
  void communicationsNumberGhostCells(int lvl);

private:
    
  bool* m_isNeighbour;
  std::vector<TypeMeshContainer<Cell*>> m_elementsToSend;
  std::vector<TypeMeshContainer<Cell*>> m_elementsToReceive;
  int* m_numberElementsToSendToNeighbour;
  int* m_numberElementsToReceiveFromNeighbour;
  int* m_numberSlopesToSendToNeighbour;
  int* m_numberSlopesToReceiveFromNeighbour;
  int m_numberPrimitiveVariables;          /*Number of primitive variables to send (phases + mixture + transports)*/
  int m_numberSlopeVariables;              /*Number of slope variables to send (phases + mixture + transports)*/
  int m_numberTransportVariables;          /*Number of transport variables to send*/

  std::vector<double**> m_bufferReceive;
  std::vector<double**> m_bufferSend;
  std::vector<double**> m_bufferReceiveSlopes;
  std::vector<double**> m_bufferSendSlopes;
  std::vector<double**> m_bufferReceiveScalar;
  std::vector<double**> m_bufferSendScalar;
  std::vector<double**> m_bufferReceiveVector;
  std::vector<double**> m_bufferSendVector;
  std::vector<double**> m_bufferReceiveTransports;
  std::vector<double**> m_bufferSendTransports;
  std::vector<double**> m_bufferReceiveXi;
  std::vector<double**> m_bufferSendXi;
  std::vector<bool**> m_bufferReceiveSplit;
  std::vector<bool**> m_bufferSendSplit;
  int* m_bufferNumberElementsToSendToNeighbor;
  int* m_bufferNumberElementsToReceiveFromNeighbour;
  int* m_bufferNumberSlopesToSendToNeighbor;
  int* m_bufferNumberSlopesToReceiveFromNeighbour;
  
  std::vector<MPI_Request**> m_reqSend;
  std::vector<MPI_Request**> m_reqReceive;
  std::vector<MPI_Request**> m_reqSendSlopes;
  std::vector<MPI_Request**> m_reqReceiveSlopes;
  std::vector<MPI_Request**> m_reqSendScalar;
  std::vector<MPI_Request**> m_reqReceiveScalar;
  std::vector<MPI_Request**> m_reqSendVector;
  std::vector<MPI_Request**> m_reqReceiveVector;
  std::vector<MPI_Request**> m_reqSendTransports;
  std::vector<MPI_Request**> m_reqReceiveTransports;
  std::vector<MPI_Request**> m_reqSendXi;
  std::vector<MPI_Request**> m_reqReceiveXi;
  std::vector<MPI_Request**> m_reqSendSplit;
  std::vector<MPI_Request**> m_reqReceiveSplit;
  MPI_Request** m_reqNumberElementsToSendToNeighbor;
  MPI_Request** m_reqNumberElementsToReceiveFromNeighbour;
  MPI_Request** m_reqNumberSlopesToSendToNeighbor;
  MPI_Request** m_reqNumberSlopesToReceiveFromNeighbour;

};

extern Parallel parallel;
extern int rankCpu;
extern int Ncpu;

#endif // PARALLEL_H