//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

//! \file      Parallel.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      November 11 2019

#include "Parallel.h"
#include "../Eos/Eos.h"

//Variables linked to parallel computation
Parallel parallel;
int rankCpu, Ncpu;

//***********************************************************************

Parallel::Parallel(): m_stateCPU(1) {}

//***********************************************************************

Parallel::~Parallel(){}

//***********************************************************************

void Parallel::initialization(int &argc, char* argv[])
{
  if (Ncpu == 1) return; //The following is not necessary in the case of monoCPU

  m_isNeighbour = new bool[Ncpu];
  m_elementsToSend.resize(Ncpu);
  m_elementsToReceive.resize(Ncpu);
  m_numberElementsToSendToNeighbour = new int[Ncpu];
  m_numberElementsToReceiveFromNeighbour = new int[Ncpu];
  m_numberSlopesToSendToNeighbour = new int[Ncpu];
  m_numberSlopesToReceiveFromNeighbour = new int[Ncpu];

  m_bufferSend.push_back(new double*[Ncpu]);
  m_bufferReceive.push_back(new double*[Ncpu]);
  m_bufferSendSlopes.push_back(new double*[Ncpu]);
  m_bufferReceiveSlopes.push_back(new double*[Ncpu]);
  m_bufferSendScalar.push_back(new double*[Ncpu]);
  m_bufferReceiveScalar.push_back(new double*[Ncpu]);
  m_bufferSendVector.push_back(new double*[Ncpu]);
  m_bufferReceiveVector.push_back(new double*[Ncpu]);
  m_bufferSendTransports.push_back(new double*[Ncpu]);
  m_bufferReceiveTransports.push_back(new double*[Ncpu]);
  m_bufferSendXi.push_back(new double*[Ncpu]);
  m_bufferReceiveXi.push_back(new double*[Ncpu]);
  m_bufferSendSplit.push_back(new bool*[Ncpu]);
  m_bufferReceiveSplit.push_back(new bool*[Ncpu]);
  m_bufferNumberElementsToSendToNeighbor = new int[Ncpu];
  m_bufferNumberElementsToReceiveFromNeighbour = new int[Ncpu];
  m_bufferNumberSlopesToSendToNeighbor = new int[Ncpu];
  m_bufferNumberSlopesToReceiveFromNeighbour = new int[Ncpu];

  m_reqSend.push_back(new MPI_Request*[Ncpu]);
  m_reqReceive.push_back(new MPI_Request*[Ncpu]);
  m_reqSendSlopes.push_back(new MPI_Request*[Ncpu]);
  m_reqReceiveSlopes.push_back(new MPI_Request*[Ncpu]);
  m_reqSendScalar.push_back(new MPI_Request*[Ncpu]);
  m_reqReceiveScalar.push_back(new MPI_Request*[Ncpu]);
  m_reqSendVector.push_back(new MPI_Request*[Ncpu]);
  m_reqReceiveVector.push_back(new MPI_Request*[Ncpu]);
  m_reqSendTransports.push_back(new MPI_Request*[Ncpu]);
  m_reqReceiveTransports.push_back(new MPI_Request*[Ncpu]);
  m_reqSendXi.push_back(new MPI_Request*[Ncpu]);
  m_reqReceiveXi.push_back(new MPI_Request*[Ncpu]);
  m_reqSendSplit.push_back(new MPI_Request*[Ncpu]);
  m_reqReceiveSplit.push_back(new MPI_Request*[Ncpu]);
  m_reqNumberElementsToSendToNeighbor = new MPI_Request*[Ncpu];
  m_reqNumberElementsToReceiveFromNeighbour = new MPI_Request*[Ncpu];
  m_reqNumberSlopesToSendToNeighbor = new MPI_Request*[Ncpu];
  m_reqNumberSlopesToReceiveFromNeighbour = new MPI_Request*[Ncpu];

  for (int i = 0; i < Ncpu; i++) {
    m_isNeighbour[i] = false;
    m_numberElementsToSendToNeighbour[i] = 0;
    m_numberElementsToReceiveFromNeighbour[i] = 0;
    m_numberSlopesToSendToNeighbour[i] = 0;
    m_numberSlopesToReceiveFromNeighbour[i] = 0;
    m_bufferSend[0][i] = NULL;
    m_bufferReceive[0][i] = NULL;
    m_reqSend[0][i] = NULL;
    m_reqReceive[0][i] = NULL;
    m_bufferSendSlopes[0][i] = NULL;
    m_bufferReceiveSlopes[0][i] = NULL;
    m_reqSendSlopes[0][i] = NULL;
    m_reqReceiveSlopes[0][i] = NULL;
    m_bufferSendScalar[0][i] = NULL;
    m_bufferReceiveScalar[0][i] = NULL;
    m_reqSendScalar[0][i] = NULL;
    m_reqReceiveScalar[0][i] = NULL;
    m_bufferSendVector[0][i] = NULL;
    m_bufferReceiveVector[0][i] = NULL;
    m_reqSendVector[0][i] = NULL;
    m_reqReceiveVector[0][i] = NULL;
    m_bufferSendTransports[0][i] = NULL;
    m_bufferReceiveTransports[0][i] = NULL;
    m_reqSendTransports[0][i] = NULL;
    m_reqReceiveTransports[0][i] = NULL;
    m_bufferSendXi[0][i] = NULL;
    m_bufferReceiveXi[0][i] = NULL;
    m_reqSendXi[0][i] = NULL;
    m_reqReceiveXi[0][i] = NULL;
    m_bufferSendSplit[0][i] = NULL;
    m_bufferReceiveSplit[0][i] = NULL;
    m_reqSendSplit[0][i] = NULL;
    m_reqReceiveSplit[0][i] = NULL;
    m_bufferNumberElementsToSendToNeighbor[i] = 0;
    m_bufferNumberElementsToReceiveFromNeighbour[i] = 0;
    m_bufferNumberSlopesToSendToNeighbor[i] = 0;
    m_bufferNumberSlopesToReceiveFromNeighbour[i] = 0;
    m_reqNumberElementsToSendToNeighbor[i] = NULL;
    m_reqNumberElementsToReceiveFromNeighbour[i] = NULL;
    m_reqNumberSlopesToSendToNeighbor[i] = NULL;
    m_reqNumberSlopesToReceiveFromNeighbour[i] = NULL;
  }
}

//***********************************************************************

void Parallel::setNeighbour(const int neighbour)
{ 
  m_isNeighbour[neighbour] = true;
}

//***********************************************************************

void Parallel::addElementToSend(int neighbour, Cell* cell)
{
  m_elementsToSend[neighbour].push_back(cell);
  m_numberElementsToSendToNeighbour[neighbour]=m_elementsToSend[neighbour].size();
}

//***********************************************************************

void Parallel::addElementToReceive(int neighbour, Cell* cell)
{
  m_elementsToReceive[neighbour].push_back(cell);
  m_numberElementsToReceiveFromNeighbour[neighbour]=m_elementsToReceive[neighbour].size();
}

//***********************************************************************

void Parallel::addSlopesToSend(int neighbour)
{
  m_numberSlopesToSendToNeighbour[neighbour] += 1;
}

//***********************************************************************

void Parallel::addSlopesToReceive(int neighbour)
{
  m_numberSlopesToReceiveFromNeighbour[neighbour] += 1;
}

//***********************************************************************

void Parallel::clearElementsAndSlopesToSendAndReceivePLusNeighbour()
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    m_elementsToSend[neighbour].clear();
    m_elementsToReceive[neighbour].clear();
    m_numberSlopesToSendToNeighbour[neighbour] = 0;
    m_numberSlopesToReceiveFromNeighbour[neighbour] = 0;
    m_isNeighbour[neighbour] = false;
  }
}

//***********************************************************************

const TypeMeshContainer<Cell*> &Parallel::getElementsToSend(int neighbour) const
{
  return m_elementsToSend[neighbour];
}

//***********************************************************************

TypeMeshContainer<Cell*> &Parallel::getElementsToSend(int neighbour)
{
  return m_elementsToSend[neighbour];
}

//***********************************************************************

TypeMeshContainer<Cell*> &Parallel::getElementsToReceive(int neighbour)
{
  return m_elementsToReceive[neighbour];
}

//***********************************************************************

void Parallel::initializePersistentCommunications(const int &numberPrimitiveVariables, const int &numberSlopeVariables, const int &numberTransportVariables, const int &dim)
{
  if (Ncpu > 1) {
    m_numberPrimitiveVariables = numberPrimitiveVariables;
    m_numberSlopeVariables = numberSlopeVariables;
    m_numberTransportVariables = numberTransportVariables;
    //Initialization of communications of primitive variables from resolved model
    parallel.initializePersistentCommunicationsPrimitives();
    //Initialization of communications of slopes for second order
    parallel.initializePersistentCommunicationsSlopes();
    //Initialization of communications necessary for additional physics (vectors of dim=3)
    parallel.initializePersistentCommunicationsVector(dim);
    //Initialization of communications of transported variables
    parallel.initializePersistentCommunicationsTransports();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::computeDt(double &dt)
{
  double dt_temp = dt;
  MPI_Allreduce(&dt_temp, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::computePMax(double &pMax, double &pMaxWall)
{
  double pMax_temp(pMax), pMaxWall_temp(pMaxWall);
  MPI_Allreduce(&pMax_temp, &pMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&pMaxWall_temp, &pMaxWall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::computeMassTotal(double &mass)
{
  double massBuff(mass);
  MPI_Allreduce(&massBuff, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::finalize(const int &lvlMax)
{
  if (Ncpu > 1) {
    this->finalizePersistentCommunicationsPrimitives(lvlMax);
    this->finalizePersistentCommunicationsSlopes(lvlMax);
    this->finalizePersistentCommunicationsVector(lvlMax);
    this->finalizePersistentCommunicationsTransports(lvlMax);
    delete[] m_isNeighbour;
    delete[] m_numberElementsToSendToNeighbour;
    delete[] m_numberElementsToReceiveFromNeighbour;
    delete[] m_numberSlopesToSendToNeighbour;
    delete[] m_numberSlopesToReceiveFromNeighbour;
    for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
        m_elementsToSend[neighbour].clear();
        m_elementsToReceive[neighbour].clear();
    }
    m_elementsToSend.clear();
    m_elementsToReceive.clear();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::stopRun()
{
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
}

//***********************************************************************

void Parallel::verifyStateCPUs()
{
  //Gathering of errors
  int nbErr_temp(0);
  int nbErr(errors.size());
  MPI_Allreduce(&nbErr, &nbErr_temp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
  //Stop if error on one CPU
  if (nbErr_temp) {
    Errors::arretCodeApresError(errors);
  }
}

//****************************************************************************
//**************** Methods for all the primitive variables *******************
//****************************************************************************

void Parallel::initializePersistentCommunicationsPrimitives()
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Determination of the number of variables to communicate
      int numberSend = m_numberPrimitiveVariables*m_numberElementsToSendToNeighbour[neighbour];
      int numberReceive = m_numberPrimitiveVariables*m_numberElementsToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSend[0][neighbour] = new MPI_Request;
      m_bufferSend[0][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSend[0][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSend[0][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceive[0][neighbour] = new MPI_Request;
      m_bufferReceive[0][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceive[0][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceive[0][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsPrimitives(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++)	{
      if (m_isNeighbour[neighbour]) {
        MPI_Request_free(m_reqSend[lvl][neighbour]);
        MPI_Request_free(m_reqReceive[lvl][neighbour]);
        delete m_reqSend[lvl][neighbour];
        delete[] m_bufferSend[lvl][neighbour];
        delete m_reqReceive[lvl][neighbour];
        delete[] m_bufferReceive[lvl][neighbour];
      }
    }
    delete[] m_reqSend[lvl];
    delete[] m_bufferSend[lvl];
    delete[] m_reqReceive[lvl];
    delete[] m_bufferReceive[lvl];
  }
  m_reqSend.clear();
  m_bufferSend.clear();
  m_reqReceive.clear();
  m_bufferReceive.clear();
}

//***********************************************************************

void Parallel::communicationsPrimitives(Eos **eos, int lvl, Prim type)
{
  int count(0);
  MPI_Status status;

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Prepation of sendings
      count = -1;
      for (int i = 0; i < m_numberElementsToSendToNeighbour[neighbour]; i++) {
        m_elementsToSend[neighbour][i]->fillBufferPrimitives(m_bufferSend[lvl][neighbour], count, lvl, neighbour, type);
      }

      //Sending request
      MPI_Start(m_reqSend[lvl][neighbour]);
      //Receiving request
      MPI_Start(m_reqReceive[lvl][neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqSend[lvl][neighbour], &status);
      MPI_Wait(m_reqReceive[lvl][neighbour], &status);

      //Receivings
      count = -1;
      for (int i = 0; i < m_numberElementsToReceiveFromNeighbour[neighbour]; i++) {
        m_elementsToReceive[neighbour][i]->getBufferPrimitives(m_bufferReceive[lvl][neighbour], count, lvl, eos, type);
      }
    }
  }
}

//****************************************************************************
//********************** Methods for all the slopes **************************
//****************************************************************************

void Parallel::initializePersistentCommunicationsSlopes()
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++)	{
    if (m_isNeighbour[neighbour]) {
      //Determination of the number of variables to communicate
      int numberSend = m_numberSlopeVariables*m_numberSlopesToSendToNeighbour[neighbour];
      int numberReceive = m_numberSlopeVariables*m_numberSlopesToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSendSlopes[0][neighbour] = new MPI_Request;
      m_bufferSendSlopes[0][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendSlopes[0][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendSlopes[0][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveSlopes[0][neighbour] = new MPI_Request;
      m_bufferReceiveSlopes[0][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveSlopes[0][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveSlopes[0][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsSlopes(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++)	{
      if (m_isNeighbour[neighbour]) {
        MPI_Request_free(m_reqSendSlopes[lvl][neighbour]);
        MPI_Request_free(m_reqReceiveSlopes[lvl][neighbour]);
        delete m_reqSendSlopes[lvl][neighbour];
        delete[] m_bufferSendSlopes[lvl][neighbour];
        delete m_reqReceiveSlopes[lvl][neighbour];
        delete[] m_bufferReceiveSlopes[lvl][neighbour];
      }
    }
    delete[] m_reqSendSlopes[lvl];
    delete[] m_bufferSendSlopes[lvl];
    delete[] m_reqReceiveSlopes[lvl];
    delete[] m_bufferReceiveSlopes[lvl];
  }
  m_reqSendSlopes.clear();
  m_bufferSendSlopes.clear();
  m_reqReceiveSlopes.clear();
  m_bufferReceiveSlopes.clear();
}

//***********************************************************************

void Parallel::communicationsSlopes(int lvl)
{
  int count(0);
  MPI_Status status;
  
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Prepation of sendings
      count = -1;
      for (int i = 0; i < m_numberElementsToSendToNeighbour[neighbour]; i++) {
        m_elementsToSend[neighbour][i]->fillBufferSlopes(m_bufferSendSlopes[lvl][neighbour], count, lvl, neighbour);
      }

      //Sending request
      MPI_Start(m_reqSendSlopes[lvl][neighbour]);
      //Receiving request
      MPI_Start(m_reqReceiveSlopes[lvl][neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqSendSlopes[lvl][neighbour], &status);
      MPI_Wait(m_reqReceiveSlopes[lvl][neighbour], &status);

      //Receivings
      count = -1;
      for (int i = 0; i < m_numberElementsToReceiveFromNeighbour[neighbour]; i++) {
        m_elementsToReceive[neighbour][i]->getBufferSlopes(m_bufferReceiveSlopes[lvl][neighbour], count, lvl);
      }
    }
  }
}

//****************************************************************************
//********************* Methods for a scalar variable ************************
//****************************************************************************

void Parallel::initializePersistentCommunicationsScalar()
{
  int number;

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Determination of the number of variables to communicate
      number = 1; //1 scalar variable
      int numberSend = number*m_numberElementsToSendToNeighbour[neighbour];
      int numberReceive = number*m_numberElementsToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSendScalar[0][neighbour] = new MPI_Request;
      m_bufferSendScalar[0][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendScalar[0][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendScalar[0][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveScalar[0][neighbour] = new MPI_Request;
      m_bufferReceiveScalar[0][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveScalar[0][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveScalar[0][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsScalar(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++)	{
      if (m_isNeighbour[neighbour]) {
        MPI_Request_free(m_reqSendScalar[lvl][neighbour]);
        MPI_Request_free(m_reqReceiveScalar[lvl][neighbour]);
        delete m_reqSendScalar[lvl][neighbour];
        delete[] m_bufferSendScalar[lvl][neighbour];
        delete m_reqReceiveScalar[lvl][neighbour];
        delete[] m_bufferReceiveScalar[lvl][neighbour];
      }
    }
    delete[] m_reqSendScalar[lvl];
    delete[] m_bufferSendScalar[lvl];
    delete[] m_reqReceiveScalar[lvl];
    delete[] m_bufferReceiveScalar[lvl];
  }
  m_reqSendScalar.clear();
  m_bufferSendScalar.clear();
  m_reqReceiveScalar.clear();
  m_bufferReceiveScalar.clear();
}

//****************************************************************************
//*********************** Methods for the vectors ****************************
//****************************************************************************

void Parallel::initializePersistentCommunicationsVector(const int &dim)
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Determination of the number of variables to communicate, as much variables as the dimension (1,2 or 3)
      int numberSend = dim*m_numberElementsToSendToNeighbour[neighbour];
      int numberReceive = dim*m_numberElementsToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSendVector[0][neighbour] = new MPI_Request;
      m_bufferSendVector[0][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendVector[0][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendVector[0][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveVector[0][neighbour] = new MPI_Request;
      m_bufferReceiveVector[0][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveVector[0][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveVector[0][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsVector(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++)	{
      if (m_isNeighbour[neighbour]) {
        MPI_Request_free(m_reqSendVector[lvl][neighbour]);
        MPI_Request_free(m_reqReceiveVector[lvl][neighbour]);
        delete m_reqSendVector[lvl][neighbour];
        delete[] m_bufferSendVector[lvl][neighbour];
        delete m_reqReceiveVector[lvl][neighbour];
        delete[] m_bufferReceiveVector[lvl][neighbour];
      }
    }
    delete[] m_reqSendVector[lvl];
    delete[] m_bufferSendVector[lvl];
    delete[] m_reqReceiveVector[lvl];
    delete[] m_bufferReceiveVector[lvl];
  }
  m_reqSendVector.clear();
  m_bufferSendVector.clear();
  m_reqReceiveVector.clear();
  m_bufferReceiveVector.clear();
}

//***********************************************************************

void Parallel::communicationsVector(Variable nameVector, const int &dim, int lvl, int num, int index)
{
  int count(0);
  MPI_Status status;

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Prepation of sendings
      count = -1;
      for (int i = 0; i < m_numberElementsToSendToNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferSendVector function of gradient coordinates
        m_elementsToSend[neighbour][i]->fillBufferVector(m_bufferSendVector[lvl][neighbour], count, lvl, neighbour, dim, nameVector, num, index);
      }

      //Sending request
      MPI_Start(m_reqSendVector[lvl][neighbour]);
      //Receiving request
      MPI_Start(m_reqReceiveVector[lvl][neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqSendVector[lvl][neighbour], &status);
      MPI_Wait(m_reqReceiveVector[lvl][neighbour], &status);
      //Receivings
      count = -1;
      for (int i = 0; i < m_numberElementsToReceiveFromNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferReceiveVector function of gradient coordinates
        m_elementsToReceive[neighbour][i]->getBufferVector(m_bufferReceiveVector[lvl][neighbour], count, lvl, dim, nameVector, num, index);
      }
    }
  }
}

//****************************************************************************
//************ Methodes pour toutes les variables transportees ***************
//****************************************************************************

void Parallel::initializePersistentCommunicationsTransports()
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Determination of the number of variables to communicate
      int numberSend = m_numberTransportVariables*m_numberElementsToSendToNeighbour[neighbour];
      int numberReceive = m_numberTransportVariables*m_numberElementsToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSendTransports[0][neighbour] = new MPI_Request;
      m_bufferSendTransports[0][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendTransports[0][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendTransports[0][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveTransports[0][neighbour] = new MPI_Request;
      m_bufferReceiveTransports[0][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveTransports[0][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveTransports[0][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsTransports(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
      if (m_isNeighbour[neighbour]) {
        MPI_Request_free(m_reqSendTransports[lvl][neighbour]);
        MPI_Request_free(m_reqReceiveTransports[lvl][neighbour]);
        delete m_reqSendTransports[lvl][neighbour];
        delete[] m_bufferSendTransports[lvl][neighbour];
        delete m_reqReceiveTransports[lvl][neighbour];
        delete[] m_bufferReceiveTransports[lvl][neighbour];
      }
    }
    delete[] m_reqSendTransports[lvl];
    delete[] m_bufferSendTransports[lvl];
    delete[] m_reqReceiveTransports[lvl];
    delete[] m_bufferReceiveTransports[lvl];
  }
  m_reqSendTransports.clear();
  m_bufferSendTransports.clear();
  m_reqReceiveTransports.clear();
  m_bufferReceiveTransports.clear();

}

//***********************************************************************

void Parallel::communicationsTransports(int lvl)
{
  int count(0);
  MPI_Status status;

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Prepation of sendings
      count = -1;
      for (int i = 0; i < m_numberElementsToSendToNeighbour[neighbour]; i++) {
        m_elementsToSend[neighbour][i]->fillBufferTransports(m_bufferSendTransports[lvl][neighbour], count, lvl, neighbour);
      }

      //Sending request
      MPI_Start(m_reqSendTransports[lvl][neighbour]);
      //Receiving request
      MPI_Start(m_reqReceiveTransports[lvl][neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqSendTransports[lvl][neighbour], &status);
      MPI_Wait(m_reqReceiveTransports[lvl][neighbour], &status);

      //Receivings
      count = -1;
      for (int i = 0; i < m_numberElementsToReceiveFromNeighbour[neighbour]; i++) {
        m_elementsToReceive[neighbour][i]->getBufferTransports(m_bufferReceiveTransports[lvl][neighbour], count, lvl);
      }
    }
  }
}

//****************************************************************************
//******************** Methodes pour les variables AMR ***********************
//****************************************************************************

void Parallel::initializePersistentCommunicationsAMR(const int &numberPrimitiveVariables, const int &numberSlopeVariables, const int &numberTransportVariables, const int &dim, const int &lvlMax)
{
  if (Ncpu > 1) {
    m_numberPrimitiveVariables = numberPrimitiveVariables;
    m_numberSlopeVariables = numberSlopeVariables;
    m_numberTransportVariables = numberTransportVariables;
    //Initialization of communications of primitive variables from resolved model
    parallel.initializePersistentCommunicationsPrimitives();
    //Initialization of communications of slopes for second order
    parallel.initializePersistentCommunicationsSlopes();
    //Initialization of communications necessary for additional physics (vectors of dim=3)
    parallel.initializePersistentCommunicationsVector(dim);
    //Initialization of communications of transported variables
    parallel.initializePersistentCommunicationsTransports();
    //Initialization of communications for AMR variables
    parallel.initializePersistentCommunicationsXi();
    parallel.initializePersistentCommunicationsSplit();
    parallel.initializePersistentCommunicationsNumberGhostCells();
    //Initialization of communications for the levels superior to 0
    parallel.initializePersistentCommunicationsLvlAMR(lvlMax);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::initializePersistentCommunicationsLvlAMR(const int &lvlMax)
{
  //Extension of parallel variables to the maximum AMR level. We starts at 1, the level 0 being already initialized
  for (int lvl = 1; lvl <= lvlMax; lvl++) {
    m_bufferSend.push_back(new double*[Ncpu]);
    m_bufferReceive.push_back(new double*[Ncpu]);
    m_bufferSendSlopes.push_back(new double*[Ncpu]);
    m_bufferReceiveSlopes.push_back(new double*[Ncpu]);
    m_bufferSendVector.push_back(new double*[Ncpu]);
    m_bufferReceiveVector.push_back(new double*[Ncpu]);
    m_bufferSendTransports.push_back(new double*[Ncpu]);
    m_bufferReceiveTransports.push_back(new double*[Ncpu]);
    m_bufferSendXi.push_back(new double*[Ncpu]);
    m_bufferReceiveXi.push_back(new double*[Ncpu]);
    m_bufferSendSplit.push_back(new bool*[Ncpu]);
    m_bufferReceiveSplit.push_back(new bool*[Ncpu]);

    m_reqSend.push_back(new MPI_Request*[Ncpu]);
    m_reqReceive.push_back(new MPI_Request*[Ncpu]);
    m_reqSendSlopes.push_back(new MPI_Request*[Ncpu]);
    m_reqReceiveSlopes.push_back(new MPI_Request*[Ncpu]);
    m_reqSendVector.push_back(new MPI_Request*[Ncpu]);
    m_reqReceiveVector.push_back(new MPI_Request*[Ncpu]);
    m_reqSendTransports.push_back(new MPI_Request*[Ncpu]);
    m_reqReceiveTransports.push_back(new MPI_Request*[Ncpu]);
    m_reqSendXi.push_back(new MPI_Request*[Ncpu]);
    m_reqReceiveXi.push_back(new MPI_Request*[Ncpu]);
    m_reqSendSplit.push_back(new MPI_Request*[Ncpu]);
    m_reqReceiveSplit.push_back(new MPI_Request*[Ncpu]);

    for (int i = 0; i < Ncpu; i++) {
      m_bufferSend[lvl][i] = NULL;
      m_bufferReceive[lvl][i] = NULL;
      m_bufferSendSlopes[lvl][i] = NULL;
      m_bufferReceiveSlopes[lvl][i] = NULL;
      m_bufferSendVector[lvl][i] = NULL;
      m_bufferReceiveVector[lvl][i] = NULL;
      m_bufferSendTransports[lvl][i] = NULL;
      m_bufferReceiveTransports[lvl][i] = NULL;
      m_bufferSendXi[lvl][i] = NULL;
      m_bufferReceiveXi[lvl][i] = NULL;
      m_bufferSendSplit[lvl][i] = NULL;
      m_bufferReceiveSplit[lvl][i] = NULL;

      m_reqSend[lvl][i] = NULL;
      m_reqReceive[lvl][i] = NULL;
      m_reqSendSlopes[lvl][i] = NULL;
      m_reqReceiveSlopes[lvl][i] = NULL;
      m_reqSendVector[lvl][i] = NULL;
      m_reqReceiveVector[lvl][i] = NULL;
      m_reqSendTransports[lvl][i] = NULL;
      m_reqReceiveTransports[lvl][i] = NULL;
      m_reqSendXi[lvl][i] = NULL;
      m_reqReceiveXi[lvl][i] = NULL;
      m_reqSendSplit[lvl][i] = NULL;
      m_reqReceiveSplit[lvl][i] = NULL;
    }
  }

  //Initialization of sendings and receivings for the couples of neighboring CPU and for each AMR level
  int numberSend(0);
  int numberReceive(0);

  for (int lvl = 1; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
      if (m_isNeighbour[neighbour]) {
        //Primitive variables
        //-------------------
        //New sending request and its associated buffer
        m_reqSend[lvl][neighbour] = new MPI_Request;
        m_bufferSend[lvl][neighbour] = new double[numberSend];
        MPI_Send_init(m_bufferSend[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSend[lvl][neighbour]);

        //New receiving request and its associated buffer
        m_reqReceive[lvl][neighbour] = new MPI_Request;
        m_bufferReceive[lvl][neighbour] = new double[numberReceive];
        MPI_Recv_init(m_bufferReceive[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceive[lvl][neighbour]);

        //Slope variables
        //---------------
        //New sending request and its associated buffer
        m_reqSendSlopes[lvl][neighbour] = new MPI_Request;
        m_bufferSendSlopes[lvl][neighbour] = new double[numberSend];
        MPI_Send_init(m_bufferSendSlopes[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendSlopes[lvl][neighbour]);

        //New receiving request and its associated buffer
        m_reqReceiveSlopes[lvl][neighbour] = new MPI_Request;
        m_bufferReceiveSlopes[lvl][neighbour] = new double[numberReceive];
        MPI_Recv_init(m_bufferReceiveSlopes[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveSlopes[lvl][neighbour]);

        //Vector variables
        //----------------
        //New sending request and its associated buffer
        m_reqSendVector[lvl][neighbour] = new MPI_Request;
        m_bufferSendVector[lvl][neighbour] = new double[numberSend];
        MPI_Send_init(m_bufferSendVector[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendVector[lvl][neighbour]);

        //New receiving request and its associated buffer
        m_reqReceiveVector[lvl][neighbour] = new MPI_Request;
        m_bufferReceiveVector[lvl][neighbour] = new double[numberReceive];
        MPI_Recv_init(m_bufferReceiveVector[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveVector[lvl][neighbour]);

        //Transported variables
        //---------------------
        //New sending request and its associated buffer
        m_reqSendTransports[lvl][neighbour] = new MPI_Request;
        m_bufferSendTransports[lvl][neighbour] = new double[numberSend];
        MPI_Send_init(m_bufferSendTransports[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendTransports[lvl][neighbour]);

        //New receiving request and its associated buffer
        m_reqReceiveTransports[lvl][neighbour] = new MPI_Request;
        m_bufferReceiveTransports[lvl][neighbour] = new double[numberReceive];
        MPI_Recv_init(m_bufferReceiveTransports[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveTransports[lvl][neighbour]);

        //Xi variable
        //-----------
        //New sending request and its associated buffer
        m_reqSendXi[lvl][neighbour] = new MPI_Request;
        m_bufferSendXi[lvl][neighbour] = new double[numberSend];
        MPI_Send_init(m_bufferSendXi[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendXi[lvl][neighbour]);

        //New receiving request and its associated buffer
        m_reqReceiveXi[lvl][neighbour] = new MPI_Request;
        m_bufferReceiveXi[lvl][neighbour] = new double[numberReceive];
        MPI_Recv_init(m_bufferReceiveXi[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveXi[lvl][neighbour]);

        //Split variable
        //--------------
        //New sending request and its associated buffer
        m_reqSendSplit[lvl][neighbour] = new MPI_Request;
        m_bufferSendSplit[lvl][neighbour] = new bool[numberSend];
        MPI_Send_init(m_bufferSendSplit[lvl][neighbour], numberSend, MPI_C_BOOL, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendSplit[lvl][neighbour]);

        //New receiving request and its associated buffer
        m_reqReceiveSplit[lvl][neighbour] = new MPI_Request;
        m_bufferReceiveSplit[lvl][neighbour] = new bool[numberReceive];
        MPI_Recv_init(m_bufferReceiveSplit[lvl][neighbour], numberReceive, MPI_C_BOOL, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveSplit[lvl][neighbour]);
      }
    }
  }
}

//***********************************************************************

void Parallel::clearRequestsAndBuffers(int lvl)
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_reqSend[lvl][neighbour] != NULL) {
      MPI_Request_free(m_reqSend[lvl][neighbour]);
      MPI_Request_free(m_reqReceive[lvl][neighbour]);
      MPI_Request_free(m_reqSendSlopes[lvl][neighbour]);
      MPI_Request_free(m_reqReceiveSlopes[lvl][neighbour]);
      MPI_Request_free(m_reqSendVector[lvl][neighbour]);
      MPI_Request_free(m_reqReceiveVector[lvl][neighbour]);
      MPI_Request_free(m_reqSendTransports[lvl][neighbour]);
      MPI_Request_free(m_reqReceiveTransports[lvl][neighbour]);
      MPI_Request_free(m_reqSendXi[lvl][neighbour]);
      MPI_Request_free(m_reqReceiveXi[lvl][neighbour]);
      MPI_Request_free(m_reqSendSplit[lvl][neighbour]);
      MPI_Request_free(m_reqReceiveSplit[lvl][neighbour]);

      delete m_reqSend[lvl][neighbour];
      delete m_reqReceive[lvl][neighbour];
      delete m_reqSendSlopes[lvl][neighbour];
      delete m_reqReceiveSlopes[lvl][neighbour];
      delete m_reqSendVector[lvl][neighbour];
      delete m_reqReceiveVector[lvl][neighbour];
      delete m_reqSendTransports[lvl][neighbour];
      delete m_reqReceiveTransports[lvl][neighbour];
      delete m_reqSendXi[lvl][neighbour];
      delete m_reqReceiveXi[lvl][neighbour];
      delete m_reqSendSplit[lvl][neighbour];
      delete m_reqReceiveSplit[lvl][neighbour];

      m_reqSend[lvl][neighbour] = NULL;
      m_reqReceive[lvl][neighbour] = NULL;
      m_reqSendSlopes[lvl][neighbour] = NULL;
      m_reqReceiveSlopes[lvl][neighbour] = NULL;
      m_reqSendVector[lvl][neighbour] = NULL;
      m_reqReceiveVector[lvl][neighbour] = NULL;
      m_reqSendTransports[lvl][neighbour] = NULL;
      m_reqReceiveTransports[lvl][neighbour] = NULL;
      m_reqSendXi[lvl][neighbour] = NULL;
      m_reqReceiveXi[lvl][neighbour] = NULL;
      m_reqSendSplit[lvl][neighbour] = NULL;
      m_reqReceiveSplit[lvl][neighbour] = NULL;

      delete[] m_bufferSend[lvl][neighbour];
      delete[] m_bufferReceive[lvl][neighbour];
      delete[] m_bufferSendSlopes[lvl][neighbour];
      delete[] m_bufferReceiveSlopes[lvl][neighbour];
      delete[] m_bufferSendVector[lvl][neighbour];
      delete[] m_bufferReceiveVector[lvl][neighbour];
      delete[] m_bufferSendTransports[lvl][neighbour];
      delete[] m_bufferReceiveTransports[lvl][neighbour];
      delete[] m_bufferSendXi[lvl][neighbour];
      delete[] m_bufferReceiveXi[lvl][neighbour];
      delete[] m_bufferSendSplit[lvl][neighbour];
      delete[] m_bufferReceiveSplit[lvl][neighbour];

      m_bufferSend[lvl][neighbour] = NULL;
      m_bufferReceive[lvl][neighbour] = NULL;
      m_bufferSendSlopes[lvl][neighbour] = NULL;
      m_bufferReceiveSlopes[lvl][neighbour] = NULL;
      m_bufferSendVector[lvl][neighbour] = NULL;
      m_bufferReceiveVector[lvl][neighbour] = NULL;
      m_bufferSendTransports[lvl][neighbour] = NULL;
      m_bufferReceiveTransports[lvl][neighbour] = NULL;
      m_bufferSendXi[lvl][neighbour] = NULL;
      m_bufferReceiveXi[lvl][neighbour] = NULL;
      m_bufferSendSplit[lvl][neighbour] = NULL;
      m_bufferReceiveSplit[lvl][neighbour] = NULL;
    }
  }
}

//***********************************************************************

void Parallel::updatePersistentCommunicationsAMR(const int &dim)
{
  //We first empty the sending and receiving variables of level 0 (from previous domain decomposition)
  this->clearRequestsAndBuffers(0);

  //Initialization of communications of primitive variables from resolved model
  parallel.initializePersistentCommunicationsPrimitives();
  //Initialization of communications of slopes for second order
  parallel.initializePersistentCommunicationsSlopes();
  //Initialization of communications necessary for additional physics (vectors of dim=3)
  parallel.initializePersistentCommunicationsVector(dim);
  //Initialization of communications of transported variables
  parallel.initializePersistentCommunicationsTransports();
  //Initialization of communications for AMR variables
  parallel.initializePersistentCommunicationsXi();
  parallel.initializePersistentCommunicationsSplit();

  MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::updatePersistentCommunicationsLvlAMR(int lvl, const int &dim)
{
  //We first empty the sending and receiving variables of level lvl (from previous domain decomposition)
  this->clearRequestsAndBuffers(lvl);

  //We write the new sending and receiving variables
  int numberSend(0), numberReceive(0);
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Primitive variables
      //-------------------
      numberSend = m_numberPrimitiveVariables*m_bufferNumberElementsToSendToNeighbor[neighbour];
      numberReceive = m_numberPrimitiveVariables*m_bufferNumberElementsToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSend[lvl][neighbour] = new MPI_Request;
      m_bufferSend[lvl][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSend[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSend[lvl][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceive[lvl][neighbour] = new MPI_Request;
      m_bufferReceive[lvl][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceive[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceive[lvl][neighbour]);

      //Slope variables
      //---------------
      numberSend = m_numberSlopeVariables*m_bufferNumberSlopesToSendToNeighbor[neighbour];
      numberReceive = m_numberSlopeVariables*m_bufferNumberSlopesToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSendSlopes[lvl][neighbour] = new MPI_Request;
      m_bufferSendSlopes[lvl][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendSlopes[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendSlopes[lvl][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveSlopes[lvl][neighbour] = new MPI_Request;
      m_bufferReceiveSlopes[lvl][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveSlopes[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveSlopes[lvl][neighbour]);

      //Vector variables
      //----------------
      numberSend = dim*m_bufferNumberElementsToSendToNeighbor[neighbour];
      numberReceive = dim*m_bufferNumberElementsToReceiveFromNeighbour[neighbour];
      //New sending request and its associated buffer
      m_reqSendVector[lvl][neighbour] = new MPI_Request;
      m_bufferSendVector[lvl][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendVector[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendVector[lvl][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveVector[lvl][neighbour] = new MPI_Request;
      m_bufferReceiveVector[lvl][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveVector[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveVector[lvl][neighbour]);

      //Transported variables
      //---------------------
      numberSend = m_numberTransportVariables*m_bufferNumberElementsToSendToNeighbor[neighbour];
      numberReceive = m_numberTransportVariables*m_bufferNumberElementsToReceiveFromNeighbour[neighbour];
      //New sending request and its associated buffer
      m_reqSendTransports[lvl][neighbour] = new MPI_Request;
      m_bufferSendTransports[lvl][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendTransports[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendTransports[lvl][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveTransports[lvl][neighbour] = new MPI_Request;
      m_bufferReceiveTransports[lvl][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveTransports[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveTransports[lvl][neighbour]);

      //Xi variable
      //-----------
      numberSend = m_bufferNumberElementsToSendToNeighbor[neighbour];
      numberReceive = m_bufferNumberElementsToReceiveFromNeighbour[neighbour];
      //New sending request and its associated buffer
      m_reqSendXi[lvl][neighbour] = new MPI_Request;
      m_bufferSendXi[lvl][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendXi[lvl][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendXi[lvl][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveXi[lvl][neighbour] = new MPI_Request;
      m_bufferReceiveXi[lvl][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveXi[lvl][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveXi[lvl][neighbour]);

      //Split variable
      //--------------
      //New sending request and its associated buffer
      m_reqSendSplit[lvl][neighbour] = new MPI_Request;
      m_bufferSendSplit[lvl][neighbour] = new bool[numberSend];
      MPI_Send_init(m_bufferSendSplit[lvl][neighbour], numberSend, MPI_C_BOOL, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendSplit[lvl][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveSplit[lvl][neighbour] = new MPI_Request;
      m_bufferReceiveSplit[lvl][neighbour] = new bool[numberReceive];
      MPI_Recv_init(m_bufferReceiveSplit[lvl][neighbour], numberReceive, MPI_C_BOOL, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveSplit[lvl][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizeAMR(const int &lvlMax)
{
  if (Ncpu > 1) {
    this->finalizePersistentCommunicationsXi(lvlMax);
    this->finalizePersistentCommunicationsSplit(lvlMax);
    this->finalizePersistentCommunicationsNumberGhostCells();
  }
  this->finalize(lvlMax);
  MPI_Barrier(MPI_COMM_WORLD);
}

//***********************************************************************

void Parallel::initializePersistentCommunicationsXi()
{
  int number(1);

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Determination of the number of variables to communicate
      int numberSend = number*m_numberElementsToSendToNeighbour[neighbour];
      int numberReceive = number*m_numberElementsToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSendXi[0][neighbour] = new MPI_Request;
      m_bufferSendXi[0][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendXi[0][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendXi[0][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveXi[0][neighbour] = new MPI_Request;
      m_bufferReceiveXi[0][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveXi[0][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveXi[0][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsXi(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
      if (m_isNeighbour[neighbour]) {
        MPI_Request_free(m_reqSendXi[lvl][neighbour]);
        MPI_Request_free(m_reqReceiveXi[lvl][neighbour]);
        delete m_reqSendXi[lvl][neighbour];
        delete[] m_bufferSendXi[lvl][neighbour];
        delete m_reqReceiveXi[lvl][neighbour];
        delete[] m_bufferReceiveXi[lvl][neighbour];
      }
    }
    delete[] m_reqSendXi[lvl];
    delete[] m_bufferSendXi[lvl];
    delete[] m_reqReceiveXi[lvl];
    delete[] m_bufferReceiveXi[lvl];
  }
  m_reqSendXi.clear();
  m_bufferSendXi.clear();
  m_reqReceiveXi.clear();
  m_bufferReceiveXi.clear();
}

//***********************************************************************

void Parallel::communicationsXi(int lvl)
{
  int count(0);
  MPI_Status status;

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Prepation of sendings
      count = -1;
      for (int i = 0; i < m_numberElementsToSendToNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferSendXi
        m_elementsToSend[neighbour][i]->fillBufferXi(m_bufferSendXi[lvl][neighbour], count, lvl, neighbour);
      }

      //Sending request
      MPI_Start(m_reqSendXi[lvl][neighbour]);
      //Receiving request
      MPI_Start(m_reqReceiveXi[lvl][neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqSendXi[lvl][neighbour], &status);
      MPI_Wait(m_reqReceiveXi[lvl][neighbour], &status);

      //Receivings
      count = -1;
      for (int i = 0; i < m_numberElementsToReceiveFromNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferReceiveXi
        m_elementsToReceive[neighbour][i]->getBufferXi(m_bufferReceiveXi[lvl][neighbour], count, lvl);
      }
    }
  }
}

//***********************************************************************

void Parallel::initializePersistentCommunicationsSplit()
{
  int number(1);

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Determination of the number of variables to communicate
      int numberSend = number*m_numberElementsToSendToNeighbour[neighbour];
      int numberReceive = number*m_numberElementsToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSendSplit[0][neighbour] = new MPI_Request;
      m_bufferSendSplit[0][neighbour] = new bool[numberSend];
      MPI_Send_init(m_bufferSendSplit[0][neighbour], numberSend, MPI_C_BOOL, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendSplit[0][neighbour]);
      
      //New receiving request and its associated buffer
      m_reqReceiveSplit[0][neighbour] = new MPI_Request;
      m_bufferReceiveSplit[0][neighbour] = new bool[numberReceive];
      MPI_Recv_init(m_bufferReceiveSplit[0][neighbour], numberReceive, MPI_C_BOOL, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveSplit[0][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsSplit(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
      if (m_isNeighbour[neighbour]) {
        MPI_Request_free(m_reqSendSplit[lvl][neighbour]);
        MPI_Request_free(m_reqReceiveSplit[lvl][neighbour]);
        delete m_reqSendSplit[lvl][neighbour];
        delete[] m_bufferSendSplit[lvl][neighbour];
        delete m_reqReceiveSplit[lvl][neighbour];
        delete[] m_bufferReceiveSplit[lvl][neighbour];
      }
    }
    delete[] m_reqSendSplit[lvl];
    delete[] m_bufferSendSplit[lvl];
    delete[] m_reqReceiveSplit[lvl];
    delete[] m_bufferReceiveSplit[lvl];
  }
  m_reqSendSplit.clear();
  m_bufferSendSplit.clear();
  m_reqReceiveSplit.clear();
  m_bufferReceiveSplit.clear();
}

//***********************************************************************

void Parallel::communicationsSplit(int lvl)
{
  int count(0);
  MPI_Status status;

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Prepation of sendings
      count = -1;
      for (int i = 0; i < m_numberElementsToSendToNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferSendSplit
        m_elementsToSend[neighbour][i]->fillBufferSplit(m_bufferSendSplit[lvl][neighbour], count, lvl, neighbour);
      }

      //Sending request
      MPI_Start(m_reqSendSplit[lvl][neighbour]);
      //Receiving request
      MPI_Start(m_reqReceiveSplit[lvl][neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqSendSplit[lvl][neighbour], &status);
      MPI_Wait(m_reqReceiveSplit[lvl][neighbour], &status);

      //Receivings
      count = -1;
      for (int i = 0; i < m_numberElementsToReceiveFromNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferReceiveSplit
        m_elementsToReceive[neighbour][i]->getBufferSplit(m_bufferReceiveSplit[lvl][neighbour], count, lvl);
      }
    }
  }
}

//***********************************************************************

void Parallel::initializePersistentCommunicationsNumberGhostCells()
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    //Determination of the number of variables to communicate
    int numberSend = 1;
    int numberReceive = 1;

    //New sending request and its associated buffer
    m_reqNumberElementsToSendToNeighbor[neighbour] = new MPI_Request;
    m_bufferNumberElementsToSendToNeighbor[neighbour] = 0;
    MPI_Send_init(&m_bufferNumberElementsToSendToNeighbor[neighbour], numberSend, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD, m_reqNumberElementsToSendToNeighbor[neighbour]);

    //New receiving request and its associated buffer
    m_reqNumberElementsToReceiveFromNeighbour[neighbour] = new MPI_Request;
    m_bufferNumberElementsToReceiveFromNeighbour[neighbour] = 0;
    MPI_Recv_init(&m_bufferNumberElementsToReceiveFromNeighbour[neighbour], numberReceive, MPI_INT, neighbour, rankCpu, MPI_COMM_WORLD, m_reqNumberElementsToReceiveFromNeighbour[neighbour]);

    //New sending request and its associated buffer
    m_reqNumberSlopesToSendToNeighbor[neighbour] = new MPI_Request;
    m_bufferNumberSlopesToSendToNeighbor[neighbour] = 0;
    MPI_Send_init(&m_bufferNumberSlopesToSendToNeighbor[neighbour], numberSend, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD, m_reqNumberSlopesToSendToNeighbor[neighbour]);

    //New receiving request and its associated buffer
    m_reqNumberSlopesToReceiveFromNeighbour[neighbour] = new MPI_Request;
    m_bufferNumberSlopesToReceiveFromNeighbour[neighbour] = 0;
    MPI_Recv_init(&m_bufferNumberSlopesToReceiveFromNeighbour[neighbour], numberReceive, MPI_INT, neighbour, rankCpu, MPI_COMM_WORLD, m_reqNumberSlopesToReceiveFromNeighbour[neighbour]);
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsNumberGhostCells()
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      MPI_Request_free(m_reqNumberElementsToSendToNeighbor[neighbour]);
      MPI_Request_free(m_reqNumberElementsToReceiveFromNeighbour[neighbour]);
      MPI_Request_free(m_reqNumberSlopesToSendToNeighbor[neighbour]);
      MPI_Request_free(m_reqNumberSlopesToReceiveFromNeighbour[neighbour]);
      delete m_reqNumberElementsToSendToNeighbor[neighbour];
      delete m_reqNumberElementsToReceiveFromNeighbour[neighbour];
      delete m_reqNumberSlopesToSendToNeighbor[neighbour];
      delete m_reqNumberSlopesToReceiveFromNeighbour[neighbour];
    }
  }
  delete[] m_reqNumberElementsToSendToNeighbor;
  delete[] m_reqNumberElementsToReceiveFromNeighbour;
  delete[] m_reqNumberSlopesToSendToNeighbor;
  delete[] m_reqNumberSlopesToReceiveFromNeighbour;
  delete[] m_bufferNumberElementsToSendToNeighbor;
  delete[] m_bufferNumberElementsToReceiveFromNeighbour;
  delete[] m_bufferNumberSlopesToSendToNeighbor;
  delete[] m_bufferNumberSlopesToReceiveFromNeighbour;
}

//***********************************************************************

void Parallel::communicationsNumberGhostCells(int lvl)
{
  MPI_Status status;

  for (int neighbour = 0; neighbour < Ncpu; neighbour++)
  {
    if (m_isNeighbour[neighbour]) {
      //Prepation de l'envoi
      m_bufferNumberElementsToSendToNeighbor[neighbour] = 0;
      m_bufferNumberSlopesToSendToNeighbor[neighbour] = 0;
      for (int i = 0; i < m_numberElementsToSendToNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferNumberElementsToSendToNeighbor and m_bufferNumberSlopesToSendToNeighbor
        m_elementsToSend[neighbour][i]->fillNumberElementsToSendToNeighbour(m_bufferNumberElementsToSendToNeighbor[neighbour], m_bufferNumberSlopesToSendToNeighbor[neighbour], lvl, neighbour, 0);
      }

      //For elements
      //Sending request
      MPI_Start(m_reqNumberElementsToSendToNeighbor[neighbour]);
      //Receiving request
      MPI_Start(m_reqNumberElementsToReceiveFromNeighbour[neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++)
  {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqNumberElementsToSendToNeighbor[neighbour], &status);
      MPI_Wait(m_reqNumberElementsToReceiveFromNeighbour[neighbour], &status);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++)
  {
    if (m_isNeighbour[neighbour]) {
      //For slopes
      //Sending request
      MPI_Start(m_reqNumberSlopesToSendToNeighbor[neighbour]);
      //Receiving request
      MPI_Start(m_reqNumberSlopesToReceiveFromNeighbour[neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++)
  {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqNumberSlopesToSendToNeighbor[neighbour], &status);
      MPI_Wait(m_reqNumberSlopesToReceiveFromNeighbour[neighbour], &status);
    }
  }
}

//***********************************************************************
