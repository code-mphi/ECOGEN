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

#include "timeStats.h"

//***********************************************************************

timeStats::timeStats() {}

//***********************************************************************

timeStats::~timeStats() {}

//***********************************************************************

void timeStats::initialize()
{
  m_InitialTime = MPI_Wtime();

  m_computationTime   = 0;
  m_AMRTime           = 0;
  m_communicationTime = 0;
  m_IOTime            = 0;
}

//***********************************************************************

void timeStats::updateComputationTime()
{
  m_computationTime += (MPI_Wtime() - m_InitialTime);
  m_InitialTime      = MPI_Wtime();
}

//***********************************************************************

void timeStats::startAMRTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_AMRRefTime = MPI_Wtime();
}

//***********************************************************************

void timeStats::endAMRTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_AMRTime += (MPI_Wtime() - m_AMRRefTime);
}

//***********************************************************************

void timeStats::startCommunicationTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_communicationRefTime = MPI_Wtime();
}

//***********************************************************************

void timeStats::endCommunicationTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_communicationTime += (MPI_Wtime() - m_communicationRefTime);
}

//***********************************************************************

void timeStats::startIOTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_IORefTime = MPI_Wtime();
}

//***********************************************************************

void timeStats::endIOTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_IOTime += (MPI_Wtime() - m_IORefTime);
}

//***********************************************************************

void timeStats::setCompTime(const double& compTime, const double& AMRTime, const double& comTime, const double& IOTime)
{
  m_computationTime   = compTime;
  m_AMRTime           = AMRTime;
  m_communicationTime = comTime;
  m_IOTime            = IOTime;
}

//***********************************************************************

void timeStats::printScreenStats(const int& numTest) const
{
  printScreenTime(m_computationTime, "Elapsed time", numTest);
  printScreenTime(m_AMRTime, "AMR time", numTest);
  printScreenTime(m_communicationTime, "Communication time", numTest);
  printScreenTime(m_IOTime, "Output time", numTest);

  //Estimation ot remaining time
  //To do...
  std::cout << "T" << numTest << " | -------------------------------------------" << std::endl;
}

//***********************************************************************

void timeStats::printScreenTime(const double& time, std::string chaine, const int& numTest) const
{
  //Managing string size
  std::string timeName(" |     " + chaine.substr(0, 18));
  for (unsigned int i = 0; i < 19 - chaine.size(); i++) {
    timeName += " ";
  }
  timeName += " = ";

  //printing time
  if (time < 60) {
    std::cout << "T" << numTest << timeName << time << " s " << std::endl;
  }
  else {
    int minute(static_cast<int>(std::trunc(time)) / 60);
    double second = time - minute * 60;
    if (minute < 60) {
      std::cout << "T" << numTest << timeName << minute << " min " << second << " s " << std::endl;
    }
    else {
      int hour(minute / 60);
      minute = minute % 60;
      std::cout << "T" << numTest << timeName << hour << " h " << minute << " min " << second << " s " << std::endl;
    }
  }
}

//***********************************************************************
