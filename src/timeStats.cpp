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

//! \file      timeStats.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "timeStats.h"
#include <iostream>

//***********************************************************************

timeStats::timeStats(){}

//***********************************************************************

timeStats::~timeStats(){}

//***********************************************************************

void timeStats::initialize()
{
  m_InitialTime = clock();
  m_computationTime = 0;
  m_AMRTime = 0;
  m_communicationTime = 0;
}

//***********************************************************************

void timeStats::updateComputationTime()
{
  m_computationTime += (clock() - m_InitialTime);
  m_InitialTime = clock();
}

//***********************************************************************

void timeStats::startAMRTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_AMRRefTime = clock();
}

//***********************************************************************

void timeStats::endAMRTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_AMRTime += (clock() - m_AMRRefTime);
}

//***********************************************************************

void timeStats::startCommunicationTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_communicationRefTime = clock();
}

//***********************************************************************

void timeStats::endCommunicationTime()
{
  MPI_Barrier(MPI_COMM_WORLD);
  m_communicationTime += (clock() - m_communicationRefTime);
}

//***********************************************************************

void timeStats::setCompTime(const clock_t &compTime, const clock_t &AMRTime, const clock_t &comTime)
{
  m_computationTime = compTime;
  m_AMRTime = AMRTime;
  m_communicationTime = comTime;
}

//***********************************************************************

void timeStats::printScreenStats(const int &numTest) const
{
  printScreenTime(m_computationTime, "Elapsed time", numTest);
  printScreenTime(m_AMRTime, "AMR time", numTest);
  printScreenTime(m_communicationTime, "Communication time", numTest);

  //Estimation temps restant
  //A faire...
  std::cout << "T" << numTest << " | -------------------------------------------" << std::endl;
}

//***********************************************************************

void timeStats::printScreenTime(const clock_t &time, std::string chaine, const int &numTest) const
{
  //Managing string size
  std::string timeName(" |     " + chaine.substr(0,18));
  for (unsigned int i = 0; i < 19 - chaine.size(); i++) {
    timeName += " ";
  }
  timeName += " = ";

  //printing time
  double convDouble = static_cast<double>(time) / CLOCKS_PER_SEC;
  int convTime = static_cast<int>(convDouble);
  int seconde(convTime);
  if (seconde < 60)
  {
    std::cout << "T" << numTest << timeName << convDouble << " s " << std::endl;
  }
  else
  {
    int minute(seconde / 60);
    seconde = seconde % 60;
    if (minute <60)
    {
      std::cout << "T" << numTest << timeName << minute << " min " << seconde << " s " << std::endl;
    }
    else
    {
      int heure(minute / 60);
      minute = minute % 60;
      std::cout << "T" << numTest << timeName << heure << " h " << minute << " min " << seconde << " s " << std::endl;
    }
  }
}

//***********************************************************************