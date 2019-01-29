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
//! \author    F. Petitpas
//! \version   1.0
//! \date      September 07 2018

#include "timeStats.h"
#include <iostream>

using namespace std;

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
  m_AMRRefTime = clock();
}

//***********************************************************************

void timeStats::endAMRTime()
{
  m_AMRTime += (clock() - m_AMRRefTime);
}

//***********************************************************************

void timeStats::setCompTime(const clock_t &time) { m_computationTime = time; }

//***********************************************************************

clock_t timeStats::getComputationTime() const { return m_computationTime; }

//***********************************************************************

void timeStats::printScreenStats(const int &numTest) const
{
  printScreenTime(m_computationTime, "Elapsed time", numTest);
  printScreenTime(m_AMRTime, "AMR time", numTest);

  //Estimation temps restant
  //A faire...
  cout << "T" << numTest << " | ------------------------------------------" << endl;
}

//***********************************************************************

void timeStats::printScreenTime(const clock_t &time, string chaine, const int &numTest) const
{
  //Managing string size
  string timeName(" |     " + chaine.substr(0,18));
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
    cout << "T" << numTest << timeName << convDouble << " s " << endl;
  }
  else
  {
    int minute(seconde / 60);
    seconde = seconde % 60;
    if (minute <60)
    {
      cout << "T" << numTest << timeName << minute << " min " << seconde << " s " << endl;
    }
    else
    {
      int heure(minute / 60);
      minute = minute % 60;
      cout << "T" << numTest << timeName << heure << " h " << minute << " min " << seconde << " s " << endl;
    }
  }
}

//***********************************************************************