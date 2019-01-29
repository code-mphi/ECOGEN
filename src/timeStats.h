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

#ifndef TIMESTATS_H
#define TIMESTATS_H

//! \file      timeStats.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      September 07 2018

#include <ctime>
#include <string>

class timeStats
{
  public:
    timeStats();
    virtual ~timeStats();

    void initialize();
    void updateComputationTime();

    void startAMRTime();
    void endAMRTime();

    void setCompTime(const clock_t &time);
    clock_t getComputationTime() const;
    void printScreenStats(const int &numTest) const;
    void printScreenTime(const clock_t &time, std::string chaine, const int &numTest) const;

  private:
    //Time analysis - Attributes are stored in miliseconds (to be divided by CLOCKS_PER_SEC)
    clock_t m_InitialTime;
    clock_t m_computationTime;            //!<Computation time
    
    clock_t m_AMRRefTime;
    clock_t m_AMRTime;                    //!<AMR additional time

};

#endif // TIMESTATS_H