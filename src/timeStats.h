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

#ifndef TIMESTATS_H
#define TIMESTATS_H

#include <iostream>
#include <string>
#include <mpi.h>
#include <cmath>

class timeStats
{
  public:
    timeStats();
    ~timeStats();

    void initialize();
    void updateComputationTime();

    void startAMRTime();
    void endAMRTime();

    void startCommunicationTime();
    void endCommunicationTime();

    void startIOTime();
    void endIOTime();

    void setCompTime(const double& compTime, const double& AMRTime, const double& comTime, const double& IOTime);
    double getComputationTime() const { return m_computationTime; };
    double getAMRTime() const { return m_AMRTime; };
    double getCommunicationTime() const { return m_communicationTime; };
    double getIOTime() const { return m_IOTime; };
    void printScreenStats(const int& numTest) const;
    void printScreenTime(const double& time, std::string chaine, const int& numTest) const;

  private:
    //Time analysis - Attributes are stored in seconds
    // clang-format off
    double m_InitialTime;
    double m_computationTime;   //!<Computational time

    double m_AMRRefTime;
    double m_AMRTime;           //!<AMR time among computational time

    double m_communicationRefTime;
    double m_communicationTime; //!<Communication time among computational time

    double m_IORefTime;
    double m_IOTime;            //!<IO time among compuational time

    // clang-format on
};

#endif // TIMESTATS_H
