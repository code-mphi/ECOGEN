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

//! \file      stretchZone.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      September 06 2018

#include "stretchZone.h"
#include "../Errors.h"
#include "../Parallel/Parallel.h"

//***********************************************************************

stretchZone::stretchZone() {}

//***********************************************************************

stretchZone::stretchZone(double startAt, double endAt, double factor, int numberCells) :
  m_startAt(startAt), m_endAt(endAt), m_factor(factor), m_numberCells(numberCells)
{}

//***********************************************************************

stretchZone::~stretchZone() {}

//***********************************************************************

int stretchZone::stretching(std::vector<double> &dX, std::vector<double> &posX)
{
  double dX0((m_endAt - m_startAt) / m_numberCells);
  if (std::fabs(m_factor - 1.) < 1e-6) {
    for (int i = 0; i < m_numberCells; i++) {
      dX.push_back(dX0);
      posX.push_back(m_startAt + (i + 0.5)*dX0);
    }
  }
  else {
    dX0 = (m_endAt - m_startAt)*(1. - m_factor) / (1. - std::pow(m_factor, m_numberCells));
    dX.push_back(dX0);
    posX.push_back(m_startAt + 0.5*dX0);
    int indicePos(posX.size());
    for (int i = 1; i < m_numberCells; i++) {
      dX.push_back(dX0*std::pow(m_factor, i));
      posX.push_back(posX[indicePos - 1] + 0.5*(dX[indicePos - 1] + dX[indicePos]));
      indicePos++;
    }
  }
  return m_numberCells;
}

//***********************************************************************

int stretchZone::verifyStretching(std::vector<stretchZone> &tabStretch, const double l, std::string fileName)
{
  if (tabStretch.size() == 0) return 0;
  try {
    if (tabStretch[0].m_startAt != 0.) throw ErrorXMLStretching(fileName, __FILE__, __LINE__);
    if (tabStretch[tabStretch.size() - 1].m_endAt != l) throw ErrorXMLStretching(fileName, __FILE__, __LINE__);
    for (unsigned int i = 1; i < tabStretch.size(); i++){
      if (tabStretch[i].m_startAt != tabStretch[i-1].m_endAt) throw ErrorXMLStretching(fileName, __FILE__, __LINE__);
      //Verifying size factor between two neighbouring cells
      if (rankCpu == 0) {
        double dX0L((tabStretch[i - 1].m_endAt - tabStretch[i - 1].m_startAt)*(1. - tabStretch[i - 1].m_factor) / (1. - std::pow(tabStretch[i - 1].m_factor, tabStretch[i - 1].m_numberCells)));
        double dXnL(dX0L*std::pow(tabStretch[i - 1].m_factor, tabStretch[i - 1].m_numberCells-1));
        double dX0R((tabStretch[i].m_endAt - tabStretch[i].m_startAt)*(1. - tabStretch[i].m_factor) / (1. - std::pow(tabStretch[i].m_factor, tabStretch[i].m_numberCells)));
        if (dXnL / dX0R > 1.1 || dXnL / dX0R < 0.9) std::cout << "WARNING: cell factor is " << dXnL / dX0R << " between zones " << i - 1 << " and " << i << std::endl;
      }
    }
  }
  catch (ErrorXML &) { throw; } // Renvoi au niveau suivant
  return 0;
}

//***********************************************************************