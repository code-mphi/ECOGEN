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

#ifndef STRETCHZONE_H
#define STRETCHZONE_H

//! \file      stretchZone.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      September 06 2018

#include <vector>
#include <string>

//! \class     stretchZone
//! \brief     managing a stretched zone in Cartesian meshes

class stretchZone
{
  public:
    //! \brief     basic stretched zone constructor
    stretchZone();
    //! \brief     improved stretched zone constructor
    stretchZone(double startAt, double endtAt, double factor, int numberCells);
    virtual ~stretchZone();

    int stretching(std::vector<double> &dX, std::vector<double> &posX);

    static int verifyStretching(std::vector<stretchZone> &tabStretch, const double l, std::string fileName = "");

  private:
    double m_startAt;   //!< zone starting position along corresponding axis
    double m_endAt;     //!< zone ending position along corresponding axis
    double m_factor;    //!< stretching factor (1. is no stretching)
    int m_numberCells;  //!< number of cells in the stretched zone
};

#endif // STRETCHZONE_H