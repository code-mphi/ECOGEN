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

#ifndef GOLINE_H
#define GOLINE_H

//! \file      GOLine.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      January 5 2018

#include "GeometricObject.h"

//! \class     GOLine
//! \brief     Class for a line geometric object
class GOLine : public GeometricObject
{
public:
  GOLine();
  GOLine(const Coord &vertex, const Coord &vecDir);
  virtual ~GOLine();

  virtual double distancePoint(const Coord &vertex) const;
  virtual Coord projectionPoint(const Coord &vertex) const;

private:
  Coord m_point;           //! Point from the line
  Coord m_vecDir;          //! Director vector of the line

};

#endif //GOLINE_H
