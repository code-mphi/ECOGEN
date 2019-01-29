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

#ifndef GOPLAN_H
#define GOPLAN_H

//! \file      GOPlan.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      January 5 2018

#include "GeometricObject.h"

//! \class     GOPlan
//! \brief     Class for a plan geometric object
class GOPlan : public GeometricObject
{
public:
  GOPlan();
  GOPlan(const Coord &vertex, const Coord &normal);
  virtual ~GOPlan();

  virtual double distancePoint(const Coord &vertex) const;
  virtual Coord projectionPoint(const Coord &vertex) const;

private:
  Coord m_point;          //! Point from the plan
  Coord m_normal;         //! Normal vector of the plan
  Coord m_tangent;        //! Tangent vector of the plan
  Coord m_binormal;       //! Binormal vector of the plan

  //! \brief     Create the plan base
  void createBase();
};

#endif //GOPLAN_H