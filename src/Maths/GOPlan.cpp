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

//! \file      GOPlan.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      January 5 2018

#include "GOPlan.h"

//***********************************************************************

GOPlan::GOPlan(){}

//***********************************************************************

GOPlan::GOPlan(const Coord &vertex, const Coord &normal) :
  GeometricObject(PLAN), m_point(vertex), m_normal(normal)
{
  if (normal.norm() < 1e-6) { 
    Errors::errorMessage("GOPlan::GOPlan impossible to create plan, normal vector null"); 
  }
  else {
    m_normal.normalized();
  }
  createBase();
}


//***********************************************************************

GOPlan::~GOPlan(){}

//***********************************************************************

double GOPlan::distancePoint(const Coord &vertex) const
{
  Coord vec; vec.setFromSubtractedVectors(m_point, vertex);
  return  std::fabs(vec.scalar(m_normal));
}

//***********************************************************************

Coord GOPlan::projectionPoint(const Coord &vertex) const
{
  Coord projection;
  Coord vec; vec.setFromSubtractedVectors(m_point, vertex);
  projection.setXYZ(m_tangent.scalar(vec), m_binormal.scalar(vec), 0.);
  return projection;
}

//***********************************************************************

void GOPlan::createBase()
{
  Coord M, N;
  M = m_point;

  //Determination of a tangent
  if (std::fabs(m_normal.getZ()) >= 1e-6) {
    N.setX(0); N.setY(0);
    N.setZ((M.getX()*m_normal.getX() + M.getY()*m_normal.getY()) / m_normal.getZ() + M.getZ());
  }
  else if (std::fabs(m_normal.getY()) >= 1e-6) {
    N.setX(0); N.setZ(0);
    N.setY((M.getX()*m_normal.getX() + M.getZ()*m_normal.getZ()) / m_normal.getY() + M.getY());
  }
  else if (std::fabs(m_normal.getX()) >= 1e-6) {
    N.setY(0); N.setZ(0);
    N.setX((M.getY()*m_normal.getY() + M.getZ()*m_normal.getZ()) / m_normal.getX() + M.getX());
  }
  else {
    Errors::errorMessage("GOPlan::createBase impossible, normal vector is problematic");
  }
  m_tangent.setFromSubtractedVectors(M, N);
  m_tangent.normalized();
  m_binormal = Coord::crossProduct(m_normal, m_tangent);
}

//***********************************************************************