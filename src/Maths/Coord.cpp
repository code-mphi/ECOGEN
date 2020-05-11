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

//! \file      Coord.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include <algorithm>
#include <cmath>
#include "Coord.h"
#include <iostream>

//*********************************************************************

Coord::Coord() : m_x(0.), m_y(0.), m_z(0.){}

//*********************************************************************

Coord::Coord(const double &x, const double &y, const double &z) :
m_x(x), m_y(y), m_z(z)
{}

//*********************************************************************

Coord::~Coord(){}

//*********************************************************************

const Coord& Coord::coord() const
{
    return defaultCoord;
}

//*********************************************************************

Coord& Coord::coord()
{
    return defaultCoordNonConst;
}

//*********************************************************************

const Coord Coord::defaultCoord = Coord();

//*********************************************************************

Coord Coord::defaultCoordNonConst = Coord();

//*********************************************************************

void Coord::setXYZ(const double &x, const double & y, const double &z)
{
  m_x = x;
  m_y = y;
  m_z = z;
}

//*********************************************************************

void Coord::setX(const double &x){m_x = x;}

//*********************************************************************

void Coord::setY(const double &y){m_y = y;}

//*********************************************************************

void Coord::setZ(const double &z){m_z = z;}

//*********************************************************************

double Coord::norm() const
{
  return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
}

//*********************************************************************

double Coord::squaredNorm() const
{
  return m_x*m_x + m_y*m_y + m_z*m_z;
}

//*********************************************************************

Coord Coord::abs() const
{
  Coord vec;
  vec.m_x = std::fabs(m_x);
  vec.m_y = std::fabs(m_y);
  vec.m_z = std::fabs(m_z);
  return vec;
}

//*********************************************************************

double Coord::scalar(const Coord &a) const
{
  return m_x*a.m_x + m_y*a.m_y + m_z*a.m_z;
}

//*********************************************************************

Coord Coord::cross(const Coord &a) const
{
  Coord vec;
  vec.m_x = m_y*a.m_z - m_z*a.m_y;
  vec.m_y = m_z*a.m_x - m_x*a.m_z;
  vec.m_z = m_x*a.m_y - m_y*a.m_x;
  return vec;
}

//*********************************************************************

void Coord::localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  Coord vecteurProjete;
  vecteurProjete.m_x = this->scalar(normal);
  vecteurProjete.m_y = this->scalar(tangent);
  vecteurProjete.m_z = this->scalar(binormal);
  *this = vecteurProjete;
}

//*********************************************************************

void Coord::reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  Coord vecteurProjete;
  vecteurProjete.m_x = normal.m_x*m_x + tangent.m_x*m_y + binormal.m_x*m_z;
  vecteurProjete.m_y = normal.m_y*m_x + tangent.m_y*m_y + binormal.m_y*m_z;
  vecteurProjete.m_z = normal.m_z*m_x + tangent.m_z*m_y + binormal.m_z*m_z;
  *this = vecteurProjete;
}

//*********************************************************************

void Coord::setFromSubtractedVectors(const Coord &a, const Coord &b)
{
  m_x = b.m_x - a.m_x;
  m_y = b.m_y - a.m_y;
  m_z = b.m_z - a.m_z;
}

//*********************************************************************

double Coord::scalarProduct(const Coord &v1, const Coord &v2)
{
  return v1.m_x*v2.m_x + v1.m_y*v2.m_y + v1.m_z*v2.m_z;
}

//*********************************************************************

Coord Coord::crossProduct(const Coord &v1, const Coord &v2)
{
  Coord vec;
  vec.m_x = v1.m_y*v2.m_z - v1.m_z*v2.m_y;
  vec.m_y = v1.m_z*v2.m_x - v1.m_x*v2.m_z;
  vec.m_z = v1.m_x*v2.m_y - v1.m_y*v2.m_x;
  return vec;
}

//*********************************************************************

void Coord::changeSign()
{
  m_x = -m_x;
  m_y = -m_y;
  m_z = -m_z;
}

//*********************************************************************

void Coord::normalized()
{
  *this = *this / this->norm();
}

//*********************************************************************

double Coord::determinant(const Coord &v1, const Coord &v2, const Coord &v3)
{
  double det(0.); 
  det =  (v1.getX()*v2.getY()*v3.getZ()) + (v1.getY()*v2.getZ()*v3.getX()) + (v1.getZ()*v2.getX()*v3.getY());
  det -= (v3.getX()*v2.getY()*v1.getZ()) + (v3.getY()*v2.getZ()*v1.getX()) + (v3.getZ()*v2.getX()*v1.getY());
  return det;
}

//*********************************************************************

double Coord::cos(const Coord &v1, const Coord &v2)
{
  double cos(0.);
  double nv1(v1.norm()), nv2(v2.norm());
  if (nv1 > 1e-6 && nv2 > 1e-6) { cos = v1.scalar(v2) / (nv1*nv2); }
  return cos;
}

//*********************************************************************

Coord Coord::sin(const Coord &v1, const Coord &v2)
{
  Coord sin;
  double nv1(v1.norm()), nv2(v2.norm());
  if (nv1 > 1e-6 && nv2 > 1e-6) { sin = v1.cross(v2) / (nv1*nv2); }
  return sin;
}

//*********************************************************************

void Coord::printInfo() const
{
  std::cout << " X = " << m_x << " Y = " << m_y << " Z = " << m_z << std::endl;
}

//*******************SURCHARGES OPERATEURS*****************************

Coord& Coord::operator=(const double &scalar)
{
  m_x = scalar;
  m_y = scalar;
  m_z = scalar;
  return *this;
}

//*********************************************************************

Coord& Coord::operator*= (const double &scalar)
{
  m_x *= scalar;
  m_y *= scalar;
  m_z *= scalar;
  return *this;
}

//*********************************************************************

Coord& Coord::operator/= (const double &scalar)
{
  m_x /= scalar;
  m_y /= scalar;
  m_z /= scalar;
  return *this;
}

//*********************************************************************

Coord Coord::operator* (const double &scalar) const
{
  Coord copie(*this);
  copie *= scalar;
  return copie;
}


//*********************************************************************

Coord Coord::operator/ (const double &scalar) const
{
  Coord copie(*this);
  copie /= scalar;
  return copie;
}

//*********************************************************************


Coord& Coord::operator+= (const Coord& a)
{
  m_x += a.m_x;
  m_y += a.m_y;
  m_z += a.m_z;
  return *this;
}

//*********************************************************************

Coord& Coord::operator-= (const Coord &a)
{
  m_x -= a.m_x;
  m_y -= a.m_y;
  m_z -= a.m_z;
  return *this;
}

//*********************************************************************
//Surcharge operateur externe a la classe car prends deux arguments

Coord operator* (const double &scalar, const Coord &a)
{
  Coord copie(a);
  copie *= scalar;
  return copie;
}

Coord operator+ (const Coord &a, const Coord &b)
{
  Coord copie(a);
  copie += b;
  return copie;
}

Coord operator- (const Coord &a, const Coord &b)
{
  Coord copie(a);
  copie -= b;
  return copie;
}
