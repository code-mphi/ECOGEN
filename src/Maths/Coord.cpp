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

#include "Coord.h"

Coord coordBuff;
Coord velocity;
Coord vFaceToElt;

//*********************************************************************

Coord::Coord() : m_x(0.), m_y(0.), m_z(0.){}

//*********************************************************************

Coord::Coord(const double& x, const double& y, const double& z) :
m_x(x), m_y(y), m_z(z)
{}

//*********************************************************************

Coord::~Coord(){}

//*********************************************************************

const Coord Coord::defaultCoord = Coord();

//*********************************************************************

Coord Coord::defaultCoordNonConst = Coord();

//*********************************************************************

void Coord::setXYZ(const double& x, const double&  y, const double& z)
{
  m_x = x;
  m_y = y;
  m_z = z;
}

//*********************************************************************

void Coord::setX(const double& x){m_x = x;}

//*********************************************************************

void Coord::setY(const double& y){m_y = y;}

//*********************************************************************

void Coord::setZ(const double& z){m_z = z;}

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
  coordBuff.m_x = std::fabs(m_x);
  coordBuff.m_y = std::fabs(m_y);
  coordBuff.m_z = std::fabs(m_z);
  return coordBuff;
}

//*********************************************************************

double Coord::scalar(const Coord& a) const
{
  return m_x*a.m_x + m_y*a.m_y + m_z*a.m_z;
}

//*********************************************************************

Coord Coord::scalar(const Tensor& t) const
{
  coordBuff.m_x = m_x * t.getXX() + m_y * t.getYX() + m_z * t.getZX();
  coordBuff.m_y = m_x * t.getXY() + m_y * t.getYY() + m_z * t.getZY();
  coordBuff.m_z = m_x * t.getXZ() + m_y * t.getYZ() + m_z * t.getZZ();
  return coordBuff;
}

//*********************************************************************

Coord Coord::cross(const Coord& a) const
{
  coordBuff.m_x = m_y*a.m_z - m_z*a.m_y;
  coordBuff.m_y = m_z*a.m_x - m_x*a.m_z;
  coordBuff.m_z = m_x*a.m_y - m_y*a.m_x;
  return coordBuff;
}

//*********************************************************************

void Coord::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  //coordBuff here is the projected vector
  coordBuff.m_x = this->scalar(normal);
  coordBuff.m_y = this->scalar(tangent);
  coordBuff.m_z = this->scalar(binormal);
  *this = coordBuff;
}

//*********************************************************************

void Coord::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  //coordBuff here is the projected vector
  coordBuff.m_x = normal.m_x*m_x + tangent.m_x*m_y + binormal.m_x*m_z;
  coordBuff.m_y = normal.m_y*m_x + tangent.m_y*m_y + binormal.m_y*m_z;
  coordBuff.m_z = normal.m_z*m_x + tangent.m_z*m_y + binormal.m_z*m_z;
  *this = coordBuff;
}

//*********************************************************************

void Coord::setFromSubtractedVectors(const Coord& a, const Coord& b)
{
  m_x = b.m_x - a.m_x;
  m_y = b.m_y - a.m_y;
  m_z = b.m_z - a.m_z;
}

//*********************************************************************

double Coord::scalarProduct(const Coord& v1, const Coord& v2)
{
  return v1.m_x*v2.m_x + v1.m_y*v2.m_y + v1.m_z*v2.m_z;
}

//*********************************************************************

Coord Coord::crossProduct(const Coord& v1, const Coord& v2)
{
  coordBuff.m_x = v1.m_y*v2.m_z - v1.m_z*v2.m_y;
  coordBuff.m_y = v1.m_z*v2.m_x - v1.m_x*v2.m_z;
  coordBuff.m_z = v1.m_x*v2.m_y - v1.m_y*v2.m_x;
  return coordBuff;
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

double Coord::determinant(const Coord& v1, const Coord& v2, const Coord& v3)
{
  double det(0.); 
  det =  (v1.m_x*v2.m_y*v3.m_z) + (v1.m_y*v2.m_z*v3.m_x) + (v1.m_z*v2.m_x*v3.m_y);
  det -= (v3.m_x*v2.m_y*v1.m_z) + (v3.m_y*v2.m_z*v1.m_x) + (v3.m_z*v2.m_x*v1.m_y);
  return det;
}

//*********************************************************************

double Coord::cos(const Coord& v1, const Coord& v2)
{
  double cos(0.);
  double nv1(v1.norm()), nv2(v2.norm());
  if (nv1 > 1e-6 && nv2 > 1e-6) { cos = v1.scalar(v2) / (nv1*nv2); }
  return cos;
}

//*********************************************************************

Coord Coord::sin(const Coord& v1, const Coord& v2)
{
  Coord sin;
  double nv1(v1.norm()), nv2(v2.norm());
  if (nv1 > 1e-6 && nv2 > 1e-6) { sin = v1.cross(v2) / (nv1*nv2); }
  return sin;
}

//*********************************************************************

void Coord::buildRelativeVelForRiemannMRF(const Coord &omega, const Coord& normal, const Coord& tangent, const Coord& binormal, const Coord& position)
{
  // Absolute velocity projection into the global frame
  this->reverseProjection(normal, tangent, binormal);
  // Define relative velocity
  *this = *this - omega.cross(position);
  // Reprojection of relative velocity in the face frame
  this->localProjection(normal, tangent, binormal);
}

//***********************************************************************

void Coord::printInfo() const
{
  std::cout << " X = " << m_x << " Y = " << m_y << " Z = " << m_z << std::endl;
}

//*******************SURCHARGES OPERATEURS*****************************

Coord& Coord::operator=(const double& scalar)
{
  m_x = scalar;
  m_y = scalar;
  m_z = scalar;
  return *this;
}

//*********************************************************************

Coord& Coord::operator+= (const double& scalar)
{
  m_x += scalar;
  m_y += scalar;
  m_z += scalar;
  return *this;
}

//*********************************************************************

Coord& Coord::operator-= (const double& scalar)
{
  m_x -= scalar;
  m_y -= scalar;
  m_z -= scalar;
  return *this;
}

//*********************************************************************

Coord& Coord::operator*= (const double& scalar)
{
  m_x *= scalar;
  m_y *= scalar;
  m_z *= scalar;
  return *this;
}

//*********************************************************************

Coord& Coord::operator/= (const double& scalar)
{
  m_x /= scalar;
  m_y /= scalar;
  m_z /= scalar;
  return *this;
}

//*********************************************************************

Coord Coord::operator* (const double& scalar) const
{
  Coord copie(*this);
  copie *= scalar;
  return copie;
}

//*********************************************************************

Coord Coord::operator/ (const double& scalar) const
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

Coord& Coord::operator-= (const Coord& a)
{
  m_x -= a.m_x;
  m_y -= a.m_y;
  m_z -= a.m_z;
  return *this;
}

//*********************************************************************
//Surcharge operateur externe a la classe car prends deux arguments

Coord operator* (const double& scalar, const Coord& a)
{
  Coord copy(a);
  copy *= scalar;
  return copy;
}

//*********************************************************************

Coord operator+ (const Coord& a, const Coord& b)
{
  Coord copy(a);
  copy += b;
  return copy;
}

//*********************************************************************

Coord operator- (const Coord& a, const Coord& b)
{
  Coord copy(a);
  copy -= b;
  return copy;
}

//*********************************************************************
