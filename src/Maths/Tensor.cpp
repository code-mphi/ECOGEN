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

#include "Tensor.h"

Tensor projectedTensor;

Tensor::Tensor() : m_xx(0.), m_xy(0.), m_xz(0.), m_yx(0.), m_yy(0.), m_yz(0.), m_zx(0.), m_zy(0.), m_zz(0.)
{}

Tensor::Tensor(const double& xx, const double& xy, const double& xz, const double& yx, const double& yy, const double& yz, const double& zx, const double& zy, const double& zz)
  : m_xx(xx), m_xy(xy), m_xz(xz), m_yx(yx), m_yy(yy), m_yz(yz), m_zx(zx), m_zy(zy), m_zz(zz)
{}

Tensor::Tensor(const Coord& x, const Coord& y, const Coord& z)
{
  m_xx = x.getX();
  m_xy = x.getY();
  m_xz = x.getZ();

  m_yx = y.getX();
  m_yy = y.getY();
  m_yz = y.getZ();

  m_zx = z.getX();
  m_zy = z.getY();
  m_zz = z.getZ();
}

Tensor::~Tensor() 
{}

void Tensor::scalar(const Coord& a) const
{
  coordBuff.setX(m_xx * a.getX() + m_xy * a.getY() + m_xz * a.getZ());
  coordBuff.setY(m_yx * a.getX() + m_yy * a.getY() + m_yz * a.getZ());
  coordBuff.setZ(m_zx * a.getX() + m_zy * a.getY() + m_zz * a.getZ());
}

void Tensor::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{ 
  this->scalar(normal);
  projectedTensor.m_xx = coordBuff.getX();
  projectedTensor.m_xy = coordBuff.getY();
  projectedTensor.m_xz = coordBuff.getZ();

  this->scalar(tangent);
  projectedTensor.m_yx = coordBuff.getX();
  projectedTensor.m_yy = coordBuff.getY();
  projectedTensor.m_yz = coordBuff.getZ();

  this->scalar(binormal);
  projectedTensor.m_zx = coordBuff.getX();
  projectedTensor.m_zy = coordBuff.getY();
  projectedTensor.m_zz = coordBuff.getZ();

  *this = projectedTensor;
}

void Tensor::setTensorByLines(const Coord& x, const Coord& y, const Coord& z)
{
  m_xx = x.getX();
  m_xy = x.getY();
  m_xz = x.getZ();

  m_yx = y.getX();
  m_yy = y.getY();
  m_yz = y.getZ();

  m_zx = z.getX();
  m_zy = z.getY();
  m_zz = z.getZ();
}

void Tensor::setTensorByColumns(const Coord& x, const Coord& y, const Coord& z)
{
  m_xx = x.getX();
  m_yx = x.getY();
  m_zx = x.getZ();

  m_xy = y.getX();
  m_yy = y.getY();
  m_zy = y.getZ();

  m_xz = z.getX();
  m_yz = z.getY();
  m_zz = z.getZ();
}

void Tensor::tensorToCoords(Coord& x, Coord& y, Coord& z)
{
  x.setXYZ(m_xx, m_xy, m_xz);
  y.setXYZ(m_yx, m_yy, m_yz);
  z.setXYZ(m_zx, m_zy, m_zz);
}