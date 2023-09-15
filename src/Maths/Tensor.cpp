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

Tensor tensorBuff;
Tensor tensorIdentity;
Tensor tensorCobase;
Tensor tensorNonConsCobase;
Tensor tensorF;
Tensor tensorG;
Tensor tensorG2;
Tensor tensorA;
Tensor tensorEigenvalues;
Tensor tensorP;
Tensor tensorPinverse;
Tensor tensorD;

Tensor::Tensor()
{
  for (int i = 0; i < 9; i++) {
    m_array[i] = 0.;
  }
}

//*********************************************************************

Tensor::Tensor(const Tensor& tensor)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] = tensor.m_array[i];
  }
}

//*********************************************************************

Tensor::Tensor(const double& xx, const double& xy, const double& xz, const double& yx, const double& yy, const double& yz, const double& zx, const double& zy, const double& zz)
{
  m_array[XX] = xx;
  m_array[XY] = xy;
  m_array[XZ] = xz;

  m_array[YX] = yx;
  m_array[YY] = yy;
  m_array[YZ] = yz;

  m_array[ZX] = zx;
  m_array[ZY] = zy;
  m_array[ZZ] = zz;
}

//*********************************************************************

Tensor::Tensor(const Coord& x, const Coord& y, const Coord& z)
{
  m_array[XX] = x.getX();
  m_array[XY] = x.getY();
  m_array[XZ] = x.getZ();

  m_array[YX] = y.getX();
  m_array[YY] = y.getY();
  m_array[YZ] = y.getZ();

  m_array[ZX] = z.getX();
  m_array[ZY] = z.getY();
  m_array[ZZ] = z.getZ();
}

//*********************************************************************

Tensor::~Tensor() {}

//*********************************************************************

const Tensor Tensor::defaultTensor = Tensor();

//*********************************************************************

Tensor Tensor::defaultTensorNonConst = Tensor();

//*********************************************************************

void Tensor::setXX(const double& xx) { m_array[XX] = xx; }

//*********************************************************************

void Tensor::setXY(const double& xy) { m_array[XY] = xy; }

//*********************************************************************

void Tensor::setXZ(const double& xz) { m_array[XZ] = xz; }

//*********************************************************************

void Tensor::setYX(const double& yx) { m_array[YX] = yx; }

//*********************************************************************

void Tensor::setYY(const double& yy) { m_array[YY] = yy; }

//*********************************************************************

void Tensor::setYZ(const double& yz) { m_array[YZ] = yz; }

//*********************************************************************

void Tensor::setZX(const double& zx) { m_array[ZX] = zx; }

//*********************************************************************

void Tensor::setZY(const double& zy) { m_array[ZY] = zy; }

//*********************************************************************

void Tensor::setZZ(const double& zz) { m_array[ZZ] = zz; }

//*********************************************************************

void Tensor::setTensor(const Tensor& tensor)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] = tensor.m_array[i];
  }
}

//*********************************************************************

void Tensor::setTensor(const double& value)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] = value;
  }
}

//*********************************************************************

void Tensor::identity()
{
  m_array[XX] = 1.;
  m_array[XY] = 0.;
  m_array[XZ] = 0.;

  m_array[YX] = 0.;
  m_array[YY] = 1.;
  m_array[YZ] = 0.;

  m_array[ZX] = 0.;
  m_array[ZY] = 0.;
  m_array[ZZ] = 1.;
}

//*********************************************************************

void Tensor::correctZeros()
{
  for (int i = 0; i < 9; i++) {
    if (std::fabs(m_array[i]) < 1.e-20) m_array[i] = 0.;
  }
}

//*********************************************************************

Coord Tensor::scalar(const Coord& a) const
{
  coordBuff.setX(m_array[XX] * a.getX() + m_array[XY] * a.getY() + m_array[XZ] * a.getZ());
  coordBuff.setY(m_array[YX] * a.getX() + m_array[YY] * a.getY() + m_array[YZ] * a.getZ());
  coordBuff.setZ(m_array[ZX] * a.getX() + m_array[ZY] * a.getY() + m_array[ZZ] * a.getZ());
  return coordBuff;
}

//*********************************************************************

void Tensor::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  //T = T * L(n,t,b)
  tensorBuff.m_array[XX] = m_array[XX] * normal.getX()   + m_array[XY] * normal.getY()   + m_array[XZ] * normal.getZ();
  tensorBuff.m_array[XY] = m_array[XX] * tangent.getX()  + m_array[XY] * tangent.getY()  + m_array[XZ] * tangent.getZ();
  tensorBuff.m_array[XZ] = m_array[XX] * binormal.getX() + m_array[XY] * binormal.getY() + m_array[XZ] * binormal.getZ();

  tensorBuff.m_array[YX] = m_array[YX] * normal.getX()   + m_array[YY] * normal.getY()   + m_array[YZ] * normal.getZ();
  tensorBuff.m_array[YY] = m_array[YX] * tangent.getX()  + m_array[YY] * tangent.getY()  + m_array[YZ] * tangent.getZ();
  tensorBuff.m_array[YZ] = m_array[YX] * binormal.getX() + m_array[YY] * binormal.getY() + m_array[YZ] * binormal.getZ();

  tensorBuff.m_array[ZX] = m_array[ZX] * normal.getX()   + m_array[ZY] * normal.getY()   + m_array[ZZ] * normal.getZ();
  tensorBuff.m_array[ZY] = m_array[ZX] * tangent.getX()  + m_array[ZY] * tangent.getY()  + m_array[ZZ] * tangent.getZ();
  tensorBuff.m_array[ZZ] = m_array[ZX] * binormal.getX() + m_array[ZY] * binormal.getY() + m_array[ZZ] * binormal.getZ();

  *this = tensorBuff;
  
  //T = L(n,t,b)^T * T
  tensorBuff.m_array[XX] = m_array[XX] * normal.getX() + m_array[YX] * normal.getY() + m_array[ZX] * normal.getZ();
  tensorBuff.m_array[XY] = m_array[XY] * normal.getX() + m_array[YY] * normal.getY() + m_array[ZY] * normal.getZ();
  tensorBuff.m_array[XZ] = m_array[XZ] * normal.getX() + m_array[YZ] * normal.getY() + m_array[ZZ] * normal.getZ();

  tensorBuff.m_array[YX] = m_array[XX] * tangent.getX() + m_array[YX] * tangent.getY() + m_array[ZX] * tangent.getZ();
  tensorBuff.m_array[YY] = m_array[XY] * tangent.getX() + m_array[YY] * tangent.getY() + m_array[ZY] * tangent.getZ();
  tensorBuff.m_array[YZ] = m_array[XZ] * tangent.getX() + m_array[YZ] * tangent.getY() + m_array[ZZ] * tangent.getZ();

  tensorBuff.m_array[ZX] = m_array[XX] * binormal.getX() + m_array[YX] * binormal.getY() + m_array[ZX] * binormal.getZ();
  tensorBuff.m_array[ZY] = m_array[XY] * binormal.getX() + m_array[YY] * binormal.getY() + m_array[ZY] * binormal.getZ();
  tensorBuff.m_array[ZZ] = m_array[XZ] * binormal.getX() + m_array[YZ] * binormal.getY() + m_array[ZZ] * binormal.getZ();

  *this = tensorBuff;
}

//*********************************************************************
void Tensor::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  tensorBuff.m_array[XX] = m_array[XX] * normal.getX() + m_array[XY] * tangent.getX() + m_array[XZ] * binormal.getX();
  tensorBuff.m_array[XY] = m_array[XX] * normal.getY() + m_array[XY] * tangent.getY() + m_array[XZ] * binormal.getY();
  tensorBuff.m_array[XZ] = m_array[XX] * normal.getZ() + m_array[XY] * tangent.getZ() + m_array[XZ] * binormal.getZ();

  tensorBuff.m_array[YX] = m_array[YX] * normal.getX() + m_array[YY] * tangent.getX() + m_array[YZ] * binormal.getX();
  tensorBuff.m_array[YY] = m_array[YX] * normal.getY() + m_array[YY] * tangent.getY() + m_array[YZ] * binormal.getY();
  tensorBuff.m_array[YZ] = m_array[YX] * normal.getZ() + m_array[YY] * tangent.getZ() + m_array[YZ] * binormal.getZ();

  tensorBuff.m_array[ZX] = m_array[ZX] * normal.getX() + m_array[ZY] * tangent.getX() + m_array[ZZ] * binormal.getX();
  tensorBuff.m_array[ZY] = m_array[ZX] * normal.getY() + m_array[ZY] * tangent.getY() + m_array[ZZ] * binormal.getY();
  tensorBuff.m_array[ZZ] = m_array[ZX] * normal.getZ() + m_array[ZY] * tangent.getZ() + m_array[ZZ] * binormal.getZ();

  *this = tensorBuff;

  tensorBuff.m_array[XX] = m_array[XX] * normal.getX() + m_array[YX] * tangent.getX() + m_array[ZX] * binormal.getX();
  tensorBuff.m_array[XY] = m_array[XY] * normal.getX() + m_array[YY] * tangent.getX() + m_array[ZY] * binormal.getX();
  tensorBuff.m_array[XZ] = m_array[XZ] * normal.getX() + m_array[YZ] * tangent.getX() + m_array[ZZ] * binormal.getX();

  tensorBuff.m_array[YX] = m_array[XX] * normal.getY() + m_array[YX] * tangent.getY() + m_array[ZX] * binormal.getY();
  tensorBuff.m_array[YY] = m_array[XY] * normal.getY() + m_array[YY] * tangent.getY() + m_array[ZY] * binormal.getY();
  tensorBuff.m_array[YZ] = m_array[XZ] * normal.getY() + m_array[YZ] * tangent.getY() + m_array[ZZ] * binormal.getY();

  tensorBuff.m_array[ZX] = m_array[XX] * normal.getZ() + m_array[YX] * tangent.getZ() + m_array[ZX] * binormal.getZ();
  tensorBuff.m_array[ZY] = m_array[XY] * normal.getZ() + m_array[YY] * tangent.getZ() + m_array[ZY] * binormal.getZ();
  tensorBuff.m_array[ZZ] = m_array[XZ] * normal.getZ() + m_array[YZ] * tangent.getZ() + m_array[ZZ] * binormal.getZ();

  *this = tensorBuff;
}

//*********************************************************************

void Tensor::setTensorByLines(const Coord& x, const Coord& y, const Coord& z)
{
  m_array[XX] = x.getX();
  m_array[XY] = x.getY();
  m_array[XZ] = x.getZ();

  m_array[YX] = y.getX();
  m_array[YY] = y.getY();
  m_array[YZ] = y.getZ();

  m_array[ZX] = z.getX();
  m_array[ZY] = z.getY();
  m_array[ZZ] = z.getZ();
}

//*********************************************************************

void Tensor::setTensorByColumns(const Coord& x, const Coord& y, const Coord& z)
{
  m_array[XX] = x.getX();
  m_array[YX] = x.getY();
  m_array[ZX] = x.getZ();

  m_array[XY] = y.getX();
  m_array[YY] = y.getY();
  m_array[ZY] = y.getZ();

  m_array[XZ] = z.getX();
  m_array[YZ] = z.getY();
  m_array[ZZ] = z.getZ();
}

//*********************************************************************

void Tensor::tensorToCoords(Coord& x, Coord& y, Coord& z) const
{
  x.setXYZ(m_array[XX], m_array[XY], m_array[XZ]);
  y.setXYZ(m_array[YX], m_array[YY], m_array[YZ]);
  z.setXYZ(m_array[ZX], m_array[ZY], m_array[ZZ]);
}

//*********************************************************************

void Tensor::tensorToArray(double* array) const
{
  for (int i = 0; i < 9; i++) {
    array[i] = m_array[i];
  }
}

//*********************************************************************

void Tensor::arrayToTensor(const double* array)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] = array[i];
  }
}

//*********************************************************************

void Tensor::transpose(Tensor& transposedTensor) const
{
  transposedTensor.m_array[XX] = m_array[XX];
  transposedTensor.m_array[XY] = m_array[YX];
  transposedTensor.m_array[XZ] = m_array[ZX];

  transposedTensor.m_array[YX] = m_array[XY];
  transposedTensor.m_array[YY] = m_array[YY];
  transposedTensor.m_array[YZ] = m_array[ZY];

  transposedTensor.m_array[ZX] = m_array[XZ];
  transposedTensor.m_array[ZY] = m_array[YZ];
  transposedTensor.m_array[ZZ] = m_array[ZZ];
}

//*********************************************************************

void Tensor::matrixProduct(const Tensor& tensor2, Tensor& resultingTensor) const
{
  //Details if made by a loop (high computational cost considering the structure of Tensor instead of a two-dimensional array)
  // for (int i = 0; i < 3; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     resultingTensor(i, j) = 0.;
  //     for (int k = 0; k < 3; k++) {
  //       resultingTensor(i, j) += *this->getElement(i, k) * tensor2.getElement(k, j);
  //     }
  //   } 
  // }
  resultingTensor.m_array[XX] = m_array[XX]*tensor2.m_array[XX] + m_array[XY]*tensor2.m_array[YX] + m_array[XZ]*tensor2.m_array[ZX];
  resultingTensor.m_array[XY] = m_array[XX]*tensor2.m_array[XY] + m_array[XY]*tensor2.m_array[YY] + m_array[XZ]*tensor2.m_array[ZY];
  resultingTensor.m_array[XZ] = m_array[XX]*tensor2.m_array[XZ] + m_array[XY]*tensor2.m_array[YZ] + m_array[XZ]*tensor2.m_array[ZZ];

  resultingTensor.m_array[YX] = m_array[YX]*tensor2.m_array[XX] + m_array[YY]*tensor2.m_array[YX] + m_array[YZ]*tensor2.m_array[ZX];
  resultingTensor.m_array[YY] = m_array[YX]*tensor2.m_array[XY] + m_array[YY]*tensor2.m_array[YY] + m_array[YZ]*tensor2.m_array[ZY];
  resultingTensor.m_array[YZ] = m_array[YX]*tensor2.m_array[XZ] + m_array[YY]*tensor2.m_array[YZ] + m_array[YZ]*tensor2.m_array[ZZ];

  resultingTensor.m_array[ZX] = m_array[ZX]*tensor2.m_array[XX] + m_array[ZY]*tensor2.m_array[YX] + m_array[ZZ]*tensor2.m_array[ZX];
  resultingTensor.m_array[ZY] = m_array[ZX]*tensor2.m_array[XY] + m_array[ZY]*tensor2.m_array[YY] + m_array[ZZ]*tensor2.m_array[ZY];
  resultingTensor.m_array[ZZ] = m_array[ZX]*tensor2.m_array[XZ] + m_array[ZY]*tensor2.m_array[YZ] + m_array[ZZ]*tensor2.m_array[ZZ];
}

//*********************************************************************

double Tensor::trace() const
{
  return m_array[XX] + m_array[YY] + m_array[ZZ];
}

//*********************************************************************

double Tensor::determinant() const
{
  return   m_array[XX] * (m_array[YY]*m_array[ZZ] - m_array[YZ]*m_array[ZY])
         - m_array[YX] * (m_array[XY]*m_array[ZZ] - m_array[XZ]*m_array[ZY])
         + m_array[ZX] * (m_array[XY]*m_array[YZ] - m_array[XZ]*m_array[YY]);
}

//*********************************************************************

void Tensor::inverse(Tensor& inverseTensor) const
{
  inverseTensor.m_array[XX] = m_array[YY] * m_array[ZZ] - m_array[YZ] * m_array[ZY];
  inverseTensor.m_array[XY] = m_array[XZ] * m_array[ZY] - m_array[XY] * m_array[ZZ];
  inverseTensor.m_array[XZ] = m_array[XY] * m_array[YZ] - m_array[XZ] * m_array[YY];

  inverseTensor.m_array[YX] = m_array[YZ] * m_array[ZX] - m_array[YX] * m_array[ZZ];
  inverseTensor.m_array[YY] = m_array[XX] * m_array[ZZ] - m_array[XZ] * m_array[ZX];
  inverseTensor.m_array[YZ] = m_array[XZ] * m_array[YX] - m_array[XX] * m_array[YZ];

  inverseTensor.m_array[ZX] = m_array[YX] * m_array[ZY] - m_array[YY] * m_array[ZX];
  inverseTensor.m_array[ZY] = m_array[XY] * m_array[ZX] - m_array[XX] * m_array[ZY];
  inverseTensor.m_array[ZZ] = m_array[XX] * m_array[YY] - m_array[XY] * m_array[YX];

  double det = this->determinant();
  if (std::fabs(det) > 1.e-10) inverseTensor /= det;
  else Errors::errorMessage("Matrix/Tensor is not inversible (det = 0)");
}

//*********************************************************************

bool Tensor::isIdentity() const
{
  if (std::fabs(m_array[XX]) - 1. < 1.e-15 &&
      std::fabs(m_array[XY])      < 1.e-15 &&
      std::fabs(m_array[XZ])      < 1.e-15 &&

      std::fabs(m_array[YX])      < 1.e-15 &&
      std::fabs(m_array[YY]) - 1. < 1.e-15 &&
      std::fabs(m_array[YZ])      < 1.e-15 &&

      std::fabs(m_array[ZX])      < 1.e-15 &&
      std::fabs(m_array[ZY])      < 1.e-15 &&
      std::fabs(m_array[ZZ]) - 1. < 1.e-15)
    return true;
  else
    return false;
}

//*********************************************************************

void Tensor::eigen(Tensor& eigenvalues, Tensor& eigenvectors) const
{
  double abserr(1.e-9); //Abs tolerance [sum of (off-diagonal elements)^2]

  //Initialize eigenvalues and eigenvectors
  eigenvalues = *this;
  eigenvectors.identity();
    
  //Sum of all off-diagonal elements (squared)
  double b2(0.);
  b2 = m_array[XY]*m_array[XY] + m_array[XZ]*m_array[XZ] + m_array[YZ]*m_array[YZ] +
       m_array[YX]*m_array[YX] + m_array[ZX]*m_array[ZX] + m_array[ZY]*m_array[ZY];
  if (b2 <= abserr) return;

  //Average for off-diagonal elements divided by 2: 
  double bar = 0.5 * b2 / 9.;

  //Compute eigenvalues and eigenvectors
  int ii(0), jj(0), ji(0), ik(0), jk(0), ki(0), kj(0);
  double beta(0.), coeff(0.), c(0.), s(0.), cs(0.), sc(0.);
  while (b2 > abserr) {
    for (int i = 0; i < 2; i++) {
      for (int j = i + 1; j < 3; j++) {
        //Compute element indexes for loop i and j
        ii = i*3 + i;
        jj = j*3 + j;
        ji = j*3 + i;

        //Compute new b2 and bar
        if (eigenvalues.m_array[ji] * eigenvalues.m_array[ji] <= bar) continue; //Do not touch small elements
        b2 -= 2. * eigenvalues.m_array[ji] * eigenvalues.m_array[ji];
        bar = 0.5 * b2 / 9.;

        //Compute coefficients c and s for Givens matrix
        beta = (eigenvalues.m_array[jj] - eigenvalues.m_array[ii]) / (2. * eigenvalues.m_array[ji]);
        coeff = 0.5 * beta / sqrt(1. + beta*beta);
        s = sqrt(std::max(0.5 + coeff, 0.));
        c = sqrt(std::max(0.5 - coeff, 0.));

        //Recompute rows i and j of eigenvalues
        for (int k = 0; k < 3; k++) {
          ik = i*3 + k;
          jk = j*3 + k;
          cs =  c * eigenvalues.m_array[ik] + s * eigenvalues.m_array[jk];
          sc = -s * eigenvalues.m_array[ik] + c * eigenvalues.m_array[jk];
          eigenvalues.m_array[ik] = cs;
          eigenvalues.m_array[jk] = sc;
        }

        //New eigenvalues and eigenvectors
        for (int k = 0; k < 3; k++) {
          ki = k*3 + i;
          kj = k*3 + j;
          cs =  c * eigenvalues.m_array[ki] + s * eigenvalues.m_array[kj];
          sc = -s * eigenvalues.m_array[ki] + c * eigenvalues.m_array[kj];
          eigenvalues.m_array[ki] = cs;
          eigenvalues.m_array[kj] = sc;
          cs =  c * eigenvectors.m_array[ki] + s * eigenvectors.m_array[kj];
          sc = -s * eigenvectors.m_array[ki] + c * eigenvectors.m_array[kj];
          eigenvectors.m_array[ki] = cs;
          eigenvectors.m_array[kj] = sc;
        }
      }
    }
  }
}

//*******************SURCHARGES OPERATEURS*****************************

Tensor& Tensor::operator=(const double& scalar)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] = scalar;
  }
  return *this;
}

//*********************************************************************

Tensor& Tensor::operator=(const Tensor& a)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] = a.m_array[i];
  }
  return *this;
}

//*********************************************************************

Tensor& Tensor::operator*= (const double& scalar)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] *= scalar;
  }
  return *this;
}

//*********************************************************************

Tensor& Tensor::operator/= (const double& scalar)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] /= scalar;
  }
  return *this;
}

//*********************************************************************

Tensor Tensor::operator* (const double& scalar) const
{
  Tensor copie(*this);
  copie *= scalar;
  return copie;
}

//*********************************************************************

Tensor Tensor::operator/ (const double& scalar) const
{
  Tensor copie(*this);
  copie /= scalar;
  return copie;
}

//*********************************************************************

Tensor& Tensor::operator+= (const Tensor& a)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] += a.m_array[i];
  }
  return *this;
}

//*********************************************************************

Tensor& Tensor::operator-= (const Tensor& a)
{
  for (int i = 0; i < 9; i++) {
    m_array[i] -= a.m_array[i];
  }
  return *this;
}

//*********************************************************************
//Surcharge operateur externe a la classe car prends deux arguments

Tensor operator* (const double& scalar, const Tensor& a)
{
  Tensor copy(a);
  copy *= scalar;
  return copy;
}

//*********************************************************************

Tensor operator+ (const Tensor& a, const Tensor& b)
{
  Tensor copy(a);
  copy += b;
  return copy;
}

//*********************************************************************

Tensor operator- (const Tensor& a, const Tensor& b)
{
  Tensor copy(a);
  copy -= b;
  return copy;
}

//*********************************************************************