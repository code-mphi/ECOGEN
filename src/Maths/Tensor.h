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

#ifndef TENSOR_H
#define TENSOR_H

#include "../Errors.h"
#include <array>

class Tensor;
#include "Coord.h"

//! \brief     Enumeration for the tensor elements
enum TensorElement { XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };

//! \class     Tensor
//! \brief     Class for a matrix 3x3 system object
class Tensor
{
  public:
    Tensor();
    Tensor(const Tensor& tensor);
    Tensor(const double &xx, const double &xy, const double &xz, const double &yx, const double &yy, const double &yz, const double &zx, const double &zy, const double &zz);
    Tensor(const Coord &x, const Coord &y, const Coord &z);
    ~Tensor();

    //! \brief     Default Tensor object (const version)
    //! \details   Used when returning a const Tensor&
    static const Tensor defaultTensor;
    //! \brief     Default Tensor object (non-const version)
    //! \details   Used when returning a const Tensor&
    static Tensor defaultTensorNonConst;

    //! \brief     Return the value of the specific element of the Tensor object
    const double& getElement(const int& element) const { return m_array[element]; }
    //! \brief     Return the value of the xx-element of the Tensor object
    const double& getXX() const { return m_array[XX]; }
    //! \brief     Return the value of the xy-element of the Tensor object
    const double& getXY() const { return m_array[XY]; }
    //! \brief     Return the value of the xz-element of the Tensor object
    const double& getXZ() const { return m_array[XZ]; }
    //! \brief     Return the value of the yx-element of the Tensor object
    const double& getYX() const { return m_array[YX]; }
    //! \brief     Return the value of the yy-element of the Tensor object
    const double& getYY() const { return m_array[YY]; }
    //! \brief     Return the value of the yz-element of the Tensor object
    const double& getYZ() const { return m_array[YZ]; }
    //! \brief     Return the value of the zx-element of the Tensor object
    const double& getZX() const { return m_array[ZX]; }
    //! \brief     Return the value of the zy-element of the Tensor object
    const double& getZY() const { return m_array[ZY]; }
    //! \brief     Return the value of the zz-element of the Tensor object
    const double& getZZ() const { return m_array[ZZ]; }
    //! \brief     Set the value in the xx-element of the Tensor object
    //! \param     xx                   value of the xx-element coordinate
    void setXX(const double& xx);
    //! \brief     Set the value in the xy-element of the Tensor object
    //! \param     xy                   value of the xy-element coordinate
    void setXY(const double& xy);
    //! \brief     Set the value in the xz-element of the Tensor object
    //! \param     xz                   value of the xz-element coordinate
    void setXZ(const double& xz);
    //! \brief     Set the value in the yx-element of the Tensor object
    //! \param     xx                   value of the yx-element coordinate
    void setYX(const double& xx);
    //! \brief     Set the value in the yy-element of the Tensor object
    //! \param     xy                   value of the yy-element coordinate
    void setYY(const double& xy);
    //! \brief     Set the value in the yz-element of the Tensor object
    //! \param     xz                   value of the yz-element coordinate
    void setYZ(const double& xz);
    //! \brief     Set the value in the zx-element of the Tensor object
    //! \param     xx                   value of the zx-element coordinate
    void setZX(const double& xx);
    //! \brief     Set the value in the zy-element of the Tensor object
    //! \param     xy                   value of the zy-element coordinate
    void setZY(const double& xy);
    //! \brief     Set the value in the zz-element of the Tensor object
    //! \param     xz                   value of the zz-element coordinate
    void setZZ(const double& xz);
    //! \brief     Set the Tensor object
    //! \param     tensor               tensor
    void setTensor(const Tensor& tensor);
    //! \brief     Set all the elements of the Tensor object to a scalar value
    //! \param     value                scalar value to set all elements
    void setTensor(const double& value);
    //! \brief     Set the Tensor object to identity
    void identity();
    //! \brief     Set to perfect zero each of the epsilon terms
    void correctZeros();

    //! \brief     Return the scalar product between the Tensor and the Coord
    //! \param     a                    vector (Coord)
    Coord scalar(const Coord& a) const;
    //! \brief     Projection in the local coordinate system which is defined by the transmitted normal, tangent and binormal
    //! \param     normal               normal vector (Coord)
    //! \param     tangent              tangent vector (Coord)
    //! \param     binormal             binormal vector (Coord)
    void localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal);
    //! \brief     Reverse projection in the absolute Cartesian coordinate system
    //! \param     normal               normal vector (Coord)
    //! \param     tangent              tangent vector (Coord)
    //! \param     binormal             binormal vector (Coord)
    void reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal);
    void setTensorByLines(const Coord& x, const Coord& y, const Coord& z);
    void setTensorByColumns(const Coord& x, const Coord& y, const Coord& z);

    //! \brief     Transform the Tensor into Coords (vectors)
    //! \param     x                    x vector (Coord)
    //! \param     y                    y vector (Coord)
    //! \param     z                    z vector (Coord)
    void tensorToCoords(Coord& x, Coord& y, Coord& z) const;
    //! \brief     Transform the Tensor into 1D array (pointer)
    //! \param     array                1D array (pointer)
    void tensorToArray(double* array) const;
    //! \brief     Transform the 1D array (pointer) into Tensor
    //! \param     array                1D array (pointer)
    void arrayToTensor(const double* array);

    //! \brief     Compute the transpose of the tensor
    //! \param     transposedTensor     transposed tensor (Tensor)
    void transpose(Tensor& transposedTensor) const;
    //! \brief     Compute the matrix product of the two tensors
    //! \param     tensor2              second tensor (Tensor)
    //! \param     resultingTensor      resulting tensor of the product (Tensor)
    void matrixProduct(const Tensor& tensor2, Tensor& resultingTensor) const;
    //! \brief     Return the trace of the tensor
    double trace() const;
    //! \brief     Return the determinant of the tensor
    double determinant() const;
    //! \brief     Compute the inverse of the tensor
    //! \param     inverseTensor        inverse tensor (Tensor)
    void inverse(Tensor& inverseTensor) const;

    //! \brief     Return true if equal to identity
    bool isIdentity() const;

    //! \brief     Compute the eigenvalues and eigenvectors of the tensor with Jacobi's algorithm
    //! \param     eigenvalues          tensor whose diagonal elements are eigenvalues (Tensor)
    //! \param     eigenvectors         tensor whose columns are eigenvectors (Tensor)
    void eigen(Tensor& eigenvalues, Tensor& eigenvectors) const;

    //Operator surcharges
    Tensor& operator=(const double& scalar);
    Tensor& operator=(const Tensor& a);
    Tensor& operator*= (const double& scalar);
    Tensor& operator/= (const double& scalar);
    Tensor operator* (const double& scalar) const; 
    Tensor operator/ (const double& scalar) const;
    Tensor& operator+= (const Tensor& a);
    Tensor& operator-= (const Tensor& a);

  protected:
    std::array<double, 9> m_array;       //! Elements of the 3x3 matrix
};

extern Tensor tensorBuff;
extern Tensor tensorIdentity;
extern Tensor tensorCobase;
extern Tensor tensorNonConsCobase;
extern Tensor tensorF;
extern Tensor tensorG;
extern Tensor tensorG2;
extern Tensor tensorA;
extern Tensor tensorEigenvalues;
extern Tensor tensorP;
extern Tensor tensorPinverse;
extern Tensor tensorD;

//Extern operator surcharges of the class because they take two arguments
Tensor operator* (const double& scalar, const Tensor& a);
Tensor operator+ (const Tensor& a, const Tensor& b);
Tensor operator- (const Tensor& a, const Tensor& b);

#endif // TENSOR_H
