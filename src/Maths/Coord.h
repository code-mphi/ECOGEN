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

#ifndef COORD_H
#define COORD_H

//! \file      Coord.h
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

//! \class     Coord
//! \brief     Class for a coordinate system object such as coordinates of the vertex or a vector
class Coord
{
public:
  Coord();
  //! \brief     Coord constructor
  //! \param     x                    value of the x-direction coordinate
  //! \param     y                    value of the y-direction coordinate (if it is not assigned it is set to 0.)
  //! \param     z                    value of the z-direction coordinate (if it is not assigned it is set to 0.)
  Coord(const double &x, const double &y = 0., const double &z = 0.);
  ~Coord();
  //! \brief     Method for defaultCoord (const version)
  //! \details   Used for the defaultCoord object
  const Coord& coord() const;
  //! \brief     Method for defaultCoord (non-const version)
  //! \details   Used for the defaultCoord object
  Coord& coord();
  //! \brief     Default Coord object (const version)
  //! \details   Used when returning a const Coord&
  static const Coord defaultCoord;
  //! \brief     Default Coord object (non-const version)
  //! \details   Used when returning a const Coord&
  static Coord defaultCoordNonConst;
  //! \brief     Set the values of the Coord object
  //! \param     x                    value of the x-direction coordinate
  //! \param     y                    value of the y-direction coordinate
  //! \param     z                    value of the z-direction coordinate
  void setXYZ(const double &x, const double & y, const double &z);
  //! \brief     Set the value in the x-direction of the Coord object
  //! \param     x                    value of the x-direction coordinate
  void setX(const double &x);
  //! \brief     Set the value in the y-direction of the Coord object
  //! \param     y                    value of the y-direction coordinate
  void setY(const double &y);
  //! \brief     Set the value in the z-direction of the Coord object
  //! \param     z                    value of the z-direction coordinate
  void setZ(const double &z);
  //! \brief     Return the value in the x-direction of the Coord object
  const double& getX() const { return m_x; }
  //! \brief     Return the value in the y-direction of the Coord object
  const double& getY() const { return m_y; }
  //! \brief     Return the value in the z-direction of the Coord object
  const double& getZ() const { return m_z; }
  //! \brief     Return the value of the norm of the Coord object
  double norm() const;
  //! \brief     Return the value of the squared norm of the Coord object
  double squaredNorm() const;
  //! \brief     Return a Coord object with absolute values of each component
  Coord abs() const;

  //! \brief     Scalar product between the present vector and vector a
  //! \param     a                    vector (Coord)
  double scalar(const Coord &a) const;
  //! \brief     Cross product between the present vector and vector a
  //! \param     a                    vector (Coord)
  Coord cross(const Coord &a) const;
  //! \brief     Projection in the local coordinate system which is defined by the transmitted normal, tangent and binormal
  //! \param     normal               normal vector (Coord)
  //! \param     tangent              tangent vector (Coord)
  //! \param     binormal             binormal vector (Coord)
  void localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal);
  //! \brief     Reverse projection in the absolute cartesian coordinate system
  //! \param     normal               normal vector (Coord)
  //! \param     tangent              tangent vector (Coord)
  //! \param     binormal             binormal vector (Coord)
  void reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal);

  //! \brief     Set a vector from the result of the substraction of the vector a from the vector b
  //! \param     a                    vector (Coord)
  //! \param     b                    vector (Coord)
  void setFromSubtractedVectors(const Coord &a, const Coord &b);
  //! \brief     Return the scalar product bewteen two vectors
  //! \param     v1                   vector (Coord)
  //! \param     v2                   vector (Coord)
  static double scalarProduct(const Coord &v1, const Coord &v2);
  //! \brief     Return the cross product bewteen two vectors
  //! \param     v1                   vector (Coord)
  //! \param     v2                   vector (Coord)
  static Coord crossProduct(const Coord &v1, const Coord &v2);
  //! \brief     Change the sign if the present vector
  void changeSign();
  //! \brief     Divide the present vector by its norm
  void normalized();
  //! \brief     Print the information of the vector
  void printInfo() const;

  //! \brief     Compute the determinant of the matrix formed by the vectors v1, v2 and v3
  //! \param     v1                   vector (Coord)
  //! \param     v2                   vector (Coord)
  //! \param     v3                   vector (Coord)
  static double determinant(const Coord &v1, const Coord &v2, const Coord &v3);
  
  //! \brief     Compute some sort of cosinus between two vectors
  //! \param     v1                   vector (Coord)
  //! \param     v2                   vector (Coord)
  static double cos(const Coord &v1, const Coord &v2);
  //! \brief     Compute some sort of sinus between two vectors
  //! \param     v1                   vector (Coord)
  //! \param     v2                   vector (Coord)
  static Coord sin(const Coord &v1, const Coord &v2);

  //Operator surcharges
  Coord& operator=(const double &scalar);
  Coord& operator*= (const double &scalar);
  Coord& operator/= (const double &scalar);
  Coord operator* (const double &scalar) const; 
  Coord operator/ (const double &scalar) const;
  Coord& operator+= (const Coord &a);
  Coord& operator-= (const Coord &a);

protected:
  double m_x;          //! Value in the x-direction
  double m_y;          //! Value in the y-direction
  double m_z;          //! Value in the z-direction
};

//Extern operator surcharges of the class because they take two arguments
Coord operator* (const double &scalar, const Coord &a);
Coord operator+ (const Coord &a, const Coord &b);
Coord operator- (const Coord &a, const Coord &b);

#endif // COORD_H
