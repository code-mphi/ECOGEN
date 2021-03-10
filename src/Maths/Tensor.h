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

#include "Coord.h"

class Tensor
{
  public:
    Tensor();
    Tensor(const double &xx, const double &xy, const double &xz, const double &yx, const double &yy, const double &yz, const double &zx, const double &zy, const double &zz);
    Tensor(const Coord &x, const Coord &y, const Coord &z);
    ~Tensor();

    void scalar(const Coord& a) const;
    void localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal);
    void setTensorByLines(const Coord& x, const Coord& y, const Coord& z);
    void setTensorByColumns(const Coord& x, const Coord& y, const Coord& z);

    void tensorToCoords(Coord& x, Coord& y, Coord& z);

  private:
    double m_xx;
    double m_xy;
    double m_xz;
    double m_yx; 
    double m_yy;
    double m_yz;
    double m_zx;
    double m_zy;
    double m_zz;
};

extern Tensor projectedTensor;

#endif // TENSOR_H
