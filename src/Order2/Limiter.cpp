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

#include "Limiter.h"

//***********************************************************************

Limiter::Limiter() : m_limType(LimiterType::NONE) {}

//***********************************************************************

Limiter::~Limiter(){}

//***********************************************************************

double Limiter::computeGradientLimiter(double val, double min, double max, double slope) const
{
  double eps(1.e-6);
  double phi(1.);

  if (slope > eps * val) {
    phi = (max - val) * 0.5 / slope;
  }
  else if (slope < - eps * val) {
    phi = (min - val) * 0.5 / slope;
  }

  double buff1(0.), buff2(0.);
  buff1 = std::min(double(m_limType) * phi, 1.);
  buff2 = std::min(phi, double(m_limType));
  double theta(0.);
  theta = std::max(theta, buff1);
  theta = std::max(theta, buff2);
  return theta;
}

//***********************************************************************