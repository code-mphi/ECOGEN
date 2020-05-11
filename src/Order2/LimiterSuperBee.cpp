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

//! \file      LimiterSuperBee.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "LimiterSuperBee.h"

//***********************************************************************

LimiterSuperBee::LimiterSuperBee(){}

//***********************************************************************

LimiterSuperBee::~LimiterSuperBee(){}

//***********************************************************************

double LimiterSuperBee::limiteSlope(const double& slope1, const double& slope2) 
{
	double zero(1e-6);
	double slope(0.), produit(slope1*slope2);
	if (produit>zero){
		slope = std::max(std::min(2.*std::fabs(slope1), std::fabs(slope2)), std::min(std::fabs(slope1),2.*std::fabs(slope2)));
	}
	if (slope1 < 0.) { slope = -slope; }
	return slope;
}
