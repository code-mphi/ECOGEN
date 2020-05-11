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

//! \file      Transport.cpp
//! \author    K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "Transport.h"

using namespace tinyxml2;

Transport* fluxBufferTransport;

//***********************************************************************

Transport::Transport() : m_value(0.)
{}

//***********************************************************************

Transport::~Transport(){}

//***********************************************************************

void Transport::setValue(double value)
{
  m_value = value;
}

//***********************************************************************

void Transport::solveRiemann(double transportLeft, double transportRight, double sM)
{
	if (sM > 0.) { m_value = transportLeft*sM; }
	else { m_value = transportRight*sM; }
}

//***********************************************************************

void Transport::solveRiemannWall()
{
	m_value = 0.;
}

//***********************************************************************

void Transport::solveRiemannInflow(double transportLeft, double sM, double valueTransport)
{
  if (sM > 0.) { m_value = transportLeft*sM; }
  else { m_value = valueTransport*sM; }
}

//***********************************************************************

void Transport::solveRiemannTank(double transportLeft, double sM, double valueTransport)
{
  if (sM > 0.) { m_value = transportLeft*sM; }
  else { m_value = valueTransport*sM; }
}

//***********************************************************************

void Transport::solveRiemannOutflow(double transportLeft, double sM, double valueTransport)
{
	if (sM > 0.) { m_value = transportLeft*sM; }
  else { m_value = valueTransport*sM; }
}

//***********************************************************************

void Transport::addFlux(double coefA, const int num)
{
  m_value += coefA*fluxBufferTransport[num].m_value;
}

//***********************************************************************

void Transport::subtractFlux(double coefA, const int num)
{
  m_value -= coefA* fluxBufferTransport[num].m_value;
}

//***********************************************************************

void Transport::addNonCons(double coefA, double transport, const double sM)
{
  m_value += -coefA*transport*sM;
}

//***********************************************************************

void Transport::subtractNonCons(double coefA, double transport, const double sM)
{
  m_value -= -coefA*transport*sM;
}

//***********************************************************************

void Transport::multiply(double scalar) 
{
  m_value *= scalar;
}

//***********************************************************************

void Transport::add(double scalar)
{
  m_value += scalar;
}

//***************************************************************************

void Transport::changeSign()
{
  m_value = -m_value;
}

//***************************************************************************

void Transport::computeSlopeTransport(const double valueLeft, const double valueRight, const double &distance)
{
  m_value = (valueRight - valueLeft) / distance;
}

//***********************************************************************

void Transport::extrapolate(const double &slope, const double &distance)
{
  m_value += slope * distance;
}