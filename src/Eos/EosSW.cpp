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

#include "EosSW.h"

//***********************************************************************

EosSW::EosSW(std::vector<std::string>& nameParameterEos, int& number) : Eos(number) { nameParameterEos.push_back("gravity"); }

//***********************************************************************

EosSW::~EosSW() {}

//***********************************************************************

void EosSW::assignParametersEos(std::string name, std::vector<double> parametersEos)
{
  assert(parametersEos.size() == 1);

  m_name    = name;
  m_gravity = parametersEos[0];
}

//***********************************************************************

//Constant methods
//****************

double EosSW::computePressure(const double& height) const
{
  double p_;

  // TO FILL
  // Compute pressure p
  assert(0 && "implement pressure computation");
  //p_ =
  return p_;
}

//***********************************************************************

double EosSW::computeSoundSpeed(const double& height) const
{
  double c_;

  // TO FILL
  // Compute sound speed c
  assert(0 && "implement sound speed computation");
  //c_ =
  return c_;
}

//***********************************************************************
