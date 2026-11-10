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

#include "QAPTurbulence.h"

//***********************************************************************

QAPTurbulence::QAPTurbulence(AddPhys* addPhys) : QuantitiesAddPhys(addPhys), m_grads(4), m_turbulentViscosity(0.)
{
  variableNamesTurb.resize(4);
  numPhasesTurb.resize(4);
  for (int i = 0; i < 4; ++i) {
    m_grads[i]       = 0.;
    numPhasesTurb[i] = 0;
  }
  variableNamesTurb[0] = Variable::velocityU;
  variableNamesTurb[1] = Variable::velocityV;
  variableNamesTurb[2] = Variable::velocityW;
  variableNamesTurb[3] = Variable::temperature;
}

//***********************************************************************

QAPTurbulence::~QAPTurbulence() {}

//***********************************************************************

void QAPTurbulence::computeQuantities(Cell* cell) { cell->computeGradients(m_grads, variableNamesTurb, numPhasesTurb); }

//***********************************************************************

void QAPTurbulence::setGrad(const Coord& grad, const int& num) { m_grads[num - 1] = grad; }

//***********************************************************************

void QAPTurbulence::setTurbulentViscosity(const double& turbulentViscosity) { m_turbulentViscosity = turbulentViscosity; }

//***********************************************************************

double QAPTurbulence::getTurbulentViscosity() const { return m_turbulentViscosity; }

//***********************************************************************
