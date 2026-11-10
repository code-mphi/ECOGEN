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

#include "Tools.h"

Tools* TB;

SimulationOrder simulationOrder = SimulationOrder::first;
bool turbulenceActivated        = false;
int turbulenceNumQPA            = 0;
double sum_u2v2w2               = 0.;
long nCells                     = 0.;
double initialSoundSpeed        = 0.;
double initialConductivity      = 0.;
double initialViscosity         = 0.;

//***********************************************************************

Tools::Tools(const int& numbPhases, const int& numbSolids, const int& numbTransports)
{
  ak.resize(numbPhases);
  ak4.resize(numbPhases);
  ak5.resize(numbPhases);
  Yk.resize(numbPhases);
  rhok.resize(numbPhases);
  pk.resize(numbPhases);
  ek.resize(numbPhases);
  ek4.resize(numbPhases);
  ek5.resize(numbPhases);
  Ek.resize(numbPhases);
  akS.resize(numbPhases);
  rhokS.resize(numbPhases);
  rhokStar.resize(numbPhases);
  pkStar.resize(numbPhases);
  ekStar.resize(numbPhases);
  EkStar.resize(numbPhases);
  YkStar.resize(numbPhases);
  vkStar.resize(numbPhases);
  Deltapk.resize(numbPhases);
  zk.resize(numbPhases);
  rho_cIksquare.resize(numbPhases);
  eos.resize(numbPhases);
  Hk0.resize(numbPhases);
  Yk0.resize(numbPhases);
  compactionPk_dkappadxi.resize(numbPhases);
  compactionPk_dkappadalpha.resize(numbPhases);
  dlambda.resize(numbPhases);
  dplast.resize(numbPhases);
  LogCobase.resize(numbPhases);
  LogCobase4.resize(numbPhases);
  LogCobase5.resize(numbPhases);
  LogLambda.resize(numbPhases);
  LogLambda4.resize(numbPhases);
  LogLambda5.resize(numbPhases);
  alphaNull.resize(numbPhases);
  relaxSolidPlast.resize(numbPhases);
  dLogCobaseDt1.resize(numbSolids);
  dLogCobaseDt2.resize(numbSolids);
  dLogCobaseDt3.resize(numbSolids);
  dLogCobaseDt4.resize(numbSolids);
  dLogCobaseDt5.resize(numbSolids);
  dLogCobaseDt6.resize(numbSolids);
  dLogLambdaDt1.resize(numbPhases);
  dLogLambdaDt2.resize(numbPhases);
  dLogLambdaDt3.resize(numbPhases);
  dLogLambdaDt4.resize(numbPhases);
  dLogLambdaDt5.resize(numbPhases);
  dLogLambdaDt6.resize(numbPhases);
  dakDt1.resize(numbPhases);
  dakDt2.resize(numbPhases);
  dakDt3.resize(numbPhases);
  dakDt4.resize(numbPhases);
  dakDt5.resize(numbPhases);
  dakDt6.resize(numbPhases);
  dekDt1.resize(numbPhases);
  dekDt2.resize(numbPhases);
  dekDt3.resize(numbPhases);
  dekDt4.resize(numbPhases);
  dekDt5.resize(numbPhases);
  dekDt6.resize(numbPhases);

  for (int k = 0; k < numbPhases; k++) {
    compactionPk_dkappadxi[k]    = 0.;
    compactionPk_dkappadalpha[k] = 0.;
    dlambda[k]                   = 0.;
    alphaNull[k]                 = false;
    relaxSolidPlast[k]           = false;
  }

  physicalTime = 0.;

  numberPhases     = numbPhases;
  numberSolids     = numbSolids;
  numberTransports = numbTransports;
}

//***********************************************************************

Tools::~Tools() {}

//***********************************************************************

void Tools::uppercase(std::string& string)
{
  for (unsigned int c = 0; c < string.size(); c++) {
    string[c] = toupper(string[c]);
  }
}

//***********************************************************************

void Tools::lowercase(std::string& string)
{
  for (unsigned int c = 0; c < string.size(); c++) {
    string[c] = tolower(string[c]);
  }
}

//***********************************************************************

void Tools::swap(double& a, double& b)
{
  double buff(a);
  a = b;
  b = buff;
}

//***********************************************************************

double Tools::returnNonZeroValue(double a)
{
  if (a > epsilonAlphaNull || a < -epsilonAlphaNull) {
    return a;
  }
  else {
    return sign(a) * epsilonAlphaNull;
  }
}

//***********************************************************************

double Tools::uselessDouble;
