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

Tools *TB;

//***********************************************************************

Tools::Tools(const int& numbPhases, const int& numbSolids, const int& numbTransports)
{
  ak = new double[numbPhases];
  Yk = new double[numbPhases];
  rhok = new double[numbPhases];
  pk = new double[numbPhases];
  ek = new double[numbPhases];
  Ek = new double[numbPhases];
  akS = new double[numbPhases];
  rhokS = new double[numbPhases];
  rhokStar = new double[numbPhases];
  pkStar = new double[numbPhases];
  ekStar = new double[numbPhases];
  EkStar = new double[numbPhases];
  YkStar = new double[numbPhases];
  vkStar = new double[numbPhases];
  Deltapk = new double[numbPhases];
  zk = new double[numbPhases];
  rho_cIksquare = new double[numbPhases];
  eos = new Eos*[numbPhases];
  Hk0 = new double[numbPhases];
  Yk0 = new double[numbPhases];
  compactionPk = new double[numbPhases];
  dlambda = new double[numbPhases];
  dplast = new Tensor[numbPhases];
  alphaNull = new bool[numbPhases];
  relaxSolidPlast = new bool[numbPhases];

  for (int k = 0; k < numbPhases; k++) {
    compactionPk[k] = 0.;
    alphaNull[k] = false;
    relaxSolidPlast[k] = false;
  }

  physicalTime = 0.;

  numberPhases = numbPhases;
  numberSolids = numbSolids;
  numberTransports = numbTransports;
}

//***********************************************************************

Tools::~Tools()
{
  if (ak != 0) delete[] ak;
  if (Yk != 0) delete[] Yk;
  if (rhok != 0) delete[] rhok;
  if (pk != 0) delete[] pk;
  if (ek != 0) delete[] ek;
  if (Ek != 0) delete[] Ek;
  if (akS != 0) delete[] akS;
  if (rhokS != 0) delete[] rhokS;
  if (rhokStar != 0) delete[] rhokStar;
  if (pkStar != 0) delete[] pkStar;
  if (ekStar != 0) delete[] ekStar;
  if (EkStar != 0) delete[] EkStar;
  if (YkStar != 0) delete[] YkStar;
  if (vkStar != 0) delete[] vkStar;
  if (Deltapk != 0) delete[] Deltapk;
  if (zk != 0) delete[] zk;
  if (rho_cIksquare != 0) delete[] rho_cIksquare;
  if (eos != 0) delete[] eos;
  if (Hk0 != 0) delete[] Hk0;
  if (Yk0 != 0) delete[] Yk0;
  if (compactionPk != 0) delete[] compactionPk;
  if (dlambda != 0) delete[] dlambda;
  if (dplast != 0) delete[] dplast;
  if (alphaNull != 0) delete[] alphaNull;
  if (relaxSolidPlast != 0) delete[] relaxSolidPlast;
}

//***********************************************************************

void Tools::uppercase(std::string& string)
{
  for (unsigned int c = 0; c < string.size(); c++){ string[c] = toupper(string[c]); }
}

//***********************************************************************

void Tools::lowercase(std::string& string)
{
	for (unsigned int c = 0; c < string.size(); c++) { string[c] = tolower(string[c]); }
}

//***********************************************************************

void Tools::swap(double &a, double &b)
{
  double buff(a);
  a = b;
  b = buff;
}

//***********************************************************************

double Tools::returnNonZeroValue(double a)
{
  if (a > epsilonAlphaNull || a < - epsilonAlphaNull) {
   return a;
  }
  else {
   return sign(a) * epsilonAlphaNull;
  }
}

//***********************************************************************

double Tools::uselessDouble;