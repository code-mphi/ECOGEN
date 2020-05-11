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

//! \file      Tools.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      December 6 2018

#include "Tools.h"

Tools *TB;

//***********************************************************************

Tools::Tools() : ak(0), rhok(0), pk(0), akS(0), rhokS(0), eos(0)
{}

//***********************************************************************

Tools::Tools(const int &numberPhases)
{
  m_numberPhases = numberPhases;
  ak = new double[numberPhases];
  Yk = new double[numberPhases];
  rhok = new double[numberPhases];
  pk = new double[numberPhases];
  akS = new double[numberPhases];
  rhokS = new double[numberPhases];
  rhokStar = new double[numberPhases];
  pkStar = new double[numberPhases];
  ekStar = new double[numberPhases];
  YkStar = new double[numberPhases];
  vkStar = new double[numberPhases];
  Deltaek = new double[numberPhases];
  eos = new Eos*[numberPhases];

  Hk0 = new double[numberPhases];
  Yk0 = new double[numberPhases];

}

//***********************************************************************

Tools::~Tools()
{
  if (ak!=0) delete[] ak;
  if (Yk != 0) delete[] Yk;
  if (rhok != 0) delete[] rhok;
  if (pk != 0) delete[] pk;
  if (akS != 0) delete[] akS;
  if (rhokS != 0) delete[] rhokS;
  if (rhokStar != 0) delete[] rhokStar;
  if (pkStar != 0) delete[] pkStar;
  if (ekStar != 0) delete[] ekStar;
  if (YkStar != 0) delete[] YkStar;
  if (vkStar != 0) delete[] vkStar;
  if (Deltaek != 0) delete[] Deltaek;
  if (eos != 0) delete[] eos;

  if (Hk0 != 0) delete[] Hk0;
  if (Yk0 != 0) delete[] Yk0;
}

//***********************************************************************

void Tools::uppercase(std::string &string)
{
  for (unsigned int c = 0; c < string.size(); c++){ string[c] = toupper(string[c]); }
}

//***********************************************************************

void Tools::lowercase(std::string& string)
{
	for (unsigned int c = 0; c < string.size(); c++) { string[c] = tolower(string[c]); }
}

//***********************************************************************

double Tools::pi()
{
  //return acos(-1.);
  return 3.14159;
}

//***********************************************************************