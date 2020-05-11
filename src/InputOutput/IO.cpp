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

//! \file      IO.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      May 03 2018

#include "IO.h"
#include "../Errors.h"

//***********************************************************************

IO::IO(){}

//***********************************************************************

IO::~IO(){}

//***********************************************************************

std::ostream& IO::writeb64Chaine(std::ostream &fluxSortie, char *chaine, int &tailleChaine)
{
  std::string chaineEncodee;
  //int tailleChaine = chaineAEncoder.size();
  int i = 0; int j = 0;
  unsigned char char_array_3[3];
  unsigned char char_array_4[4];
  //const char* chaine = chaineAEncoder.c_str();

  static const std::string base64_chars =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789+/";

  while (tailleChaine--)
  {
    char_array_3[i++] = *(chaine++);
    if (i == 3)
    {
      char_array_4[0] = (char_array_3[0] & 0xfc) >> 2; //Decalage 2 bits a droite
      char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
      char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
      char_array_4[3] = char_array_3[2] & 0x3f;
      for (i = 0; (i <4); i++)
        chaineEncodee += base64_chars[char_array_4[i]];
      i = 0;
    }
  }
  if (i) // Le reste si non multiple de 3
  {
    for (j = i; j < 3; j++)
      char_array_3[j] = '\0';
    char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
    char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
    char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
    char_array_4[3] = char_array_3[2] & 0x3f;
    for (j = 0; (j < i + 1); j++)
      chaineEncodee += base64_chars[char_array_4[j]];
    while ((i++ < 3))
      chaineEncodee += '=';
  }

  return fluxSortie << chaineEncodee.c_str();
}

//***********************************************************************

void IO::copieFichier(std::string file, std::string dossierSource, std::string dossierDestination)
{
  try {
    std::ifstream fichierSource;
    std::ofstream fichierDestination;

    fichierSource.open((dossierSource + file).c_str());
    fichierDestination.open((dossierDestination + file).c_str());
    //std::cout << "copie : " << dossierSource + file << " -> " << dossierDestination + file << std::endl;
    if (!fichierSource) throw ErrorECOGEN("IO::copieFichier : file non trouve \"" + dossierSource + file + "\"", __FILE__, __LINE__);
    if (!fichierDestination) throw ErrorECOGEN("IO::copieFichier : dossier non trouve \"" + dossierDestination + "\"", __FILE__, __LINE__);

    std::string ligne;
    while (getline(fichierSource, ligne)) { fichierDestination << ligne << std::endl; }
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************