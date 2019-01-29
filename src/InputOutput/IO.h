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

#ifndef IO_H
#define IO_H

//! \file      IO.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      May 03 2018

#include <fstream>
#include <algorithm>
#include <sstream>

class IO
{
public:
  IO();
  virtual ~IO();

  //Templates pour convertir donnees depuis ou vers string
  //-----------------------------------------------------
  template <typename T>
  static std::string toString(const T& i)
  {
	  std::ostringstream stream;
	  stream << i;
	  return stream.str();
  }

  template <typename T>
  static T fromString(const std::string& str)
  {
	  std::istringstream stream(str);
	  T tmp;
	  stream >> tmp;
	  return tmp;
  }

  //Templates Format Binaire pour le Legacy VTK
  //-------------------------------------------

  //Definition de template pour print au format binary
  template <typename T>
  static std::ostream& write(std::ostream &fluxSortie, T &value)
  {
    //Swap Little <-> Big endian eventuel
    IO::endswap(&value);
    return fluxSortie.write(reinterpret_cast<char*>(&value), sizeof(T));
  };

  // //Definition de template pour lecture au format binary
  // template <typename T>
  // static std::istream& read(std::istream &fluxEntree, T& value)
  // {
  //   return fluxEntree.read(reinterpret_cast<char*>(&value), sizeof(T));
  //   //Swap Little <-> Big endian eventuel
  //   IO::endswap(&value);
  //   //Bug ici je pense, a voir...
  // };

  //Templates Format Binaire Base64 pour le XML VTK
  //-----------------------------------------------

  //Template pour ajouter n importe quelle type de donnee a une chaine de caractere
  template <typename T>
  static void ajouteAlaChaine(char *chaine, int &taille, T &value)
  {
    char *conversionChaine = reinterpret_cast<char*>(&value);
    for (int octet = 0; octet < sizeof(value); octet++)
    { chaine[taille++] = conversionChaine[octet]; }  
  }

  //Definition de template pour print d un number au format binary base64
  template <typename T>
  static std::ostream& writeb64(std::ostream &fluxSortie, T &value)
  {
    //Swap Little <-> Big endian eventuel ne marche pas ??
    //IO::endswap(&value);
    //Encodage Base64
    char* chaine = reinterpret_cast<char*>(&value);
    int tailleChaine = sizeof(value);
    return IO::writeb64Chaine(fluxSortie, chaine, tailleChaine);
  };

  // //ATTENTION !!!!!!!!!!!Lecture non Fonctionnelle !!!!!!!!!!!!!!
  // //Definition de template pour lecture au format binary base64
  // template <typename T>
  // static std::istream& readb64(std::istream &fluxEntree, T& value)
  // {
  //   return fluxEntree.read(reinterpret_cast<char*>(&value), sizeof(T));
  //   //Swap Little <-> Big endian eventuel
  //   IO::endswap(&value);
  // };
  // //ATTENTION !!!!!!!!!!!Lecture non Fonctionnelle !!!!!!!!!!!!!!

  static std::ostream& writeb64Chaine(std::ostream &fluxSortie, char *chaineAEncoder, int &taille);

  static void copieFichier(std::string file, std::string dossierSource, std::string dossierDestination);

private:

  //Swap Little <-> Big Endian
  template <typename T>
  static void endswap(T *objp)
  {
    unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
    std::reverse(memp, memp + sizeof(T));
  }
  
};

#endif // IO_H


/*
----------------Mode d emploi binary-----------------
exemple avec un entier code sur 4 octets

int value = 256 487 423
Sa representation binary avec les bits de poids fort a gauche est :
0000 1111    0100 1001    1010 1111    1111 1111
En hexadecimal cela devient:
 0    F       4    9       A    F       F    F         soit 0xF49AFFF
Cet entier peut être stocke dans un tableau de 4 caracteres de 1 octet chacun par:
char *tableau = reinterpret_cast<char *>(value)
Chaque octet correspond alors a un caractere de la table ASCII (sur 8 bits)

----------BASE64-------------

Base64 : L'objectif est d'encoder chaque groupe de 24 bits successif par une chaîne de 4 caracteres simples.
Pour passer en Base64 soit sur des multiples de 24 bits (4 x 6bits), on decoupe
la representation binary en paquet de 6 bits. Dans notre cas, cela devient:
000011 110100 100110 101111    111111 11
le dernier etant incomplet on complete avec des 0 :
000011 110100 100110 101111    111111 110000

On associe ensuite a chaque code binary sur 6 bits un caractere simple parmis les 64 suivants:
"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

Ainsi, la transcription en Base64 donnera:
  D       0      m      v         /      w

Tout ceci est stocke a nouveau dans des chaînes de caractere (8 bits par caractere),
on complete alors par deux bits nuls devant chaque jeu de 6 bits :
00000011 00110100 00100110 00101111    00111111 00110000
on complete egalement pour obtenir un multiple de 4 caracteres avec le caractere '=':

Ainsi la transcription final de notre entier en Base64 sera :
  D0mv/w==

qui aura la representation binary final (Pas sur pour le = si c'est sa representation ASCII ou pas)
00000011 00110100 00100110 00101111    00111111 00110000 00111101 00111101

------------------------------------------------------

zone de tests inversion b64
CwsLCwsLCws=    //chaine binary 64 d'entiers sur 8 bits
000010 110000 101100 001011 000010 110000 101100 001011 000010 110000 101100  //Traduction binary
00001011 00001011 00001011 00001011 00001011 00001011 00001011 00001011 (00)    //Regroupement par 8
11 11 11 11 11 11 11 11 //Traduction entier sur 8 bits

CQkJCQkJCQkJCQkJCQkJCQ==
000010 010000 100100 001001 000010 010000 100100 001001 000010 010000 100100 001001 000010 010000 100100 001001 000010 010000 100100 001001 000010 010000
00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 00001001 0000
9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 


*/