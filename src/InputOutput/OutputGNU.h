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

#ifndef OUTPUTGNU_H
#define OUTPUTGNU_H

//! \file      OutputGNU.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      May 03 2018

#include "Output.h"
class OutputGNU : public Output
{
public:
  OutputGNU();
  OutputGNU(std::string casTest, std::string run, tinyxml2::XMLElement *element, std::string fileName, Input *entree);
  virtual ~OutputGNU();

  virtual void prepareSortieSpecifique() {};  //Rien a faire pour cette sortie
  virtual void ecritSolution(Mesh *mesh, std::vector<Cell *> *cellsLvl);

protected:

  void ecritScriptGnuplot(const int &dim);

  std::string creationNameFichierGNU(const char* name, int lvl = -1, int proc = -1, int numFichier = -1, std::string nameVariable = "defaut") const;
  void printBlocGnuplot(std::ofstream &fileStream, int &index, const int &dim);

  std::string m_fileNameVisu;
};

#endif //OUTPUTGNU_H