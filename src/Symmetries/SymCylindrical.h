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

#ifndef SYMCYLINDRICAL_H
#define SYMCYLINDRICAL_H

//! \file      SymmetryCylindrical.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "Symmetry.h"

//! \class     SymmetryCylindrical
//! \brief     General class for cylindrical symmetry
class SymCylindrical : public Symmetry
{
public:
  SymCylindrical();
  SymCylindrical(tinyxml2::XMLElement *element, std::string nameFile = "Unknown file");
  virtual ~SymCylindrical();

  virtual void addSymmetricTerms(Cell *cell, const int &numberPhases, Prim type = vecPhases);
  virtual void addSymmetricTermsAddPhys(Cell *cell, const int &numberPhases, AddPhys &addPhys);
};

#endif //SYMCYLINDRICAL_H