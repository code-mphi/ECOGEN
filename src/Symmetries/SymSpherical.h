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

#ifndef SYMSPHERICAL_H
#define SYMSPHERICAL_H

//! \file      SymmetrySpherical.h
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "Symmetry.h"

//! \class     SymmetrySpherical
//! \brief     General class for spherical symmetry
class SymSpherical : public Symmetry
{
public:
  SymSpherical();
  SymSpherical(tinyxml2::XMLElement *element, std::string nameFile = "Unknown file");
  virtual ~SymSpherical();

  virtual void addSymmetricTerms(Cell *cell, const int &numberPhases, Prim type = vecPhases);
  virtual void addSymmetricTermsAddPhys(Cell *cell, const int &numberPhases, AddPhys &addPhys) { Errors::errorMessage("addSymmetricTermsAddPhys not implemented for spherical symmetry with additional physics"); };
};

#endif //SYMSPHERICAL_H