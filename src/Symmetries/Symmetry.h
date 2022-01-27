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

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "../libTierces/tinyxml2.h"
#include "../Errors.h"
#include "../Tools.h"

class Symmetry; //Predeclaration of the class Symmetry to include Cell.h
#include "../Order1/Cell.h"

//! \class     Symmetry
//! \brief     General class for axial symmetries
//! \details   This is a pure virtual class: can not be instantiated
class Symmetry
{
public:
  Symmetry();
  virtual ~Symmetry();
  //! \brief     Add the symmetric terms for the cylindrical or spherical symmetry assumption
  //! \param     cell           cell to add the terms
  //! \param     type           enumeration allowing to correct either state in the cell or second order half time step state
  virtual void addSymmetricTerms(Cell* /*cell*/, Prim /*type*/ = vecPhases) {};
  //! \brief     Add the additional-physics, symmetric terms for the cylindrical or spherical symmetry assumption
  //! \param     cell           cell to add the terms
  //! \param     addPhys        additional-physics object to call the corresponding symmetry subroutine
  virtual void addSymmetricTermsAddPhys(Cell* /*cell*/, AddPhys& /*addPhys*/) {};

protected:
  int m_radialAxis;   //!< Name of the radial axis for the axi-symmetry
};

#endif // SYMMETRY_H