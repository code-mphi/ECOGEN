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

#ifndef BOUNDCONDSYMMETRYO2_H
#define BOUNDCONDSYMMETRYO2_H

//! \file      BoundCondSymmetryO2.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCondWallO2.h"


class BoundCondSymmetryO2 : public BoundCondWallO2
{
public:
  BoundCondSymmetryO2();
  BoundCondSymmetryO2(const BoundCondSymmetryO2& Source, const int lvl = 0); //Constructeur de copie (utile pour AMR)
  BoundCondSymmetryO2(int numPhysique);
  virtual ~BoundCondSymmetryO2();

  virtual void creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces);

  virtual int whoAmI() const { return 6; };

  //Pour methode AMR
  virtual void creerCellInterfaceChild();  /*!< Creer un child cell interface (non initialize) */

protected:

private:
};

#endif // BOUNDCONDSYMMETRYO2_H
