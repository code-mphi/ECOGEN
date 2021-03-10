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

#ifndef BOUNDCONDSYMMETRYO2_H
#define BOUNDCONDSYMMETRYO2_H

#include "BoundCondWallO2.h"


class BoundCondSymmetryO2 : public BoundCondWallO2
{
public:
  BoundCondSymmetryO2(const BoundCondSymmetryO2& Source, const int& lvl = 0); //Copy ctor (useful for AMR)
  BoundCondSymmetryO2(int numPhysique);
  virtual ~BoundCondSymmetryO2();

  virtual void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces);

  virtual int whoAmI() const { return SYMMETRY; };

  //For AMR method
  virtual void creerCellInterfaceChild();  /*!< Create a child cell interface (not initialized) */

protected:

private:
};

#endif // BOUNDCONDSYMMETRYO2_H
