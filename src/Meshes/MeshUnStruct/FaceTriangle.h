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

#ifndef FACETRIANGLE_H
#define FACETRIANGLE_H

//! \file      FaceTriangle.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include "FaceNS.h"

class FaceTriangle : public FaceNS
{
public:
  FaceTriangle(const int &numNoeud1, const int &numNoeud2, const int &numNoeud3, int tri=1);
  virtual ~FaceTriangle();

private:
  virtual void computeSurface(const Coord *noeuds);
  virtual void computeRepere(const Coord *noeuds, const int &numNoeudAutre, ElementNS *elementVoisin);

  static const int NOMBRENOEUDS;
};

#endif // FACETRIANGLE_H
