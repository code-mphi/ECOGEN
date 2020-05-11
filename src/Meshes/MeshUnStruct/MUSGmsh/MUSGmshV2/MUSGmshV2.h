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

#ifndef MUSGMSHV2_H
#define MUSGMSHV2_H

//! \file      MUSGmshV2.h
//! \author    J. Caze
//! \version   1.0
//! \date      November 19 2019

#include "../MUSGmsh.h"

class MUSGmshV2 : public MUSGmsh
{
public:
  MUSGmshV2(const std::string& meshFile, const std::string& meshExtension, bool switchTags = false);
  virtual ~MUSGmshV2();

  // --- MeshUnStruct virtual member functions --- 
  virtual void initGeometryMonoCPU(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string computeOrder = "FIRSTORDER");
  virtual void initGeometryParallel(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string computeOrder = "FIRSTORDER");
  virtual void preProcessMeshFileForParallel();

private:
  // --- Gmsh v2 related member functions ---
  void readMeshMonoCPU(std::vector<ElementNS*>** neighborNodes);
  void readElement(const Coord *nodesTable, std::ifstream &meshFile, ElementNS **element);
  void readMeshParallel();

  bool m_switchTags;
};


#endif // MUSGMSHV2_H