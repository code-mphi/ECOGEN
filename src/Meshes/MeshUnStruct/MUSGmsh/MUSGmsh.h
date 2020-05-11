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

#ifndef MUSGMSH_H
#define MUSGMSH_H

//! \file      MUSGmsh.h
//! \author    J. Caze
//! \version   1.0
//! \date      November 20 2019

#include <string>
#include <sstream>
#include "../../MeshUnStruct.h"
#include "HeaderElements.h"

class MUSGmsh : public MeshUnStruct
{
public:
  MUSGmsh(const std::string &meshFile, const std::string &meshExtension);
  virtual ~MUSGmsh();

  static std::string readVersion(const std::string &meshFile);

  // --- MeshUnStruct virtual member functions --- 
  virtual void initGeometryMonoCPU(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string computeOrder = "FIRSTORDER") = 0;
  virtual void initGeometryParallel(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string computeOrder = "FIRSTORDER") = 0;
  virtual void preProcessMeshFileForParallel() = 0;

private:
};

#endif // MUSGMSH_H