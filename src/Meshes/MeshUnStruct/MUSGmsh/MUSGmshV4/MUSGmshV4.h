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

#ifndef MUSGMSHV4_H
#define MUSGMSHV4_H

//! \file      MUSGmshV4.h
//! \author    J. Caze
//! \version   1.0
//! \date      November 19 2019

#include "../MUSGmsh.h"
#include <map>

class MUSGmshV4 : public MUSGmsh
{
public:
	MUSGmshV4(const std::string &meshFile, const std::string &meshExtension);
	virtual ~MUSGmshV4();

	// --- MeshUnStruct virtual member functions --- 
	virtual void initGeometryMonoCPU(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string computeOrder = "FIRSTORDER");
	//TODO// Reading partitioned mesh file for parallel computation
	virtual void initGeometryParallel(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string computeOrder = "FIRSTORDER") { Errors::errorMessage("reading partitioned Gmsh v4.x mesh file for parallel computation not avalaible in ECOGEN"); };
	virtual void preProcessMeshFileForParallel() { Errors::errorMessage("reading partitioned Gmsh v4.x mesh file for parallel computation not avalaible in ECOGEN"); };

private:
	// --- Member functions ---
	// --- Gmsh v4 related member functions ---
	void readMeshMonoCPU();
	void readElement(const Coord* nodesTable, std::ifstream& meshFile, ElementNS** element, std::vector< std::map<int, int> >& entities, int entityDim, int entityTag, int eltType);
	void readMeshParallel() { Errors::errorMessage("reading partitioned Gmsh v4.x mesh file for parallel computation not avalaible in ECOGEN"); };
};

#endif // MUSGMSHV4_H