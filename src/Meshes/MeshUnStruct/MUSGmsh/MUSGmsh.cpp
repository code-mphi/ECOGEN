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

//! \file      MUSGmsh.h
//! \author    J. Caze
//! \version   1.0
//! \date      November 20 2019

#include "MUSGmsh.h"

//***********************************************************************

MUSGmsh::MUSGmsh(const std::string &meshFile, const std::string &meshExtension) : MeshUnStruct(meshFile,meshExtension)
{}

//***********************************************************************

MUSGmsh::~MUSGmsh() {}

//***********************************************************************

std::string MUSGmsh::readVersion(const std::string &meshFile)
{
	try {
		std::string pathMeshFile("./libMeshes/" + meshFile);
		std::string currentLine;
		std::ifstream mesh(pathMeshFile.c_str(), std::ios::in);
		if (!mesh) { throw ErrorXML("mesh file not found : " + pathMeshFile, __FILE__, __LINE__); }

		getline(mesh, currentLine); // read $MeshFormat

		std::stringstream lineToTreat;
		std::string fileVersion;
		getline(mesh, currentLine);
		lineToTreat << currentLine;
		lineToTreat >> fileVersion;
		mesh.close();
		return fileVersion;
	}
	catch (ErrorXML&) { throw; }
}

//***********************************************************************
