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

#ifndef OUTPUTGLOBALGNU_H
#define OUTPUTGLOBALGNU_H

//! \file      OutputGlobalGNU.h
//! \author    J. Caze
//! \version   1.0
//! \date      January 02 2020

#include "OutputGNU.h"

class OutputGlobalGNU : public OutputGNU
{
public:
	OutputGlobalGNU();
	OutputGlobalGNU(std::string casTest, std::string run, std::string fileName, Input *entree, std::string quantity);
	virtual ~OutputGlobalGNU();

	virtual void prepareSortieSpecifique();

	virtual void ecritSolution(Mesh* mesh, std::vector<Cell*>* cellsLvl);

protected:
	double m_quantity; //!< Physical quantity recorded (mass only available right now)

	void extractTotalQuantity(std::vector<Cell*>* cellsLvl);

	void writeSpecificGnuplotScript();
};

#endif //OUTPUTGLOBALGNU_H