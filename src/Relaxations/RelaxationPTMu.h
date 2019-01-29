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

#ifndef RELAXATIONPTMU_H
#define RELAXATIONPTMU_H

//! \file      RelaxationPTMu.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      October 16 2018

#include "Relaxation.h"

//! \class     RelaxationPTMU
//! \brief     Pressure-Temperature-Chemical Potential relaxation / Phase change
class RelaxationPTMu : public Relaxation
{
public:
	RelaxationPTMu();
	//! \brief     Relaxation constructor from a XML format reading
	//! \details   Reading data from XML file under the following format:
	//!            ex: <dataPTMu liquid="SG_waterLiq.xml" vapor="IG_waterVap.xml"/>
	//! \param     element          XML element to read for source term
	//! \param     fileName         string name of readed XML file
	RelaxationPTMu(tinyxml2::XMLElement *element, std::string fileName = "Unknown file");
	virtual ~RelaxationPTMu();

	//! \brief     Stiff Thermo-Chemical relaxation method
	//! \details   call for this method computes the thermodyanmical equilibrium state in a given cell for a liquid and its vapor. Relaxed state is stored depending on the type enum
	//! \param     cell           cell to relax
	//! \param     numberPhases   number of phases
	//! \param     type           enumeration allowing to relax either state in the cell or second order half time step state
	virtual void stiffRelaxation(Cell *cell, const int &numberPhases, Prim type = vecPhases) const;

private:
	int m_liq;   //!< Liquid phase number for phase change
	int m_vap;   //!< Vapor phase number for phase change
};

#endif // RELAXATIONPTMU_H
