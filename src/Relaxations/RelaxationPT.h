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

#ifndef RELAXATIONPT_H
#define RELAXATIONPT_H

//! \file      RelaxationPT.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      October 16 2018

#include "Relaxation.h"

//! \class     RelaxationPT
//! \brief     Pressure-Temperature relaxation
class RelaxationPT : public Relaxation
{
public:
	RelaxationPT();
  virtual ~RelaxationPT();

  //! \brief     Stiff Pressure-Temperature relaxation method
  //! \details   call for this method computes the mechanical and thermal relaxed state in a given cell. Relaxed state is stored depending on the type enum
  //! \param     cell           cell to relax
  //! \param     numberPhases   number of phases
  //! \param     type           enumeration allowing to relax either state in the cell or second order half time step state
  virtual void stiffRelaxation(Cell *cell, const int &numberPhases, Prim type = vecPhases) const;

  //! \brief     Pressure determination with analytical formulae for 2 phases governed by SG or IG EOS
  //! \details   call for this method determines the pressure using the analytical formulae only valid in the specific case of 2 phases governed by SG or IG EOS
  //! \param     cell           cell to relax
  //! \param     numberPhases   number of phases
  //! \param     type           enumeration allowing to relax either state in the cell or second order half time step state
  //! \return    pressure
  double analyticalPressure(Cell *cell, const int &numberPhases, Prim type = vecPhases) const;

  //! \brief     Temperature calculus with analytical formulae for N phases governed by SG or IG EOS
  //! \details   call for this method computes the temprerature for N phases governed by SG or IG EOS using the known relaxed pressure
  //! \param     pressure       the value of relaxed pressure
  //! \param     cell           cell to relax
  //! \param     numberPhases   number of phases
  //! \param     type           enumeration allowing to relax either state in the cell or second order half time step state
  //! \return    temperature
  double analyticalTemperature(double pressure, Cell *cell, const int &numberPhases, Prim type = vecPhases) const;
};

#endif // RELAXATIONPT_H


