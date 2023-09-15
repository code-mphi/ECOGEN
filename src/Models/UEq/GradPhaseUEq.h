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

#ifndef GRADPHASEUEQ_H
#define GRADPHASEUEQ_H

#include "../GradPhase.h"
#include "PhaseUEq.h"

class Phase;

//! \class     GradPhaseUEq
//! \brief     Phase variable gradients for UEq model. Stored for 2nd-order computation on unstructured mesh (O2 NS)

class GradPhaseUEq : public GradPhase
{
public:
  GradPhaseUEq();
  virtual ~GradPhaseUEq();

  virtual void initializeGradientVectors();

  virtual void computeDistanceGradientScalarProduct(Coord const& distance, Phase* phase) const;
  virtual void limitGradients(const Phase& gradientLimiter);

  // O2 parallel
  virtual int numberOfTransmittedGradients() const;

protected:
  //! \brief     Enumeration for the phase flow variables, specific to UEq
  enum VarLocal { alpha, density, pressure };
};

#endif