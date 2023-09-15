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

#ifndef GRADPHASE_H
#define GRADPHASE_H

#include "../Maths/Coord.h"
#include "../Tools.h"

class Cell;
class Phase;

//! \class     GradPhase
//! \brief     Phase variable gradients. Stored for 2nd-order computation on unstructured mesh (O2 NS)

class GradPhase
{
public:
  GradPhase();
  virtual ~GradPhase();

  virtual void initializeGradientVectors() { Errors::errorMessage("initializeGradientVectors not available for required GradPhase"); };

  virtual void computeGradients(Cell* cell, int const& phase);
  virtual void computeDistanceGradientScalarProduct(Coord const& /*distance*/, Phase* /*phase*/) const { Errors::errorMessage("computeDistanceGradientScalarProduct not available for required GradPhase"); };
  virtual void limitGradients(const Phase& /*gradientLimiter*/) { Errors::errorMessage("limitGradients not available for required GradPhase"); };

  // -- O2 parallel --
  virtual int numberOfTransmittedGradients() const { Errors::errorMessage("numberOfTransmittedGradients not available for required GradPhase"); return Errors::defaultInt; };
  virtual void getBufferGradients(double* buffer, int& counter);
  virtual void fillBufferGradients(double* buffer, int& counter);

protected:
  void initializeGradsVariablesNamesNumerators();

  std::vector<Coord> m_grads;    //!< Vector of gradients of the phase flow variables
};

extern std::vector<Variable> variableNamesPhases;       //!< Variable names of the corresponding gradients (for phases)
extern std::vector<std::vector<int>> numeratorPhases;   //!< Numerator of the computed phase for gradients (for phases)

#endif