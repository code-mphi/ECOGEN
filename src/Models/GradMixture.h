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

#ifndef GRADMIXTURE_H
#define GRADMIXTURE_H

#include "../Tools.h"

class Cell;
class Mixture;

#include "../Maths/Coord.h"

//! \class     GradMixture
//! \brief     Mixture variable gradients. Stored for 2nd-order computation on unstructured mesh (O2 NS)

class GradMixture
{
public:
  GradMixture();
  virtual ~GradMixture();

  virtual void initializeGradientVectors() { Errors::errorMessage("initializeGradientVectors not available for required GradMixture"); };

  virtual void computeGradients(Cell* cell);
  virtual void computeDistanceGradientScalarProduct(Coord const& /*distance*/, Mixture* /*mixture*/) const { Errors::errorMessage("computeDistanceGradientScalarProduct not available for required GradMixture"); };
  virtual void limitGradients(const Mixture& /*gradientLimiter*/) { Errors::errorMessage("limitGradients not available for required GradMixture"); };

  // -- O2 parallel --
  virtual int numberOfTransmittedGradients() const { Errors::errorMessage("numberOfTransmittedGradients not available for required GradMixture"); return Errors::defaultInt; };
  virtual void getBufferGradients(double * buffer, int& counter);
  virtual void fillBufferGradients(double * buffer, int& counter);

protected:
  void initializeGradsVariablesNamesNumerators();
  
  std::vector<Coord> m_grads;    //!< Vector of gradients of the mixture flow variables
};

extern std::vector<Variable> variableNamesMixture;   //!< Variable names of the corresponding gradients (for mixture)
extern std::vector<int> numeratorMixture;            //!< Default number (not used for mixture variables)

#endif