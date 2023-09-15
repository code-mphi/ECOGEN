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

#ifndef GRADMIXEULER_H
#define GRADMIXEULER_H

#include "../GradMixture.h"

//! \class     GradMixEuler
//! \brief     Mixture variable gradients for Euler model. Stored for 2nd-order computation on unstructured mesh (O2 NS)
//! \details   No variables stored here since this model is for single-phase flows.

class GradMixEuler : public GradMixture
{
public:
  GradMixEuler() {};
  virtual ~GradMixEuler() {};

  virtual void initializeGradientVectors() {};

  virtual void computeGradients(Cell* /*cell*/) {};
  virtual void computeDistanceGradientScalarProduct(Coord const& /*distance*/, Mixture* /*mixture*/) const {};
  virtual void limitGradients(const Mixture& /*gradientLimiter*/) {}; 

  // -- O2 parallel --
  virtual int numberOfTransmittedGradients() const { return 0; };
  virtual void getBufferGradients(double * /*buffer*/, int& /*counter*/) {};
  virtual void fillBufferGradients(double * /*buffer*/, int& /*counter*/) {};
};

#endif