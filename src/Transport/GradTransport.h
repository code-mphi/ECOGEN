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

#ifndef GRADTRANSPORT_H
#define GRADTRANSPORT_H

#include "../Maths/Coord.h"
#include "../Tools.h"

//! \class     GradTransport
//! \brief     Transport variable gradients. Stored for 2nd-order computation on unstructured mesh (O2 NS)

class Cell;

class GradTransport
{
public:
  GradTransport();
  virtual ~GradTransport();

  virtual void initializeGradientVectors();

  virtual void computeGradient(Cell* cell, const int& numTransport);
  virtual void computeDistanceGradientScalarProduct(Coord const& distance, double &transport) const;
  virtual void limitGradients(const double& gradientLimiter);

  virtual int numberOfTransmittedGradients() const;
  virtual void getBufferGradients(double * buffer, int& counter);
  virtual void fillBufferGradients(double * buffer, int& counter);

protected:
  std::vector<Coord> m_grads;    //!< Gradient of the transport variable (size of 1)
};

extern std::vector<Variable> variableNamesTransports;        //!< Variable names of the corresponding gradients (for transports)
extern std::vector<std::vector<int>> numeratorTransports;    //!< Numerator of the computed transport for gradients

#endif