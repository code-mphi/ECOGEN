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

#ifndef CELLO2NS_H
#define CELLO2NS_H

#include "CellO2.h"

class CellO2NS : public CellO2
{
public:
  CellO2NS();
  virtual ~CellO2NS();

  virtual void allocateSecondOrderBuffersAndGradientVectors(Phase** phases, Mixture* mixture);
  virtual void allocate(const std::vector<AddPhys*>& addPhys);
  virtual void computeGradientsO2();
  virtual void limitGradientsO2(Limiter& globalLimiter);
  virtual void computeLocalSlopes(CellInterface& cellInterfaceRef);
  
  // For 2nd order with parallel
  virtual GradPhase* getGradPhase(const int& phaseNumber) const;
  virtual GradMixture* getGradMixture() const;
  virtual GradTransport* getGradTransport(const int& transportNumber) const;
  virtual void getBufferSlopes(double* /*buffer*/, int& /*counter*/, const int& /*lvl*/) { Errors::errorMessage("getBufferSlopes not available for CellO2NS"); };
  virtual void fillBufferSlopes(double* buffer, int& counter, const int& /*lvl*/, const int& /*neighbour*/) const;

protected:
  GradPhase** m_gradPhase;
  GradMixture* m_gradMixture;
  GradTransport* m_gradTransport;
};

extern Phase** buffPhasesMin;         //!< Stores minimum phases from neighbors of a cell
extern Phase** buffPhasesMax;         //!< Stores maximum phases from neighbors of a cell
extern Mixture* buffMixtureMin;       //!< Stores minimum mixture from neighbors of a cell
extern Mixture* buffMixtureMax;       //!< Stores maximum mixture from neighbors of a cell
extern double* buffTransportMin;      //!< Stores minimum transport from neighbors of a cell
extern double* buffTransportMax;      //!< Stores maximum transport from neighbors of a cell

#endif