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

#ifndef CELLO2CARTESIAN_H
#define CELLO2CARTESIAN_H

#include "CellO2.h"

class CellInterface;

class CellO2Cartesian : public CellO2
{
public:
  CellO2Cartesian();
  CellO2Cartesian(int lvl);
  virtual ~CellO2Cartesian();

  virtual void computeLocalSlopes(CellInterface& cellInterfaceRef, Limiter& globalLimiter, Limiter& interfaceLimiter,
      Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter,
      double& alphaCellAfterOppositeSide, double& alphaCell, double& alphaCellOtherInterfaceSide, double& epsInterface);
      
  virtual void computeLocalSlopesLimite(CellInterface& cellInterfaceRef, Limiter& globalLimiter,
      Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter,
      double& epsInterface);

  //Pour methode AMR
  virtual void createChildCell(const int& lvl); /*!< Creer une cell enfant (non initializee) */

  //Pour methodes ordre 2 parallele
  virtual void getBufferSlopes(double* /*buffer*/, int& /*counter*/, const int& /*lvl*/) { Errors::errorMessage("getBufferSlopes not available for CellO2Cartesian"); };
  virtual void fillBufferSlopes(double* buffer, int& counter, const int& lvl, const int& neighbour) const;
};

#endif