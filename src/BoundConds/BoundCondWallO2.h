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

#ifndef BOUNDCONDWALLO2_H
#define BOUNDCONDWALLO2_H

//! \file      BoundCondWallO2.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCondWall.h"


class BoundCondWallO2 : public BoundCondWall
{
public:
  BoundCondWallO2();
  BoundCondWallO2(const BoundCondWallO2& Source, const int lvl = 0); //Constructeur de copie (utile pour AMR)
  BoundCondWallO2(int numPhysique);
  virtual ~BoundCondWallO2();

  virtual void creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces);
  virtual void allocateSlopes(const int &numberPhases, const int &numberTransports, int &allocateSlopeLocal);
  virtual void computeSlopes(const int &numberPhases, const int &numberTransports, Prim type = vecPhases);
  virtual void solveRiemann(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type = vecPhases);

  virtual int whoAmI() const { return 2; };

  //Accesseurs
  virtual Phase* getSlopesPhase(const int &phaseNumber) const;
  virtual Mixture* getSlopesMixture() const;
  virtual Transport* getSlopesTransport(const int &numberTransport) const;

  //Pour methode AMR
  virtual void creerCellInterfaceChild();  /*!< Creer un child cell interface (non initialize) */

protected:
  int m_numberPhases;
  Phase **m_vecPhasesSlopes;         /*!< vecteur des slopes des phases */
  Mixture *m_mixtureSlopes;          /*!< vecteur des slopes de mixture */
  Transport *m_vecTransportsSlopes;	/*!< vecteur des slopes des transports */
private:
};

#endif // BOUNDCONDWALLO2_H
