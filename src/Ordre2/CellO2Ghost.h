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

#ifndef CELLO2GHOST_H
#define CELLO2GHOST_H

//! \file      CellO2Ghost.h
//! \author    F. Petitpas, K.Schmidmayer
//! \version   1.0
//! \date      July 30 2018

#include "CellO2.h"

class CellO2Ghost : public CellO2
{
public:
	CellO2Ghost();
	CellO2Ghost(int lvl); //Pour AMR
	virtual ~CellO2Ghost();

	virtual void allocate(const int &numberPhases, const int &numberTransports, const std::vector<AddPhys*> &addPhys, Model *model);
	virtual void computeLocalSlopes(const int &numberPhases, const int &numberTransports, CellInterface &bordRef, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, double &alphaCellAfterOppositeSide, double &alphaCell, double &alphaCellOtherInterfaceSide, double &epsInterface);
	virtual void createChildCell(const int &num, const int &lvl);
	virtual void getBufferSlopes(double *buffer, int &counter);
	virtual void getBufferSlopesAMR(double *buffer, int &counter, const int &lvl);
	virtual bool isCellO2Ghost() const { return true; };

protected:
	Phase **m_vecPhasesSlopesGhost;         /*!< pour stocker les slopes des phases */
	Mixture *m_mixtureSlopesGhost;          /*!< pour stocker les slopes de mixture */
	double *m_vecTransportsSlopesGhost;	    /*!< pour stocker les slopes des transports */
  double m_alphaCellAfterOppositeSide;    /*!< pour stocker la fraction volumique apres la cellule fantome */
	
private:
};

#endif // CELLO2GHOST_H
