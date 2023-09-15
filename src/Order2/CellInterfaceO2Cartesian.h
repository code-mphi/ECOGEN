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

#ifndef CELLINTERFACEO2CARTESIAN_H
#define CELLINTERFACEO2CARTESIAN_H

#include "CellInterfaceO2.h"

class CellInterfaceO2Cartesian : public CellInterfaceO2
{
public:
    CellInterfaceO2Cartesian();
    CellInterfaceO2Cartesian(int lvl);
    virtual ~CellInterfaceO2Cartesian();

    virtual void allocateSlopes(int& allocateSlopeLocal);
    virtual void computeSlopes(Prim type = vecPhases);
    /*!< Specific Riemann problem for 2nd order */
    virtual void solveRiemann(double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type = vecPhases); 

    //Accessors
    virtual Phase* getSlopesPhase(const int& phaseNumber) const;
    virtual Mixture* getSlopesMixture() const;
    virtual Transport* getSlopesTransport(const int& numberTransport) const;

    //For AMR method
    virtual void creerCellInterfaceChild(); /*!< Creer un child cell interface (non initialize) */
    virtual void creerCellInterfaceChildInterne(const int& lvl, std::vector<CellInterface*>* childrenInternalCellInterfaces); /*!< Creer un intern child cell interface (non initialize) */

protected:
    Phase** m_vecPhasesSlopes;         /*!< Model based array of phase slopes */
    Mixture* m_mixtureSlopes;          /*!< Model based mixture slopes */
    Transport* m_vecTransportsSlopes;  /*!< Model based transport slopes */
};

#endif