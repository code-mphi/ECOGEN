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

#ifndef CELLINTERFACEO2_H
#define CELLINTERFACEO2_H

#include "../Order1/CellInterface.h"

class CellInterfaceO2 : public CellInterface
{
  public:
    CellInterfaceO2();
    CellInterfaceO2(int lvl); //Pour AMR
    virtual ~CellInterfaceO2();

    virtual void allocateSlopes(int& allocateSlopeLocal);
    virtual void computeFlux(double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type = vecPhases);
    /*!< Specific Riemann problem for 2nd order */
    virtual void solveRiemann(double& /*dtMax*/, Limiter& /*globalLimiter*/, Limiter& /*interfaceLimiter*/, Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/, Prim /*type*/ = vecPhases) = 0;

    // -- Cartesian --
    virtual void computeSlopes(Prim /*type*/ = vecPhases) {};

    //Accessors
    virtual Phase* getSlopesPhase(const int& /*phaseNumber*/) const { return nullptr; };
    virtual Mixture* getSlopesMixture() const { return nullptr; };
    virtual Transport* getSlopesTransport(const int& /*numberTransport*/) const { return nullptr; };

    //For AMR method
    virtual void creerCellInterfaceChild() {}; /*!< Creer un child cell interface (non initialize) */
    virtual void creerCellInterfaceChildInterne(const int& /*lvl*/, std::vector<CellInterface*>* /*childrenInternalCellInterfaces*/) {}; /*!< Creer un intern child cell interface (non initialize) */
};

extern Phase** slopesPhasesLocal1;
extern Phase** slopesPhasesLocal2;
extern Mixture* slopesMixtureLocal1;
extern Mixture* slopesMixtureLocal2;
extern double* slopesTransportLocal1;
extern double* slopesTransportLocal2;

#endif // CELLINTERFACEO2_H
