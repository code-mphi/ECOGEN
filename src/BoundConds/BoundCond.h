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

#ifndef BOUNDCOND_H
#define BOUNDCOND_H

#include <iostream>

#include "../Order1/CellInterface.h"
#include "../Order2/CellInterfaceO2.h" //Ajouter pour l'AMR, a priori ne pose pas de probleme
#include "../libTierces/tinyxml2.h"
#include "../Errors.h"
#include "../Tools.h"

class BoundCond : public CellInterface
{
  public:
    BoundCond();
    BoundCond(int numPhysique);
    BoundCond(const BoundCond& Source, const int& lvl = 0);
    virtual ~BoundCond();

    virtual void createBoundary(TypeMeshContainer<CellInterface*>& /*cellInterfaces*/) { Errors::errorMessage("Impossible to create the limit in createBoundary"); };
    virtual void createBoundary(TypeMeshContainer<CellInterface*>& /*cellInterfaces*/, std::string /*ordreCalcul*/) { Errors::errorMessage("Impossible to create the limit in createBoundary"); };
    virtual void initialize(Cell* cellLeft, Cell* /*cellRight*/);

    virtual void computeFlux(const int& numberPhases, const int& numberTransports, double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type = vecPhases);
    virtual void computeFluxAddPhys(const int& numberPhases, AddPhys& addPhys);
    virtual void solveRiemann(const int& numberPhases, const int& numberTransports, double& dtMax, Limiter& /*globalLimiter*/, Limiter& /*interfaceLimiter*/, Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/, Prim type = vecPhases);
    virtual void addFlux(const int& /*numberPhases*/, const int& /*numberTransports*/, const double& /*coefAMR*/) {}; //Since it is a boundary there is nothing to add to 'right' cell
    virtual void solveRiemannBoundary(Cell& /*cellLeft*/, const int& /*numberPhases*/, const double& /*dxLeft*/, double& /*dtMax*/) { Errors::errorMessage("Warning solveRiemannLimite not provided for boundary used"); };
    virtual void solveRiemannTransportBoundary(Cell& /*cellLeft*/, const int& /*numberTransports*/) const { Errors::errorMessage("Warning solveRiemannTransportLimite not provided for boundary used"); };

    virtual int whoAmI() const { Errors::errorMessage("whoAmI not available for boundary used"); return 0; };
    virtual int whoAmIHeat() const { Errors::errorMessage("whoAmIHeat not available for boundary used"); return ADIABATIC; };
    virtual void printInfo(){};

    virtual const int& getNumPhys() const { return m_numPhysique; };
    virtual double getMassflow() const;
    virtual double getPowerFlux() const;
    virtual double getBoundaryHeatQuantity() const { Errors::errorMessage("getBoundaryHeatQuantity not available for boundary used"); return 0.; }

    //Pour methode AMR
    virtual void computeXi(const double& /*criteriaVar*/, const bool& /*varRho*/, const bool& /*varP*/, const bool& /*varU*/, const bool& /*varAlpha*/) {};
    virtual void computeFluxXi();
    virtual void raffineCellInterfaceExterne(const int& nbCellsY, const int& nbCellsZ, const double& dXParent, const double& dYParent, const double& dZParent, Cell* cellRef, const int& dim);
    virtual void deraffineCellInterfaceExterne(Cell* cellRef);

  protected:
    int m_numPhysique;
    double m_massflow;     //!< Total massflow through boundary (kg.s-1)
    double m_powerFlux;    //!< Power flux through boundary (W)

  private:
};

#endif // BOUNDCOND_H
