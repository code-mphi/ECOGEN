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

#ifndef ADDPHYS_H
#define ADDPHYS_H

//! \file      AddPhys.h
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

class AddPhys; //Predeclaration of the AddPhys class to include Cell.h and CellBoundary.h

#include "../Errors.h"
#include "../Cell.h"
#include "../CellInterface.h"
#include "QuantitiesAddPhys.h"
#include "../Parallel.h"

//! \class     AddPhys
//! \brief     General class for additional physics
//! \details   This is a pure virtual class: can not be instantiated
class AddPhys
{
  public:
    AddPhys();
    virtual ~AddPhys();
    
    //! \brief     Add the quantities for the additional physic
    //! \param     cell                 corresponding cell
    virtual void addQuantityAddPhys(Cell *cell) { Errors::errorMessage("addQuantityAddPhys not implemented for used additional physic"); };

    virtual std::string whoAmI() const { Errors::errorMessage("whoAmI not implemented for used additional physic"); return 0; };

    //! \brief     Compute and send back mass energie linked to the physic (0 if no linked energy)
    //! \param     QPA                  corresponding additional physic quantities
    virtual double computeEnergyAddPhys(QuantitiesAddPhys* QPA) { return 0.;  };
    //! \brief     Compute the additional physic flux between two cells
    //! \param     cellBound            cell boundary
    //! \param     numberPhases         number of phases
    void computeFluxAddPhys(CellInterface *cellBound, const int &numberPhases);
    //! \brief     Compute the additional physic flux at a domain boundary
    //! \param     cellBound            cell boundary
    //! \param     numberPhases         number of phases
    void computeFluxAddPhysBoundary(CellInterface *cellBound, const int &numberPhases);
    //! \brief     Add the non-conservative terms of the additional physic in a cell
    //! \param     cell                 cell
    //! \param     numberPhases         number of phases
    void addNonConsAddPhys(Cell *cell, const int &numberPhases);
    //! \brief     Solve the additional physic flux between two cells
    //! \param     cellBound            cell boundary
    //! \param     numberPhases         number of phases
    virtual void solveFluxAddPhys(CellInterface *cellBound, const int &numberPhases) { Errors::errorMessage("solveFluxAddPhys not implemented for used additional physic"); };
    //! \brief     Solve the additional physic flux at a domain boundary
    //! \param     cellBound            cell boundary
    //! \param     numberPhases         number of phases
    virtual void solveFluxAddPhysBoundary(CellInterface *cellBound, const int &numberPhases) { Errors::errorMessage("solveFluxAddPhysLimite not implemented for used additional physic"); };
    //! \brief     Add the additional physic flux between two cells at the corresponding cell
    //! \param     cellBound            cell boundary
    //! \param     numberPhases         number of phases
    //! \param     coefAMR              Adaptive Mesh Refinement coefficient for numerical scheme purposes
    virtual void addFluxAddPhys(CellInterface *cellBound, const int &numberPhases, const double &coefAMR);
    //! \brief     Subtract the additional physic flux between two cells at the corresponding cell
    //! \param     cellBound            cell boundary
    //! \param     numberPhases         number of phases
    //! \param     coefAMR              Adaptive Mesh Refinement coefficient for numerical scheme purposes
    virtual void subtractFluxAddPhys(CellInterface *cellBound, const int &numberPhases, const double &coefAMR);
    //! \brief     Add the non-conservative terms of the corresponding additional physic in a cell
    //! \param     cell                 cell
    //! \param     numberPhases         number of phases
    virtual void addNonCons(Cell *cell, const int &numberPhases) { Errors::errorMessage("addNonCons not implemented for used additional physic"); };
    //! \brief     Add the symmetrical terms of the corresponding additional physic in a cell (when radial axe is on X)
    //! \param     cell                 cell
    //! \param     numberPhases         number of phases
    virtual void addSymmetricTermsRadialAxeOnX(Cell *cell, const int &numberPhases) { Errors::errorMessage("addSymmetricTerms not implemented for used additional physic"); };
    //! \brief     Add the symmetrical terms of the corresponding additional physic in a cell (when radial axe is on Y)
    //! \param     cell                 cell
    //! \param     numberPhases         number of phases
    virtual void addSymmetricTermsRadialAxeOnY(Cell *cell, const int &numberPhases) { Errors::errorMessage("addSymmetricTerms not implemented for used additional physic"); };

    //! \brief     Reinitialize the color function for the surface-tension terms with the volume fraction equation
    //! \param     cellsLvl             cells vector for every level
    //! \param     lvl                  level
		virtual void reinitializeColorFunction(std::vector<Cell *> *cellsLvl, int &lvl) {};

    //! \brief     Communication of the additional physics quantities for parallel purposes
    //! \param     cells                cells vector
    //! \param     dim                  dimension
    virtual void communicationsAddPhys(Cell **cells, const int &dim) { Errors::errorMessage("communicationsAddPhys not implemented for used additional physic"); };
    //! \brief     Communication of the additional physics quantities for parallel purposes with Adaptive Mesh Refinement
    //! \param     cells                cells vector
    //! \param     dim                  dimension
    //! \param     lvl                  level
    virtual void communicationsAddPhysAMR(Cell **cells, const int &dim, const int &lvl) { Errors::errorMessage("communicationsAddPhysAMR not implemented for used additional physic"); };
    
    //! \brief     Return the associated number of the transport equation (only used for surface tension)
    virtual int getNumTransportAssociated() const { Errors::errorMessage("getNumTransportAssociated not implemented for used additional physic"); return 0; };

  protected:
    
  private:
};

#endif // ADDPHYS_H
