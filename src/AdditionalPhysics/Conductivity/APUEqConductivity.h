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

#ifndef APUEQCONDUCTIVITY_H
#define APUEQCONDUCTIVITY_H

#include "../APUEq.h"
#include "QAPConductivity.h"
#include "../../Eos/Eos.h"

//! \class     APUEqConductivity
//! \brief     General class for thermal conductivity for the velocity-equilibrium system of equations
class APUEqConductivity : public APUEq
{
  public:
    APUEqConductivity();
    APUEqConductivity(int& numberQPA, Eos** eos, const int& numberPhases);
    virtual ~APUEqConductivity();

    virtual void addQuantityAddPhys(Cell* cell);

    virtual void solveFluxAddPhys(CellInterface* cellInterface, const int& numberPhases);
    virtual void solveFluxAddPhysBoundary(CellInterface* cellInterface, const int& numberPhases);
    //! \brief     Solve the conductivity flux between two cells
    //! \param     gradTkLeft           temperature gradient of phase k of the left cell
    //! \param     gradTkRight          temperature gradient of phase k of the right cell
    //! \param     alphakL              volume fraction of phase k of the left cell
    //! \param     alphakR              volume fraction of phase k of the right cell
    //! \param     numPhase             number of the phase
    void solveFluxConductivityInner(const Coord& gradTkLeft, const Coord& gradTkRight, const double& alphakL, const double& alphakR, const int& numPhase) const;
    //! \brief     Solve the conductivity flux at a boundary with an non-reflecting type
    //! \param     gradTkLeft           temperature gradient of phase k of the left cell
    //! \param     alphakL              volume fraction of phase k of the left cell
    //! \param     numPhase             number of the phase
    void solveFluxConductivityNonReflecting(const Coord& gradTkLeft, const double& alphakL, const int& numPhase) const;
    //! \brief     Solve the conductivity flux at a boundary with an wall type
    void solveFluxConductivityWall() const;
    //! \brief     Solve the conductivity flux at a boundary with an outflow type
    void solveFluxConductivityOutflow() const;
    //! \brief     Solve the conductivity flux at a boundary with an inflow type
    void solveFluxConductivityInflow() const;
    //! \brief     Solve the conductivity flux at a boundary with non-defined type yet
    void solveFluxConductivityOther() const;
    virtual void addNonCons(Cell* /*cell*/, const int& /*numberPhases*/) {}; //The conductivity does not involve non-conservative terms.

    virtual void communicationsAddPhys(const int& numberPhases, const int& dim, const int& lvl);

  protected:

  private:
    double* m_lambdak;        //!< Thermal conductivity (W/(m.K)) of each phase (taken from the EOS classes) (buffer)
    int m_numQPA;             //!< Number of the associated variable for each cell (m_vecGrandeursAddPhys)

    Coord m_gradTkLeft;       //!< Left gradient of the corresponding phase temperature for the flux computation (buffer)
    Coord m_gradTkRight;      //!< Right gradient of the corresponding phase temperature for the flux computation (buffer)
    Coord m_normal;           //!< Normal vector of the corresponding face for the flux computation (buffer)
    Coord m_tangent;          //!< Tangent vector of the corresponding face for the flux computation (buffer)
    Coord m_binormal;         //!< Binormal vector of the corresponding face for the flux computation (buffer)
};

#endif // APUEQCONDUCTIVITY_H
