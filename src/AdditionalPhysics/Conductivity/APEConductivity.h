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
//  If not, see <http://www.gnu.org/licenses/>.#pragma once

#ifndef APECONDUCTIVITY_H
#define APECONDUCTIVITY_H

#include "../APEuler.h"
#include "QAPConductivity.h"
#include "../../Eos/Eos.h"

//! \class APEConductivity
//! \brief General class for thermal conductivity for the Euler model
class APEConductivity : public APEuler
{
public:
  APEConductivity(int& numberQPA, Eos** eos);
  virtual ~APEConductivity();
  
  virtual void addQuantityAddPhys(Cell* cell);

  virtual void solveFluxAddPhys(CellInterface* cellInterface, const int& numberPhases);
  
  virtual void solveFluxAddPhysBoundary(CellInterface* cellInterface, const int& numberPhases);
  //! \brief     Solve the conductivity flux between two cells
  //! \param     gradTLeft           temperature gradient of phase k of the left cell
  //! \param     gradTRight          temperature gradient of phase k of the right cell
  void solveFluxConductivityInner(const Coord& gradTLeft, const Coord& gradTRight) const;
  //! \brief     Solve the conductivity flux at a boundary with an non-reflecting type
  //! \param     gradTLeft           temperature gradient of phase k of the left cell
  void solveFluxConductivityNonReflecting(const Coord& gradTLeft) const;
  //! \brief     Solve the conductivity flux at a wall boundary with imposed temperature
  //! \param     cellInterface       allows to retrieve imposed temp. and distance cell/boundary
  void solveFluxConductivityWallImposedTemp(CellInterface* cellInterface);
  //! \brief     Solve the conductivity flux at a wall boundary with imposed flux density
  //! \param     cellInterface       allows to retrieve imposed flux density
  void solveFluxConductivityWallImposedFlux(CellInterface* cellInterface);
  //! \brief     Solve the conductivity flux at a boundary with non-defined type yet or for adiabatic wall
  void solveFluxConductivityOther() const;

  virtual void addNonCons(Cell* /*cell*/, const int& /*numberPhases*/) {}; //The conductivity does not involve non-conservative terms.
  virtual void communicationsAddPhys(const int& /*numberPhases*/, const int& dim, const int& lvl);

private:
  double m_lambda;          //!< Thermal conductivity (W/(m.K)) of phase (taken from the EOS classe) (buffer)
  int m_numQPA;             //!< Number of the associated variable for each cell (m_vecGrandeursAddPhys)

  Coord m_gradTLeft;        //!< Left gradient of the phase temperature for the flux computation (buffer)
  Coord m_gradTRight;       //!< Right gradient of the phase temperature for the flux computation (buffer)
  Coord m_normal;           //!< Normal vector of the face for the flux computation (buffer)
  Coord m_tangent;          //!< Tangent vector of the face for the flux computation (buffer)
  Coord m_binormal;         //!< Binormal vector of the face for the flux computation (buffer)
};

#endif // APECONDUCTIVITY_H