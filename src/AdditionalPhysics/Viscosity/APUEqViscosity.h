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

#ifndef APUEQVISCOSITY_H
#define APUEQVISCOSITY_H

#include "../APUEq.h"
#include "QAPViscosity.h"
#include "../../Eos/Eos.h"

//! \class     APUEqViscosity
//! \brief     General class for thermal viscosity for the velocity-equilibrium system of equations
class APUEqViscosity : public APUEq
{
  public:
    APUEqViscosity(int& numberQPA, Eos** eos, const int& numberPhases);
    virtual ~APUEqViscosity();

    virtual void addQuantityAddPhys(Cell* cell);

    virtual void solveFluxAddPhys(CellInterface* cellInterface, const int& numberPhases);
    virtual void solveFluxAddPhysBoundary(CellInterface* cellInterface, const int& numberPhases);
    //! \brief     Solve the viscosity flux between two cells
    //! \param     velocityLeft         velocity of the left cell
    //! \param     velocityRight        velocity of the right cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradURight           gradient of the velocity in the x-direction of the right cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradVRight           gradient of the velocity in the y-direction of the right cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     gradWRight           gradient of the velocity in the z-direction of the right cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     muMixRight           dynamic viscosity of the mixture of the right cell
    //! \param     numberPhases         number of phases
    void solveFluxViscosityInner(const Coord& velocityLeft, const Coord& velocityRight, const Coord& gradULeft, const Coord& gradURight, const Coord& gradVLeft, const Coord& gradVRight,
      const Coord& gradWLeft, const Coord& gradWRight, const double& muMixLeft, const double& muMixRight, const int& numberPhases) const;
    //! \brief     Solve the viscosity flux at a boundary with a non-reflecting type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     numberPhases         number of phases
    void solveFluxViscosityNonReflecting(const Coord& velocityLeft, const Coord& gradULeft, const Coord& gradVLeft, const Coord& gradWLeft, const double& muMixLeft, const int& numberPhases) const;
    //! \brief     Solve the viscosity flux at a boundary with a wall type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     distLeft             distance between the center of the left cell and its corresponding edge
    //! \param     numberPhases         number of phases
    void solveFluxViscosityWall(const Coord& velocityLeft, const double& muMixLeft, const double& distLeft, const int& numberPhases) const;
    //! \brief     Solve the viscosity flux at a boundary with a symmetry type
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     numberPhases         number of phases
    void solveFluxViscositySymmetry(const Coord& gradULeft, const Coord& gradVLeft, const Coord& gradWLeft, const double& muMixLeft, const int& numberPhases) const;
    //! \brief     Solve the viscosity flux at a boundary with non-defined type yet
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     numberPhases         number of phases
    void solveFluxViscosityOther(const int& numberPhases) const;
    virtual void addNonCons(Cell* cell, const int& numberPhases);
    virtual void addSymmetricTermsRadialAxisOnX(Cell* cell, const int& numberPhases);
    virtual void addSymmetricTermsRadialAxisOnY(Cell* cell, const int& numberPhases);

    virtual void communicationsAddPhys(const int& /*numberPhases*/, const int& dim, const int& lvl);

  protected:
  
  private:
    double* m_muk;            //!< Dynamic viscosity (kg/m/s or Pa.s) of each phase (taken from the EOS classes) (buffer)
    int m_numQPA;             //!< Number of the associated variable for each cell (m_vecGrandeursAddPhys)

    Coord m_velocityLeft;     //!< Left velocity vector for the flux computation (buffer)
    Coord m_gradULeft;        //!< Left gradient of the velocity in the x-direction for the flux computation (buffer)
    Coord m_gradVLeft;        //!< Left gradient of the velocity in the y-direction for the flux computation (buffer)
    Coord m_gradWLeft;        //!< Left gradient of the velocity in the z-direction for the flux computation (buffer)
    Tensor m_tensorLeft;      //!< Left tensor of the velocity gradients (buffer)
    Coord m_velocityRight;    //!< Right velocity vector for the flux computation (buffer)
    Coord m_gradURight;       //!< Right gradient of the velocity in the x-direction for the flux computation (buffer)
    Coord m_gradVRight;       //!< Right gradient of the velocity in the y-direction for the flux computation (buffer)
    Coord m_gradWRight;       //!< Right gradient of the velocity in the z-direction for the flux computation (buffer)
    Tensor m_tensorRight;     //!< Right tensor of the velocity gradients (buffer)
    Coord m_normal;           //!< Normal vector of the corresponding face for the flux computation (buffer)
    Coord m_tangent;          //!< Tangent vector of the corresponding face for the flux computation (buffer)
    Coord m_binormal;         //!< Binormal vector of the corresponding face for the flux computation (buffer)
};

#endif // APUEQVISCOSITY_H
