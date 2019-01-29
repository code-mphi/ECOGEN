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

#ifndef APKVISCOSITY_H
#define APKVISCOSITY_H

//! \file      APKViscosity.h
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "../APKapila.h"
#include "QAPViscosity.h"
#include "../../Eos/Eos.h"

//! \class     APKViscosity
//! \brief     General class for thermal viscosity for the Kapila model
class APKViscosity : public APKapila
{
  public:
    APKViscosity();
    APKViscosity(int& numberQPA, Eos** eos, int &numberPhases, std::string nameFile = "Unknown file");
    virtual ~APKViscosity();

    virtual void addQuantityAddPhys(Cell *cell);

    virtual void solveFluxAddPhys(CellInterface *cellBound, const int &numberPhases);
    virtual void solveFluxAddPhysBoundary(CellInterface *cellBound, const int &numberPhases);
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
    //! \param     distLeft             distance between the center of the left cell and its corresponding edge
    //! \param     distRight            distance between the center of the right cell and its corresponding edge
    //! \param     numberPhases         number of phases
    void solveFluxViscosityInner(Coord &velocityLeft, Coord &velocityRight, Coord &gradULeft, Coord &gradURight, Coord &gradVLeft, Coord &gradVRight, Coord &gradWLeft, Coord &gradWRight, double &muMixLeft, double &muMixRight, double &distLeft, double &distRight, int numberPhases) const;
    //! \brief     Solve the viscosity flux at a boundary with an absorption type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     distLeft             distance between the center of the left cell and its corresponding edge
    //! \param     numberPhases         number of phases
    void solveFluxViscosityAbs(Coord &velocityLeft, Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, double &muMixLeft, double &distLeft, int numberPhases) const;
    //! \brief     Solve the viscosity flux at a boundary with an wall type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     distLeft             distance between the center of the left cell and its corresponding edge
    //! \param     numberPhases         number of phases
    void solveFluxViscosityWall(Coord &velocityLeft, double &muMixLeft, double &distLeft, int numberPhases) const;
    //! \brief     Solve the viscosity flux at a boundary with non-defined type yet
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     distLeft             distance between the center of the left cell and its corresponding edge
    //! \param     numberPhases         number of phases
    void solveFluxViscosityOther(Coord &velocityLeft, Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, double &muMixLeft, double &distLeft, int numberPhases) const;
    virtual void addNonCons(Cell *cell, const int &numberPhases);

    virtual void communicationsAddPhys(Cell **cells, const int &dim);
		virtual void communicationsAddPhysAMR(Cell **cells, const int &dim, const int &lvl);

  protected:
  
  private:
    double *m_muk;            //!< Dynamic viscosity (kg/m/s or Pa.s) of each phase (taken from the EOS classes) (buffer)
    int m_numQPA;             //!< Number of the associated variable for each cell (m_vecGrandeursAddPhys)

    Coord m_velocityLeft;     //!< Left velocity vector for the flux computation (buffer)
    Coord m_gradULeft;        //!< Left gradient of the velocity in the x-direction for the flux computation (buffer)
    Coord m_gradVLeft;        //!< Left gradient of the velocity in the y-direction for the flux computation (buffer)
    Coord m_gradWLeft;        //!< Left gradient of the velocity in the z-direction for the flux computation (buffer)
    Coord m_velocityRight;    //!< Right velocity vector for the flux computation (buffer)
    Coord m_gradURight;       //!< Right gradient of the velocity in the x-direction for the flux computation (buffer)
    Coord m_gradVRight;       //!< Right gradient of the velocity in the y-direction for the flux computation (buffer)
    Coord m_gradWRight;       //!< Right gradient of the velocity in the z-direction for the flux computation (buffer)
    Coord m_normal;           //!< Normal vector of the corresponding face for the flux computation (buffer)
    Coord m_tangent;          //!< Tangent vector of the corresponding face for the flux computation (buffer)
    Coord m_binormal;         //!< Binormal vector of the corresponding face for the flux computation (buffer)
};

#endif // APKVISCOSITY_H
