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

#ifndef APKSURFACETENSION_H
#define APKSURFACETENSION_H

//! \file      APKSurfaceTension.h
//! \author    K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "../APKapila.h"
#include "QAPSurfaceTension.h"

//! \class     APKSurfaceTension
//! \brief     General class for surface tension for the Kapila model
class APKSurfaceTension : public APKapila
{
  public:
    APKSurfaceTension();
    APKSurfaceTension(tinyxml2::XMLElement *element, int& numberQPA, std::vector<std::string> nameTransports, std::vector<std::string> namePhases, std::string nameFichier = "Unknown file");
    virtual ~APKSurfaceTension();

    virtual void addQuantityAddPhys(Cell *cell);

    virtual double computeEnergyAddPhys(QuantitiesAddPhys* QPA);
    virtual void solveFluxAddPhys(CellInterface *cellInterface, const int &numberPhases);
    virtual void solveFluxAddPhysBoundary(CellInterface *cellInterface, const int &numberPhases);
    //! \brief     Solve the surface-tension flux between two cells
    //! \param     velocityLeft         velocity of the left cell
    //! \param     velocityRight        velocity of the right cell
    //! \param     gradCLeft            color function gradient of the left cell
    //! \param     gradCRight           color function gradient of the right cell
    void solveFluxSurfaceTensionInner(Coord &velocityLeft, Coord &velocityRight, Coord &gradCLeft, Coord &gradCRight) const;
    //! \brief     Solve the surface-tension flux at a boundary with an non-reflecting type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradCLeft            color function gradient of the left cell
    void solveFluxSurfaceTensionNonReflecting(Coord &velocityLeft, Coord &gradCLeft) const;
    //! \brief     Solve the surface-tension flux at a boundary with a wall type
    //! \param     gradCLeft            color function gradient of the left cell
    void solveFluxSurfaceTensionWall(Coord &gradCLeft) const;
    //! \brief     Solve the surface-tension flux at a boundary with an outflow type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradCLeft            color function gradient of the left cell
    void solveFluxSurfaceTensionOutflow(Coord &velocityLeft, Coord &gradCLeft) const;
    //! \brief     Solve the surface-tension flux at a boundary with an inflow (injection) type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradCLeft            color function gradient of the left cell
    void solveFluxSurfaceTensionInflow(Coord &velocityLeft, Coord &gradCLeft) const;
    //! \brief     Solve the surface-tension flux at a boundary with non-defined type yet
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradCLeft            color function gradient of the left cell
    void solveFluxSurfaceTensionOther(Coord &velocityLeft, Coord &gradCLeft) const;
    virtual void addNonCons(Cell *cell, const int &numberPhases) {}; //The surface-tension effects do not involve non-conservative terms.
    virtual void addSymmetricTermsRadialAxisOnX(Cell *cell, const int &numberPhases);
    virtual void addSymmetricTermsRadialAxisOnY(Cell *cell, const int &numberPhases);

    virtual void reinitializeColorFunction(std::vector<Cell *> *cellsLvl, int &lvl);
    virtual bool reinitializationActivated() { return m_reinitializationActivated; };

    virtual void communicationsAddPhys(int numberPhases, const int &dim, const int &lvl);
    virtual const int& getNumTransportAssociated() const { return m_numTransportAssociated; };

  protected:
  
  private:
    std::string m_nameTransportAssociated; //!< Name of the associated variable for each cell (m_vecTransports)
    int m_numTransportAssociated;          //!< Number of the associated variable for each cell (m_vecTransports)
    int m_numQPAGradC;                     //!< Number of the associated variable for each cell (m_vecGrandeursAddPhys)
    double m_sigma;                        //!< Surface tension coefficient
    bool m_reinitializationActivated;      //!< Reinitialization of the transported variable with the volume fraction of a phase
    std::string m_namePhaseAssociated;     //!< Name of the associated variable for each cell (m_vecPhases)
    int m_numPhaseAssociated;              //!< Number of the associated variable for each cell (m_vecPhases)

    Coord m_velocityLeft;                  //!< Left velocity vector for the flux computation (buffer)
    Coord m_gradCLeft;                     //!< Left gradient of the color function for the flux computation (buffer)
    Coord m_velocityRight;                 //!< Right velocity vector for the flux computation (buffer)
    Coord m_gradCRight;                    //!< Right gradient of the color function vector for the flux computation (buffer)
    Coord m_normal;                        //!< Normal vector of the corresponding face for the flux computation (buffer)
    Coord m_tangent;                       //!< Tangent vector of the corresponding face for the flux computation (buffer)
    Coord m_binormal;                      //!< Binormal vector of the corresponding face for the flux computation (buffer)
};

#endif // APKSURFACETENSION_H
