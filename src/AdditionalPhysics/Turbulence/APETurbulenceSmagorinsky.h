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

#ifndef APETURBULENCESMAGORINSKY_H
#define APETURBULENCESMAGORINSKY_H

#include "../APEuler.h"
#include "QAPTurbulence.h"
#include "../../Eos/Eos.h"

//! \class     APUEqTurbulence
//! \brief     General class for handling turbulence for the velocity-equilibrium system of equations
class APETurbulenceSmagorinsky : public APEuler
{
  public:
    APETurbulenceSmagorinsky();
    //! \brief    APETurbulenceSmagorinsky constructor
    //! \param    numberQPA             corresponding additional physic quantities number
    //! \param    eos                   equation of state
    //! \param    elementSmagorinsky    XML element
    //! \param    fileName              string name of readed XML file
    APETurbulenceSmagorinsky(int& numberQPA, Eos** eos, tinyxml2::XMLElement* elementSmagorinsky, std::string fileName);
    ~APETurbulenceSmagorinsky() override;

    void addQuantityAddPhys(Cell* cell) override;

    void solveFluxAddPhys(CellInterface* cellInterface) override;
    void solveFluxAddPhysBoundary(CellInterface* cellInterface) override;
    //! \brief     Solve the turbulence flux between two cells
    //! \param     velocityLeft         velocity of the left cell
    //! \param     velocityRight        velocity of the right cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradURight           gradient of the velocity in the x-direction of the right cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradVRight           gradient of the velocity in the y-direction of the right cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     gradWRight           gradient of the velocity in the z-direction of the right cell
    //! \param     gradTLeft            gradient of the temperature of the left cell
    //! \param     gradTRight           gradient of the temperature of the right cell
    //! \param     gamma                polytropic constant of the fluid
    //! \param     cv                   volume calorific energy of the fluid
    //! \param     mu                   dynamic molecular viscosity of the fluid
    //! \param     rhoLeft              density of the left cell
    //! \param     rhoRight             density of the right cell
    //! \param     cellVolumeLeft       volume of the left cell
    //! \param     strainRateLeft       strain rate tensor of the left cell
    //! \param     strainRateRight      strain rate tensor of the right cell
    void solveFluxTurbulenceInner(CellInterface* cellInterface,
                                  const Coord& velocityLeft,
                                  const Coord& velocityRight,
                                  const Coord& gradULeft,
                                  const Coord& gradURight,
                                  const Coord& gradVLeft,
                                  const Coord& gradVRight,
                                  const Coord& gradWLeft,
                                  const Coord& gradWRight,
                                  const Coord& gradTLeft,
                                  const Coord& gradTRight,
                                  const double& gamma,
                                  const double& cv,
                                  const double& mu,
                                  const double& lambda,
                                  const double& rhoLeft,
                                  const double& rhoRight,
                                  const double& cellVolumeLeft,
                                  const Tensor& strainRateLeft,
                                  const Tensor& strainRateRight) const;
    //! \brief     Solve the turbulence flux at a boundary with a non-reflecting type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     gradTLeft            gradient of the temperature of the left cell
    //! \param     rhoLeft              density of the left cell
    //! \param     cellVolumeLeft       volume of the left cell
    //! \param     strainRateLeft       strain rate tensor of the left cell
    void solveFluxTurbulenceNonReflecting(CellInterface* cellInterface,
                                          const Coord& velocityLeft,
                                          const Coord& gradULeft,
                                          const Coord& gradVLeft,
                                          const Coord& gradWLeft,
                                          const Coord& gradTLeft,
                                          const double& rhoLeft,
                                          const double& cellVolumeLeft,
                                          const Tensor& strainRateLeft) const;
    //! \brief     Solve the turbulence flux at a boundary with a wall type in MRF case
    //! \param     velocityLeft         velocity of the left cell
    //! \param     densityLeft          density of the left cell
    //! \param     distLeft             distance between the center of the left cell and its corresponding edge
    //! \param     velocityWall         velocity of the wall in the Moving Reference Frame context (MRF)
    // void solveFluxTurbulenceWall(const Coord& velocityLeft,
    //                             const double& rhoLeft,
    //                             const double& distLeft,
    //                             Coord& velocityWall = Coord::defaultCoordNonConst) const;
    //! \brief     Solve the turbulence flux at a boundary with a symmetry type
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     rhoLeft              density of the left cell
    //! \param     cellVolumeLeft       volume of the left cell
    //! \param     strainRateLeft       strain rate tensor of the left cell
    //! \param     mu                   dynamic molecular viscosity of the fluid
    void solveFluxTurbulenceSymmetry(CellInterface* cellInterface,
                                     const Coord& gradULeft,
                                     const Coord& gradVLeft,
                                     const Coord& gradWLeft,
                                     const double& rhoLeft,
                                     const double& cellVolumeLeft,
                                     const Tensor& strainRateLeft,
                                     const double& mu) const;
    //! \brief     Solve the turbulence flux at a boundary with non-defined type yet
    void solveFluxTurbulenceOther() const;
    void addNonCons(Cell* /*cell*/) override {};
    void addSymmetricTermsRadialAxisOnX(Cell* /*cell*/) override
    {
      Errors::errorMessage("addSymmetricTermsX not implemented for used additional physic");
    };
    void addSymmetricTermsRadialAxisOnY(Cell* /*cell*/) override
    {
      Errors::errorMessage("addSymmetricTermsY not implemented for used additional physic");
    };

    void communicationsAddPhys(const int& dim, const int& lvl) override;

    //! \brief     Print the diagnostics for turbulence decay cases
    virtual void computeTurbulenceParameters() const override;

  protected:
  private:
    bool m_Ducros;   //!< Bool to turn on/off the Ducros flow sensor
    double m_Cs;     //!< Coefficient modelling the off-diagonal terms of the SGS stress tensor (nondimensional) (buffer)
    double m_Ci;     //!< Coefficient modelling the trace of the SGS stress tensor (nondimensional) (buffer)
    double m_PrT;    //!< Turbulent Prandtl number (nondimensional) (buffer)
    double m_cv;     //!< Volume calorific energy (J/kg/K) (taken from the EOS class) (buffer)
    double m_gamma;  //!< Polytropic constant (nondimensional) (taken from the EOS class) (buffer)
    double m_mu;     //!< Dynamic viscosity (Pa.s) (molecular: taken from the EOS class) (buffer)
    double m_lambda; //!< Thermal conductivity (W/(m.K)) of phase (taken from the EOS class) (buffer)
    int m_numQPA;    //!< Number of the associated variable for each cell (m_vecGrandeursAddPhys)

    Coord m_velocityLeft;  //!< Left velocity vector for the flux computation (buffer)
    Coord m_gradULeft;     //!< Left gradient of the velocity in the x-direction for the flux computation (buffer)
    Coord m_gradVLeft;     //!< Left gradient of the velocity in the y-direction for the flux computation (buffer)
    Coord m_gradWLeft;     //!< Left gradient of the velocity in the z-direction for the flux computation (buffer)
    Tensor m_tensorLeft;   //!< Left tensor of the velocity gradients (buffer)
    Coord m_gradTLeft;     //!< Left vector of the temperature gradients (buffer)
    Coord m_velocityRight; //!< Right velocity vector for the flux computation (buffer)
    Coord m_gradURight;    //!< Right gradient of the velocity in the x-direction for the flux computation (buffer)
    Coord m_gradVRight;    //!< Right gradient of the velocity in the y-direction for the flux computation (buffer)
    Coord m_gradWRight;    //!< Right gradient of the velocity in the z-direction for the flux computation (buffer)
    Tensor m_tensorRight;  //!< Right tensor of the velocity gradients (buffer)
    Coord m_gradTRight;    //!< Right vector of the temperature gradients (buffer)
    Coord m_normal;        //!< Normal vector of the corresponding face for the flux computation (buffer)
    Coord m_tangent;       //!< Tangent vector of the corresponding face for the flux computation (buffer)
    Coord m_binormal;      //!< Binormal vector of the corresponding face for the flux computation (buffer)
};

#endif // APETURBULENCESMAGORINSKY_H
