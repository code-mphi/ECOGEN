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

#ifndef MODEXTENDEDNAVIERSTOKESKORTEWEG1D_H
#define MODEXTENDEDNAVIERSTOKESKORTEWEG1D_H

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "FluxExtendedNavierStokesKorteweg1D.h"
#include "MixExtendedNavierStokesKorteweg1D.h"

//! \class     ModExtendedNavierStokesKorteweg1D
//! \brief     Model class for Extended Navier-Stokes-Korteweg mathematical system of equations (single phase)
class ModExtendedNavierStokesKorteweg1D : public Model
{
  public:
    //! \brief     ExtendedNavierStokesKorteweg1D model constructor
    //! \param     numbTransports    number of additional transport equations
    ModExtendedNavierStokesKorteweg1D(const int& numbTransports);
    //! \brief     ExtendedNavierStokesKorteweg1D model constructor
    //! \param     numbTransports       number of additional transport equations
    //! \param     lambda               parameter lambda of Navier-Stokes-Korteweg equations
    //! \param     beta                 parameter beta of Navier-Stokes-Korteweg equations
    //! \param     tau                  small parameter controlling dissipartion of Navier-Stokes-Korteweg equations
    //! \param     eps                  small parameter controlling dispersion of Navier-Stokes-Korteweg equations
    ModExtendedNavierStokesKorteweg1D(const int& numbTransports, const double& lambda, const double& beta, const double& tau, const double& eps);
    //! \brief     ExtendedNavierStokesKorteweg1D model constructor
    //! \param     name                 model name
    //! \param     numbTransports       number of additional transport equations
    //! \param     lambda               parameter lambda of Navier-Stokes-Korteweg equations
    //! \param     beta                 parameter beta of Navier-Stokes-Korteweg equations
    //! \param     tau                  small parameter controlling dissipartion of Navier-Stokes-Korteweg equations
    //! \param     eps                  small parameter controlling dispersion of Navier-Stokes-Korteweg equations
    ModExtendedNavierStokesKorteweg1D(
      const std::string& name, const int& numbTransports, const double& lambda, const double& beta, const double& tau, const double& eps);
    ~ModExtendedNavierStokesKorteweg1D() override;

    void allocateCons(Flux** cons) override;
    void allocatePhase(Phase** phase) override;
    void allocateMixture(Mixture** mixture) override;

    //! \details    Does nothing for this model
    void fulfillState(Phase** /*phases*/, Mixture* /*mixture*/) override;

    //! \details    Does nothing for this model
    void fulfillStateResume(Phase** /*phases*/, Mixture* /*mixture*/) override {};

    //! \brief     Complete the augmented variables (such as the ones of Navier-Stokes-Korteweg model)
    //! \param     cell           cell
    void initializeAugmentedVariables(Cell* cell) override;

    //Fluid-flow Riemann solvers
    //--------------------------
    void solveRiemannIntern(Cell& cellLeft,
                            Cell& cellRight,
                            const double& dxLeft,
                            const double& dxRight,
                            double& dtMax,
                            std::vector<double>& /*boundData*/ = DEFAULT_VEC_INTERFACE_DATA) const override;

    void solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& boundData) const override;
    void solveRiemannNullFlux() const override;

    void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const override;

    //Accessors
    //---------
    //! \brief  Select a specific scalar variable
    //! \param  phases         phases array variables
    //! \param  mixture        mixture variables
    //! \param  vecTransports  vector of transport variables
    //! \param  nameVariables  Name of the variable to select
    //! \param  numPhases      Phases number's
    double selectScalar(Phase** phases, Mixture* /*mixture*/, Transport* transports, Variable nameVariable, int num = 0) const override;
    const double& getSM() override
    {
      Errors::errorMessage("getSM not available for required model");
      return Errors::defaultDouble;
    };
    const Coord& getVelocity(const Cell* cell) const override { return cell->getPhase(0)->getVelocity(); };
    Coord& getVelocity(Cell* cell) override { return cell->getPhase(0)->getVelocity(); };

    const std::string& whoAmI() const override { return m_name; };

    //Methods specific to Extended Navier-Stokes-Korteweg
    //----------------------------------
    //! \brief     Compute maximal wave speed using Davis approximation and eigenvalues of hyperbolic equations
    //! \param     cellLeft     cell left
    //! \param     cellRight    cell right
    //! \param     rhoL         density left
    //! \param     rhoR         density right
    //! \param     uL           velocity left
    //! \param     uR           velocity right
    //! \param     etaL         eta left
    //! \param     etaR         eta right
    //! \param     vecPxL       vecPxL left
    //! \param     vecPxR       vecPxR right
    //! \param     kappaL       kappaL left
    //! \param     kappaR       kappaR right
    //! \param     kappaPrimeL  kappaPrimeL left
    //! \param     kappaPrimeR  kappaPrimeR right
    //! \param     kappaSecondL kappaSecondL left
    //! \param     kappaSecondR kappaSecondR right
    //! \param     c2L          squareSoundSpeedL left
    //! \param     c2R          squareSoundSpeedR right
    //! \param     sL           RP left characteristic
    //! \param     sR           RP right characteristic
    virtual void computeMaxWaveSpeed(const double& rhoL,
                                     const double& rhoR,
                                     const double& uL,
                                     const double& uR,
                                     const double& rhoSurEtaL,
                                     const double& rhoSurEtaR,
                                     const double& vecPxL,
                                     const double& vecPxR,
                                     const double& kappaL,
                                     const double& kappaR,
                                     const double& kappaPrimeL,
                                     const double& kappaPrimeR,
                                     const double& kappaSecondL,
                                     const double& kappaSecondR,
                                     const double& c2L,
                                     const double& c2R,
                                     double& sL,
                                     double& sR) const;

  private:
    static const std::string NAME;
};

#endif // MODEXTENDEDNAVIERSTOKESKORTEWEG1D_H
