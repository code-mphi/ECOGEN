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

#ifndef MODNEQ_H
#define MODNEQ_H

#include "../../Model.h"
#include "../../../Order1/Cell.h"
#include "../MixNEq.h"
#include "../PhaseNEq.h"
#include "../GradPhaseNEq.h"
#include "../GradMixNEq.h"

class ModNEq;

#include "../FluxNEq/FluxNEq.h"

//! \class     ModNEq
//! \brief     Model class for the non-equilibrium system of equations
class ModNEq : public Model
{
  public:
    //! \brief     NEq model constructor
    //! \param     numbTransports    number of additional transport equations
    //! \param     numbPhases        number of phases
    //! \param     solver            name of riemann problem solver
    ModNEq(const int& numbTransports, const int& numbPhases);
    //! \brief     Generic model constructor (used by derived classes)
    //! \param     name              model name
    //! \param     numbTransports    number of additional transport equations
    //! \param     solver            name of riemann problem solver
    ModNEq(const std::string& name, const int& numbTransports);
    ~ModNEq() override;

    void allocateCons(Flux** cons) override;
    void allocatePhase(Phase** phase) override;
    void allocateMixture(Mixture** mixture) override;
    void allocatePhaseGradient(GradPhase** phase) override;
    void allocateMixtureGradient(GradMixture** mixture) override;

    //! \details    Complete multiphase state from volume fractions, pressures, densities and velocity
    void fulfillState(Phase** phases, Mixture* mixture) override;

    //! \details    Does nothing for this model
    void fulfillStateResume(Phase** /*phases*/, Mixture* /*mixture*/) override {};

    //! \details    Does nothing for this model
    void initializeAugmentedVariables(Cell* /*cell*/) override {};

    //Fluid-flow Riemann solvers
    //--------------------------
    void solveRiemannIntern(Cell& cellLeft,
                            Cell& cellRight,
                            const double& dxLeft,
                            const double& dxRight,
                            double& dtMax,
                            std::vector<double>& /*boundData*/ = DEFAULT_VEC_INTERFACE_DATA) const override; // Riemann between two computed cells

    //Riemann solvers from litterature
    //--------------------------------
    //! \brief   Call appropriate solver
    virtual void solver(Cell& /*cellLeft*/,
                        Cell& /*cellRight*/,
                        const double& /*dxLeft*/,
                        const double& /*dxRight*/,
                        double& /*dtMax*/,
                        std::vector<double>& /*boundData*/ = DEFAULT_VEC_INTERFACE_DATA) const
    {
      Errors::errorMessage("RiemannSolver not available");
    };

    //Transports Riemann solvers
    //--------------------------
    // virtual void solveRiemannTransportIntern(Cell& cellLeft, Cell& cellRight);
    // virtual void solveRiemannTransportWall();
    // virtual void solveRiemannTransportMovingWall(Cell& cellLeft, double wallVelocity);
    // virtual void solveRiemannTransportInletTank(Cell& cellLeft, double* valueTransports);
    // virtual void solveRiemannTransportInletInjStagState(Cell& cellLeft, double* valueTransports);
    // virtual void solveRiemannTransportOutletPressure(Cell& cellLeft, double* valueTransport);

    void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const override;

    //Low-Mach preconditioning
    //------------------------
    // virtual void lowMachSoundSpeed(double& machRef, const double& uL, double& cL, const double& uR = Errors::defaultDouble, double& cR = Tools::uselessDouble) const;

    //Moving Reference Frame
    //----------------------
    // virtual void addNonConsMrfFlux(Phase** phases);
    // virtual void reverseProjectionMrfFlux(const Coord normal, const Coord tangent, const Coord binormal) const;

    //Accessors
    //---------
    //! \brief  Select a specific scalar variable
    //! \param  phases         phases array variables
    //! \param  mixture        mixture variables
    //! \param  vecTransports  vector of transport variables
    //! \param  nameVariables  Name of the variable to select
    //! \param  numPhases      Phases number's
    double selectScalar(Phase** phases, Mixture* mixture, Transport* transports, Variable nameVariable, int num = 0) const override;
    const double& getSM() override;
    const Coord& getVelocity(const Cell* cell) const override { return cell->getMixture()->getVelocity(); };
    Coord& getVelocity(Cell* cell) override { return cell->getMixture()->getVelocity(); };

    const std::string& whoAmI() const override { return m_name; };
    // virtual void setSmoothCrossSection1d(const bool& applySmooth) { m_smoothCrossSection1d = applySmooth; };

    friend class FluxNEq;

  private:
    static const std::string NAME;
};

#endif // MODNEQ_H
