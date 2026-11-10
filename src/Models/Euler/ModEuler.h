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

#ifndef MODEULER_H
#define MODEULER_H

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "FluxEuler.h"
#include "MixEuler.h"
#include "PhaseEuler.h"
#include "GradPhaseEuler.h"
#include "GradMixEuler.h"

//! \class     ModEuler
//! \brief     Model class for Euler mathematical system of equations (single phase)
class ModEuler : public Model
{
  public:
    //! \brief     Euler model constructor
    //! \param     numbTransports      number of additional transport equations
    ModEuler(const int& numbTransports);
    ~ModEuler() override;

    void allocateCons(Flux** cons) override;
    void allocatePhase(Phase** phase) override;
    void allocateMixture(Mixture** mixture) override;
    void allocatePhaseGradient(GradPhase** phase) override;
    void allocateMixtureGradient(GradMixture** mixture) override;

    //! \details    Complete single fluid state from pressure, density and velocity
    void fulfillState(Phase** phases, Mixture* /*mixture*/) override;

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
                            std::vector<double>& boundData = DEFAULT_VEC_INTERFACE_DATA) const override;
    void solveRiemannInternMRF(Cell& cellLeft,
                               Cell& cellRight,
                               const double& dxLeft,
                               const double& dxRight,
                               double& dtMax,
                               const Coord& omega,
                               const Coord& normal,
                               const Coord& tangent,
                               const Coord& binormal,
                               const Coord& position) const override;
    void solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& boundData) const override;
    void
    solveRiemannMovingWall(Cell& cellLeft, const double& dxLeft, double& dtMax, double wallVelocity, std::vector<double>& boundData) const override;
    void solveRiemannInletTank(Cell& cellLeft,
                               const double& dxLeft,
                               double& dtMax,
                               const double* /*ak0*/,
                               const double* rhok0,
                               const double& p0,
                               const double& /*T0*/,
                               std::vector<double>& boundData) const override;
    void solveRiemannInletInjStagState(Cell& cellLeft,
                                       const double& dxLeft,
                                       double& dtMax,
                                       const double m0,
                                       const double* /*ak0*/,
                                       const double* rhok0,
                                       const double* pk0,
                                       std::vector<double>& boundData) const override;
    void solveRiemannInletInjTemp(Cell& cellLeft,
                                  const double& dxLeft,
                                  double& dtMax,
                                  const double m0,
                                  const double* Tk0,
                                  const double* /*ak0*/,
                                  std::vector<double>& boundData) const override;
    void
    solveRiemannOutletPressure(Cell& cellLeft, const double& dxLeft, double& dtMax, const double p0, std::vector<double>& boundData) const override;
    void
    solveRiemannOutletMassflow(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, std::vector<double>& boundData) const override;
    void solveRiemannNullFlux() const override;

    //Transports Riemann solvers
    //--------------------------
    void solveRiemannTransportIntern(Cell& cellLeft, Cell& cellRight) override;
    void solveRiemannTransportWall() override;
    void solveRiemannTransportMovingWall(Cell& cellLeft, double wallVelocity) override;
    void solveRiemannTransportInletTank(Cell& cellLeft, double* valueTransports) override;
    void solveRiemannTransportInletInjStagState(Cell& cellLeft, double* valueTransports) override;
    void solveRiemannTransportOutletPressure(Cell& cellLeft, double* valueTransport) override;

    void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const override;

    //low-Mach preconditioning
    //------------------------
    void lowMachSoundSpeed(
      double& machRef, const double& uL, double& cL, const double& uR = Errors::defaultDouble, double& cR = Tools::uselessDouble) const override;

    //Moving Reference Frame
    //----------------------
    void addNonConsMrfFlux(Phase** phases) override;
    void reverseProjectionMrfFlux(const Coord normal, const Coord tangent, const Coord binormal) const override;

    //Accessors
    //---------
    //! \brief  Select a specific scalar variable
    //! \param  phases         phases array variables
    //! \param  mixture        mixture variables
    //! \param  vecTransports  vector of transport variables
    //! \param  nameVariables  Name of the variable to select
    //! \param  numPhases      Phases number's
    double selectScalar(Phase** phases, Mixture* /*mixture*/, Transport* transports, Variable nameVariable, int num = 0) const override;
    const double& getSM() override;
    const Coord& getVelocity(const Cell* cell) const override { return cell->getPhase(0)->getVelocity(); };
    Coord& getVelocity(Cell* cell) override { return cell->getPhase(0)->getVelocity(); };

    const std::string& whoAmI() const override { return m_name; };

    void setSmoothCrossSection1d(const bool& applySmooth) override { m_smoothCrossSection1d = applySmooth; };

  protected:
  private:
    static const std::string NAME;
};

#endif // MODEULER_H
