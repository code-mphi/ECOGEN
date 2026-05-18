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

#ifndef MODSHALLOWWATER_H
#define MODSHALLOWWATER_H

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "FluxShallowWater.h"
#include "MixShallowWater.h"
#include "PhaseShallowWater.h"
#include "GradPhaseShallowWater.h"

//! \class     ModShallowWater
//! \brief     Model class for ShallowWater mathematical system of equations (single phase)
class ModShallowWater : public Model
{
  public:
    //! \brief     ShallowWater model constructor
    //! \param     numbTransports      number of additional transport equations
    ModShallowWater(const int& numbTransports);
    ~ModShallowWater() override;

    void allocateCons(Flux** cons) override;
    void allocatePhase(Phase** phase) override;
    void allocateMixture(Mixture** mixture) override;

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
    void solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& boundData) const override;

    //Transports Riemann solvers
    //--------------------------

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

    const std::string& whoAmI() const override { return m_name; };

  protected:
  private:
    static const std::string NAME;
};

#endif // MODSHALLOWWATER_H
