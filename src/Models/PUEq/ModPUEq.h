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

#ifndef MODPUEQ_H
#define MODPUEQ_H

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "MixPUEq.h"

class ModPUEq;

#include "FluxPUEq.h"

//! \class     ModPUEq
//! \brief     Model class for the pressure-velocity-equilibrium (mechanical equilibrium) system of equations (Kapila)
class ModPUEq : public ModUEq
{
  public:
    //! \brief     PUEq model constructor
    //! \param     numbTransports    number of additional transport equations
    //! \param     numbPhases        number of phases
    //! \param     thresholdPressure threshold for pressure relaxation
    ModPUEq(const int& numbTransports, const int& numbPhases, const double& thresholdPressure);
    ~ModPUEq() override;

    void allocateCons(Flux** cons) override;
    void allocatePhase(Phase** phase) override;
    void allocateMixture(Mixture** mixture) override;

    //! \details    Complete pressures when resuming a simulation
    void fulfillStateResume(Phase** phases, Mixture* mixture) override;

    //! \brief     Energy correction after relaxations (only for PUEq model)
    //! \param     cell              cell
    //! \param     type              primitive-variable type
    void correctionEnergy(Cell* cell, Prim type = vecPhases) const override;

    //Accessors
    //---------
    //! \brief  Select a specific scalar variable
    //! \param  phases         phases array variables
    //! \param  mixture        mixture variables
    //! \param  vecTransports  vector of transport variables
    //! \param  nameVariables  Name of the variable to select
    //! \param  numPhases      Phases number's
    double selectScalar(Phase** phases, Mixture* mixture, Transport* transports, Variable nameVariable, int num = 0) const override;

  private:
    static const std::string NAME;

    friend class FluxPUEq;
};

#endif // MODPUEQ_H
