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

#ifndef MODPTUEQ_H
#define MODPTUEQ_H

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "FluxPTUEq.h"
#include "MixPTUEq.h"

//! \class     ModPTUEq
//! \brief     Model class for pressure-temperature-velocity (mechanical and thermal equilibrium) system of equations
class ModPTUEq : public Model
{
  public:
    //! \brief     Thermal equilibrium model constructor
    //! \param     numbTransports    number of additional transport equations
    //! \param     numbPhases        number of phases
    ModPTUEq(const int& numbTransports, const int& numbPhases);
    virtual ~ModPTUEq();

    virtual void allocateCons(Flux** cons);
    virtual void allocatePhase(Phase** phase);
    virtual void allocateMixture(Mixture** mixture);

    virtual void fulfillState(Phase** phases, Mixture* mixture);

    //! \details    Does nothing for this model
    virtual void fulfillStateRestart(Phase** /*phases*/, Mixture* /*mixture*/) {};

    //! \details    Does nothing for this model
    virtual void initializeAugmentedVariables(Cell* /*cell*/) {};

    //Hydrodynamic Riemann solvers
    //----------------------------
    virtual void solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double> &boundData = DEFAULT_VEC_INTERFACE_DATA) const; // Riemann between two computed cells
    virtual void solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& boundData) const; // Riemann between left cell and wall
    virtual void solveRiemannTank(Cell& cellLeft, const double& dxLeft, double& dtMax, const double* ak0, const double* rhok0, const double& p0, const double& T0, std::vector<double> &boundData) const; // Riemann for tank
    virtual void solveRiemannOutflow(Cell& cellLeft, const double& dxLeft, double& dtMax, const double p0, std::vector<double> &boundData) const; // Riemann for outflow with imposed pressure

    virtual void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const;

    //Accessors
    //---------
    virtual const double& getSM();
    virtual const Coord& getVelocity(const Cell* cell) const { return cell->getMixture()->getVelocity(); };
    virtual Coord& getVelocity(Cell* cell) { return cell->getMixture()->getVelocity(); };

    virtual const std::string& whoAmI() const { return m_name; };

  private:
    static const std::string NAME;
};

#endif // MODPTUEQ_H
