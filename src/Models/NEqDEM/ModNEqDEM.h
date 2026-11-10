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

#ifndef MODNEQDEM_H
#define MODNEQDEM_H

#include "../NEq/ModNEq/ModNEq.h"
#include "RiemannSolvers/HeaderRiemannSolvers.h"

class ModNEqDEM;

#include "FluxNEqDEM.h"
#include "../Euler/FluxEuler.h"

//! \class     ModNEqDEM
//! \brief     Model class for the non-equilibrium system of equations using the DEM
class ModNEqDEM : public ModNEq
{
  public:
    //! \brief     NEq model constructor
    //! \param     numbTransports    number of additional transport equations
    //! \param     numbPhases        number of phases
    //! \param     solverType        bool: false = EXACT, true = HLLC
    //! \param     typeContactAreas  type of contact areas. int: 1 = symmetric, 2 = nonsymmetric
    ModNEqDEM(const int& numbTransports, const int& numbPhases, const bool& solverType, TypeContactAreas typeContactAreas);
    //! \brief     Generic model constructor (used by derived classes)
    //! \param     name              model name
    //! \param     numbTransports    number of additional transport equations
    ModNEqDEM(const std::string& name, const int& numbTransports);
    ~ModNEqDEM() override;

    void allocateCons(Flux** cons) override;

    //Fluid-flow Riemann solvers
    //--------------------------
    void solveRiemannIntern(Cell& cellLeft,
                            Cell& cellRight,
                            const double& dxLeft,
                            const double& dxRight,
                            double& dtMax,
                            std::vector<double>& /*boundData*/ = DEFAULT_VEC_INTERFACE_DATA) const override; // Riemann between two computed cells
    //! \brief  Call Riemann solver with wall BC
    void solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& /*boundData*/) const override;
    // virtual void solveRiemannInternMRF(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, const Coord& omega,
    //                                    const Coord& normal, const Coord& tangent, const Coord& binormal, const Coord& position) const;
    // virtual void solveRiemannInletTank(Cell& cellLeft, const double& dxLeft, double& dtMax, const double* ak0, const double* rhok0, const double& p0,
    //                                    const double& /*T0*/, std::vector<double>& boundData) const;

    // virtual void solveRiemannInletInjStagState(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, const double* ak0,
    //                                            const double* rhok0, const double* pk0, std::vector<double>& boundData) const;

    // virtual void solveRiemannInletInjTemp(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, const double* Tk0, const double* ak0,
    //                                       std::vector<double>& boundData) const;

    // virtual void solveRiemannOutletPressure(Cell& cellLeft, const double& dxLeft, double& dtMax, const double p0, std::vector<double>& boundData) const;

    // virtual void solveRiemannOutletMassflow(Cell& cellLeft, const double& dxLeft, double& dtMax, const double m0, std::vector<double>& boundData) const;

    // virtual void solveRiemannNullFlux() const;

    void setPhaseSlope(Phase& phaseSlope, int k) override;
    void setSlopeDistance(const double& slopeDistance) override;

  protected:
    RiemannSolver* m_RiemannSolver; //!< Riemann solver pointer to either specific RP Euler solvers such as the HLLC or the exact (for Euler model).
    mutable FluxEuler m_fluxEuler;  //!< Buffer to temporary store each phase flux in EulerHLLC and EulerExact methods.
                                    //!< Declared as a member to not pay the price of the instanciation at each call of the Riemann solver.
    mutable std::vector<PhaseNEq>
      m_phaseSlope; //!< Buffer to store the in-cell phase (primitive) slopes of each phase for second-order extrapolation in Riemann solvers, to compute in-cell Lagrangian term
    mutable double m_slopeDistance = 0.;      //!< Geometric distance over which the in-cell slopes extrapolate (cell width along the face normal)
    mutable std::vector<double> m_alphaLeft;  //!< Volume fractions from the left cell (used to compute contact areas)
    mutable std::vector<double> m_alphaRight; //!< Volume fractions from the right cell (used to compute contact areas)

  private:
    static const std::string NAME;
};

#endif // MODNEQDEM_H
