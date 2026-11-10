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

#ifndef FLUXNEQ_H
#define FLUXNEQ_H

#include "../../Flux.h"

class FluxNEq;

#include "../ModNEq/ModNEq.h"

//! \class     FluxNEq
//! \brief     Flux class for the non-equilibrium system of equations
class FluxNEq : public Flux
{
  public:
    FluxNEq(const int& numbPhases);
    ~FluxNEq() override;

    void printFlux() const override;
    void addFlux(double coefA) override;
    void addFlux(Flux* flux) override;
    void subtractFlux(double coefA) override;
    void addFluxRotatingRegion(double coefA) override;
    void subtractFluxRotatingRegion(double coefA) override;
    void multiply(double scalar) override;
    void setBufferFlux(Cell& cell) override;
    void buildCons(Phase** phases, Mixture* /*mixture*/) override;
    void buildPrim(Phase** phases, Mixture* mixture) override;
    virtual void setToZero() override;
    virtual void
    addNonCons(double /*coefA*/, const Cell* /*cell*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override
    {
      Errors::errorMessage("addNonCons not available for required model");
    }
    virtual void
    subtractNonCons(double /*coefA*/, const Cell* /*cell*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override
    {
      Errors::errorMessage("substractNonCons not available for required model");
    }
    //! \brief     Adds a geometric source term to simulate symmetry along one axis (cylindrical axisymmetry) or two axes (spherical axisymmetry).
    //! \param     phases    Phases of the cell.
    //! \param     mixture   Mixture variables of the cell.
    //! \param     r         Radial distance from the symmetry axis.
    //! \param     v         Coefficient determining the type of symmetry.
    virtual void addSymmetricTerms(Phase** phases, Mixture* mixture, const double& r, const double& v) override;

    // Accessors
    //----------
    const double& getAlpha(const int& numPhase) const override { return m_alpha[numPhase]; };
    const double& getMass(const int& numPhase) const override { return m_mass[numPhase]; };
    const double& getEnergy(const int& numPhase) const override { return m_energ[numPhase]; };
    const Coord& getMomentum(const int& numPhase) const override { return m_momentum[numPhase]; };
    void setCons(const Flux* cons) override;

  protected:
    double* m_alpha;   //!< volume fraction array
    double* m_mass;    //!< mass array
    double* m_energ;   //!< specific total energy array
    Coord* m_momentum; //!< momentum array
    double m_uI;       //!< interfacial velocity
    double m_pI;       //!< interfacial pressure

  private:
    friend class ModNEq;
    // To modify if needed, example: To add a class APUEqViscosity, add friend class APUEqViscosity.
    friend class APUEqSurfaceTension;
    friend class APUEqViscosity;
    friend class APUEqConductivity;
    friend class RiemannSolverHLLC;  //!< RiemannSolverHLLC need to access to flux variables to solve discontinuity problem
    friend class RiemannSolverExact; //!< RiemannSolverExact need to access to flux variables to solve discontinuity problem
    friend class ModNEqHLL;
    friend class ModNEqHLLAC;
    friend class ModNEqHLLC4Waves;
};

#endif // FLUXNEQ_H
