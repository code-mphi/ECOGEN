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

#ifndef FLUXEULER_H
#define FLUXEULER_H

#include "../Flux.h"

//! \class     FluxEuler
//! \brief     Model class for Euler Flux (single phase)
class FluxEuler : public Flux
{
  public:
    FluxEuler();
    ~FluxEuler() override;

    void printFlux() const override;
    void addFlux(double coefA) override;
    void addFlux(Flux* flux) override;
    void subtractFlux(double coefA) override;
    void addFluxRotatingRegion(double coefA) override;
    void subtractFluxRotatingRegion(double coefA) override;
    void multiply(double scalar) override;
    void setBufferFlux(Cell& cell) override;
    void buildCons(Phase** phase, Mixture* /*mixture*/) override;
    void buildPrim(Phase** phase, Mixture* /*mixture*/) override;
    void setToZero() override;
    void addNonCons(double /*coefA*/, const Cell* /*cell*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override {};
    void
    subtractNonCons(double /*coefA*/, const Cell* /*cell*/, const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override {};

    void addFluxSmooth1D(double coefA, const Coord& normal, Cell* cell) override;
    void substractFluxSmooth1D(double coefA, const Coord& normal, Cell* cell) override;

    void addSymmetricTerms(Phase** phases, Mixture* /*mixture*/, const double& r, const double& v) override;
    void prepSourceTermsGravity(const Coord& g) override;
    void prepSourceTermsHeating(const double& q) override;
    void prepSourceTermsMRF(Cell* cell, const Coord& omega) override;
    void prepSourceTermsPlaneWave(const Cell* cell, const double& sourceMomentum, const Coord& pulseDirection) override;

    //Moving Reference Frame
    void addNonConsMrfFlux(Phase** /*phases*/) override {};

    // Accessors
    //----------
    const Coord& getMomentum() const override { return m_momentum; };
    const double& getMassMix() const override { return m_mass; };
    const double& getEnergyMix() const override { return m_energ; };
    void setCons(const Flux* cons) override;

  protected:
    double m_mass;    //!< mass
    Coord m_momentum; //!< momentum
    double m_energ;   //!< total energy

  private:
    friend class ModEuler;
    // To modify if needed, example: to add a class APEViscosity, add friend class APEViscosity.
    friend class APEuler;
    friend class APEViscosity;
    friend class APEConductivity;
    friend class APETurbulenceSmagorinsky;
    friend class ModNEqDEM; //!< Because a variable of type FluxEuler is used in ModNEqDEM to store local fluxes and we need a write acces to FluxEuler members
    friend class RiemannSolverHLLC; //!< Because a variable of type FluxEuler is used in RiemannSolverHLLC to store local fluxes and we need a write acces to FluxEuler members
    friend class RiemannSolverExact; //!< Because a variable of type FluxEuler is used in RiemannSolverExact to store local fluxes and we need a write acces to FluxEuler members
};

#endif // FLUXEULER_H
