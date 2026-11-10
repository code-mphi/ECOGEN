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

#ifndef MIXNEQ_H
#define MIXNEQ_H

#include "../Mixture.h"

//! \class     MixNEq
//! \brief     Mixture variables for the non-equilibrium system of equations
class MixNEq : public Mixture
{
  public:
    MixNEq();
    ~MixNEq() override;

    void allocateAndCopyMixture(Mixture** mixture) override;
    void copyMixture(Mixture& mixture) override;
    double computeDensity(const double* alphak, const double* rhok) override;
    double computePressure(const double* alphak, const double* pk) override;
    double computeVelocity(const double* Yk, const double* uk) override;
    double computeInternalEnergy(const double* Yk, const double* ek) override;

    void computeMixtureVariables(Phase** vecPhase) override;
    void computeTotalEnergy(std::vector<QuantitiesAddPhys*>& /*vecGPA*/) override {};

    void localProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override;
    void reverseProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override;

    //Data printing
    int getNumberScalars() const override { return numberScalarsMixture; };
    int getNumberVectors() const override { return 1; };
    double returnScalar(const int& numVar) const override;
    Coord returnVector(const int& numVar) const override;
    std::string returnNameScalar(const int& numVar) const override;
    std::string returnNameVector(const int& numVar) const override;

    //Data reading
    void setScalar(const int& numVar, const double& value) override;
    void setVector(const int& numVar, const Coord& value) override;

    //Parallel
    int numberOfTransmittedVariables() const override { return 0; };
    void fillBuffer(double* /*buffer*/, int& /*counter*/) const override {};
    void fillBuffer(std::vector<double>& /*dataToSend*/) const override {};
    void getBuffer(double* /*buffer*/, int& /*counter*/) override {};
    void getBuffer(std::vector<double>& /*dataToReceive*/, int& /*counter*/) override {};

    //Second order
    void computeSlopesMixture(const Mixture& /*sLeft*/, const Mixture& /*sRight*/, const double& /*distance*/) override {};
    void setToZero() override {};
    void setToMax() override {};
    void extrapolate(const Mixture& /*slope*/, const double& /*distance*/) override {};
    void limitSlopes(const Mixture& /*slopeGauche*/, const Mixture& /*slopeDroite*/, Limiter& /*globalLimiter*/) override {};
    void setMin(const Mixture& /*mixture1*/, const Mixture& /*mixture2*/) override {};
    void setMax(const Mixture& /*mixture1*/, const Mixture& /*mixture2*/) override {};
    void computeGradientLimiter(const Limiter& /*globalLimiter*/,
                                const Mixture& /*mixture*/,
                                const Mixture& /*mixtureMin*/,
                                const Mixture& /*mixtureMax*/,
                                const Mixture& /*slope*/) override {};

    //Parallel second order
    int numberOfTransmittedSlopes() const override { return 0; };
    void fillBufferSlopes(double* /*buffer*/, int& /*counter*/) const override {};
    void getBufferSlopes(double* /*buffer*/, int& /*counter*/) override {};

    //Accessors
    const double& getDensity() const override { return m_density; };
    const double& getPressure() const override { return m_pressure; };
    const double& getU() const override { return m_velocity.getX(); };
    const double& getV() const override { return m_velocity.getY(); };
    const double& getW() const override { return m_velocity.getZ(); };
    const Coord& getVelocity() const override { return m_velocity; };
    Coord& getVelocity() override { return m_velocity; };
    const double& getEnergy() const override { return m_energy; };
    const double& getTotalEnergy() const override { return m_totalEnergy; };

    void setPressure(const double& p) override { m_pressure = p; };
    void setVelocity(const double& u, const double& v, const double& w) override { m_velocity.setXYZ(u, v, w); };
    void setVelocity(const Coord& vit) override { m_velocity = vit; };
    void setU(const double& u) override { m_velocity.setX(u); };
    void setV(const double& v) override { m_velocity.setY(v); };
    void setW(const double& w) override { m_velocity.setZ(w); };
    void setTotalEnergy(double& totalEnergy) override { m_totalEnergy = totalEnergy; };

    //Operators
    void changeSign() override {};
    void multiplyAndAdd(const Mixture& /*slopesMixtureTemp*/, const double& /*coeff*/) override {};
    void divide(const double& /*coeff*/) override {};

  protected:
    double m_density;     //!< mixture density
    double m_pressure;    //!< mixture pressure
    Coord m_velocity;     //!< mixture velocity
    double m_energy;      //!< mixture internal specific energy
    double m_totalEnergy; //!< mixture total specific energy
};

#endif // MIXNEQ_H
