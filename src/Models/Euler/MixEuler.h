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

#ifndef MIXEULER_H
#define MIXEULER_H

#include "../Mixture.h"

//! \class     MixEuler
//! \brief     Mixture variables for Euler equations (single phase)
class MixEuler : public Mixture
{
  public:
    MixEuler();
    ~MixEuler() override;

    void allocateAndCopyMixture(Mixture** mixture) override;
    void copyMixture(Mixture& /*mixture*/) override {};
    double computeDensity(const double* /*alphak*/, const double* /*rhok*/) override { return 0.; };
    double computePressure(const double* /*alphak*/, const double* /*pk*/) override { return 0.; };
    double computeInternalEnergy(const double* /*Yk*/, const double* /*ek*/) override { return 0.; };
    double computeFrozenSoundSpeed(const double* /*Yk*/, const double* /*ck*/) override { return 0.; };

    void computeMixtureVariables(Phase** /*vecPhase*/) override {};
    void computeTotalEnergy(std::vector<QuantitiesAddPhys*>& /*vecGPA*/) override {};

    void localProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override {};
    void reverseProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override {};

    //Specific methods for data printing
    //----------------------------------
    int getNumberScalars() const override { return 0; };
    int getNumberVectors() const override { return 0; };
    double returnScalar(const int& /*numVar*/) const override { return 0.; };
    Coord returnVector(const int& /*numVar*/) const override { return 0; };
    std::string returnNameScalar(const int& /*numVar*/) const override { return Errors::defaultString; };
    std::string returnNameVector(const int& /*numVar*/) const override { return Errors::defaultString; };

    //Specific methods for parallel computing
    //---------------------------------------
    int numberOfTransmittedVariables() const override { return 0; };
    void fillBuffer(double* /*buffer*/, int& /*counter*/) const override {};
    void fillBuffer(std::vector<double>& /*dataToSend*/) const override {};
    void getBuffer(double* /*buffer*/, int& /*counter*/) override {};
    void getBuffer(std::vector<double>& /*dataToReceive*/, int& /*counter*/) override {};

    //Specific methods for second order
    //---------------------------------
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

    //Specific methods for parallele computing at second order
    //--------------------------------------------------------
    int numberOfTransmittedSlopes() const override { return 0; };
    void fillBufferSlopes(double* /*buffer*/, int& /*counter*/) const override {};
    void getBufferSlopes(double* /*buffer*/, int& /*counter*/) override {};

    //Accessors
    //---------
    const double& getDensity() const override { return Errors::defaultDouble; };
    const double& getPressure() const override { return Errors::defaultDouble; };
    const double& getU() const override { return Errors::defaultDouble; };
    const double& getV() const override { return Errors::defaultDouble; };
    const double& getW() const override { return Errors::defaultDouble; };
    const Coord& getVelocity() const override { return Coord::defaultCoord; };
    Coord& getVelocity() override { return Coord::defaultCoordNonConst; };
    const double& getEnergy() const override { return Errors::defaultDouble; };
    const double& getTotalEnergy() const override { return Errors::defaultDouble; };
    const double& getFrozenSoundSpeed() const override { return Errors::defaultDouble; };
    const double& getWoodSoundSpeed() const override { return Errors::defaultDouble; };

    void setPressure(const double& /*p*/) override {};
    void setVelocity(const double& /*u*/, const double& /*v*/, const double& /*w*/) override {};
    void setVelocity(const Coord& /*vit*/) override {};
    void setU(const double& /*u*/) override {};
    void setV(const double& /*v*/) override {};
    void setW(const double& /*w*/) override {};
    void setTotalEnergy(double& /*totalEnergy*/) override {};

    //Operators
    //---------
    void changeSign() override {};
    void multiplyAndAdd(const Mixture& /*slopesMixtureTemp*/, const double& /*coeff*/) override {};
    void divide(const double& /*coeff*/) override {};

  protected:
  private:
};

#endif // MIXEULER_H
