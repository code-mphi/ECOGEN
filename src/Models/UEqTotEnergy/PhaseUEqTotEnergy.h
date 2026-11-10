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

#ifndef PHASEUEQTOTENERGY_H
#define PHASEUEQTOTENERGY_H

#include "../Phase.h"
#include "../../Eos/Eos.h"
#include <fstream>

//! \class     PhaseUEqTotEnergy
//! \brief     Phase variables for the velocity-equilibrium system of equations
class PhaseUEqTotEnergy : public Phase
{
  public:
    PhaseUEqTotEnergy();
    //! \brief     Phase constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!           ex:  <dataFluid alpha="0.5" density="1.0" pressure="1.e5"/>
    //! \param     material           XML element to read for phase data
    //! \param     eos                EOS pointer to compute thermodynamic variables
    //! \param     fileName           string name of readed XML file
    PhaseUEqTotEnergy(tinyxml2::XMLElement* material, Eos* eos, std::string fileName);
    ~PhaseUEqTotEnergy() override;

    void allocateAndCopyPhase(Phase** vecPhase) override;
    void copyPhase(Phase& phase) override;
    void extendedCalculusPhase(const Coord& velocity) override;
    void computeMassFraction(const double& density) override;

    void localProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override {};
    void reverseProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) override {};

    //Specific methods for data printing
    //----------------------------------
    int getNumberScalars() const override { return 5; };
    int getNumberVectors() const override { return 0; };
    double returnScalar(const int& numVar) const override;
    Coord returnVector(const int& /*numVar*/) const override { return 0; };
    std::string returnNameScalar(const int& numVar) const override;
    std::string returnNameVector(const int& /*numVar*/) const override { return Errors::defaultString; };

    //Specific method for reading from file
    //-------------------------------------
    void setScalar(const int& numVar, const double& value) override;

    //Specific methods for parallel computing
    //---------------------------------------
    int numberOfTransmittedVariables() const override;
    void fillBuffer(double* buffer, int& counter) const override;
    void fillBuffer(std::vector<double>& dataToSend) const override;
    void getBuffer(double* buffer, int& counter, Eos** eos) override;
    void getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos) override;

    //Specific methods for second order
    //---------------------------------
    void computeSlopesPhase(const Phase& sLeft, const Phase& sRight, const double& distance) override;
    void setToZero() override;
    void extrapolate(const Phase& slope, const double& distance) override;
    void limitSlopes(const Phase& slopeGauche, const Phase& slopeDroite, Limiter& globalLimiter, Limiter& volumeFractionLimiter) override;

    //Specific methods for parallele computing at second order
    //--------------------------------------------------------
    int numberOfTransmittedSlopes() const override;
    void fillBufferSlopes(double* buffer, int& counter) const override;
    void getBufferSlopes(double* buffer, int& counter) override;

    //Verifications
    //-------------
    void verifyPhase(const std::string& message = "") const override;
    void verifyAndCorrectPhase() override;
    void verifyAndCorrectDensityMax(const double& mass) override;
    void verifyAndCorrectDensityMax() override;

    //Accessors
    //---------
    const double& getAlpha() const override { return m_alpha; };
    const double& getDensity() const override { return m_density; };
    const double& getPressure() const override { return m_pressure; };
    const double& getY() const override { return m_Y; };
    const double& getU() const override { return Errors::defaultDouble; };
    const double& getV() const override { return Errors::defaultDouble; };
    const double& getW() const override { return Errors::defaultDouble; };
    Coord& getVelocity() override { return Coord::defaultCoordNonConst; };
    const Coord& getVelocity() const override { return Coord::defaultCoord; };
    Eos* getEos() const override { return m_eos; };
    const double& getTotalEnergy() const override { return m_totEnergy; };
    const double& getSoundSpeed() const override { return m_soundSpeed; };
    double getTemperature() const override { return m_temperature; };

    void setAlpha(double alpha) override;
    void setDensity(double density) override;
    void setPressure(double pressure) override;
    void setVelocity(const double& /*u*/, const double& /*v*/, const double& /*w*/) override {};
    void setVelocity(const Coord& /*vit*/) override {};
    void setU(const double& /*u*/) override {};
    void setV(const double& /*v*/) override {};
    void setW(const double& /*w*/) override {};
    void setEos(Eos* eos) override;
    void setTotalEnergy(double totalEnergy) override;
    void setTotalEnergy(const double& energy, const Coord& vel) override;
    void setSoundSpeed(double soundSpeed) override;

    //Operators
    //---------
    void changeSign() override;
    void multiplyAndAdd(const Phase& slopesPhasesTemp, const double& coeff) override;
    void divide(const double& coeff) override;

  protected:
    double m_alpha;       //!< phase volume fraction
    double m_density;     //!< phase specific mass
    double m_pressure;    //!< phase pressure
    double m_Y;           //!< phase mass fraction
    double m_temperature; //!< phase temperature
    Eos* m_eos;           //!< pointer to phase equation of state
    double m_totEnergy;   //!< phase specific total energy
    double m_soundSpeed;  //!< phase speed of sound
};

#endif // PHASEUEQTOTENERGY_H
