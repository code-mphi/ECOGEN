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

#ifndef PHASESHALLOWWATER_H
#define PHASESHALLOWWATER_H

#include "../Phase.h"
#include "../../Eos/Eos.h"
#include <fstream>

//! \class     PhaseShallowWater
//! \brief     Phase variables for ShallowWater equations (single phase)
class PhaseShallowWater : public Phase
{
  public:
    PhaseShallowWater();
    //! \brief     Phase constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!            ex: <dataFluid height="1.0">
    //!                  <velocity x = "0." y = "0." / >
    //!                </dataFluid>
    //! \param     material           XML element to read for phase data
    //! \param     eos                EOS type
    //! \param     fileName           string name of readed XML file
    PhaseShallowWater(tinyxml2::XMLElement* material, Eos* eos, std::string fileName);
    ~PhaseShallowWater() override;

    void allocateAndCopyPhase(Phase** vecPhase) override;
    void copyPhase(Phase& vecPhase) override;
    void extendedCalculusPhase(const Coord& /*velocity*/) override;

    void localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal) override;
    void reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal) override;

    //Specific methods for data printing
    //----------------------------------
    int getNumberScalars() const override { return 2; };
    int getNumberVectors() const override { return 1; };
    double returnScalar(const int& numVar) const override;
    Coord returnVector(const int& numVar) const override;
    std::string returnNameScalar(const int& numVar) const override;
    std::string returnNameVector(const int& numVar) const override;

    //Specific method for reading from file
    //-------------------------------------
    void setScalar(const int& numVar, const double& value) override;
    void setVector(const int& numVar, const Coord& value) override;

    //Specific methods for parallel computing
    //---------------------------------------
    int numberOfTransmittedVariables() const override;
    void fillBuffer(double* buffer, int& counter) const override;
    void fillBuffer(std::vector<double>& dataToSend) const override;
    void getBuffer(double* buffer, int& counter, Eos** eos) override;
    void getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos) override;

    //Accessors
    //---------
    const double& getHeight() const override { return m_height; };
    const double& getPressure() const override { return m_pressure; };
    const double& getU() const override { return m_velocity.getX(); };
    const double& getV() const override { return m_velocity.getY(); };
    Coord& getVelocity() override { return m_velocity; };
    Eos* getEos() const override { return m_eos; };
    const double& getSoundSpeed() const override { return m_soundSpeed; };

    void setHeight(const double& height) override;
    void setPressure(double pressure) override;
    void setVelocity(const double& u, const double& v, const double& w) override;
    void setU(const double& u) override;
    void setV(const double& v) override;
    void setEos(Eos* eos) override;
    void setSoundSpeed(double soundSpeed) override;

  protected:
    double m_height;     //!< water height
    Coord m_velocity;    //!< velocity
    double m_pressure;   //!< pressure
    double m_soundSpeed; //!< speed of sound
    Eos* m_eos;          //!< pointer to equation of state
  private:
};

#endif // PHASESHALLOWWATER_H
