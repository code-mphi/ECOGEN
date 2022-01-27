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

#ifndef PHASEUEQ_H
#define PHASEUEQ_H

#include "../Phase.h"
#include "../../Eos/Eos.h"
#include <fstream>

//! \class     PhaseUEq
//! \brief     Phase variables for the velocity-equilibrium system of equations
class PhaseUEq : public Phase
{
  public:
    PhaseUEq();
    //! \brief     Phase constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!           ex:  <dataFluid alpha="0.5" density="1.0" pressure="1.e5"/> 
    //! \param     material           XML element to read for phase data
    //! \param     eos                EOS pointer to compute thermodynamic variables
    //! \param     fileName           string name of readed XML file
    PhaseUEq(tinyxml2::XMLElement* material, Eos* eos, std::string fileName);
    virtual ~PhaseUEq();

    virtual void allocateAndCopyPhase(Phase** vecPhase);
    virtual void copyPhase(Phase& phase);
    virtual void extendedCalculusPhase(const Coord& /*velocity*/);
    virtual void computeMassFraction(const double& density);

    virtual void localProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) {};
    virtual void reverseProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) {};

    //Specific methods for data printing
    //----------------------------------
    virtual int getNumberScalars() const { return 5; };
    virtual int getNumberVectors() const { return 0; };
    virtual double returnScalar(const int& numVar) const;
    virtual Coord returnVector(const int& /*numVar*/) const { return 0; };
    virtual std::string returnNameScalar(const int& numVar) const;
    virtual std::string returnNameVector(const int& /*numVar*/) const { return 0; };

    //Specific method for reading from file
    //-------------------------------------
    virtual void setScalar(const int& numVar, const double& value);
    
    //Specific methods for parallel computing
    //---------------------------------------
    virtual int numberOfTransmittedVariables() const;
    virtual void fillBuffer(double* buffer, int& counter) const;
    virtual void fillBuffer(std::vector<double>& dataToSend) const;
    virtual void getBuffer(double* buffer, int& counter, Eos** eos);
    virtual void getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos);

    //Specific methods for second order
    //---------------------------------
    virtual void computeSlopesPhase(const Phase& sLeft, const Phase& sRight, const double& distance);
    virtual void setToZero();
    virtual void extrapolate(const Phase& slope, const double& distance);
    virtual void limitSlopes(const Phase& slopeGauche, const Phase& slopeDroite, Limiter& globalLimiter, Limiter& volumeFractionLimiter);

    //Specific methods for parallele computing at second order
    //--------------------------------------------------------
    virtual int numberOfTransmittedSlopes() const;
    virtual void fillBufferSlopes(double* buffer, int& counter) const;
    virtual void getBufferSlopes(double* buffer, int& counter);

    //Verifications
    //-------------
    virtual void verifyPhase(const std::string& message = "") const;
    virtual void verifyAndCorrectPhase();
    virtual void verifyAndCorrectDensityMax(const double& mass);
    virtual void verifyAndCorrectDensityMax();

    //Accessors
    //---------
    virtual const double& getAlpha() const { return m_alpha; };
    virtual const double& getMassFraction() const { return m_Y; };
    virtual const double& getDensity() const { return m_density; };
    virtual const double& getPressure() const { return m_pressure; };
    virtual const double& getY() const { return m_Y; };
    virtual const double& getU() const { return Errors::defaultDouble; };
    virtual const double& getV() const { return Errors::defaultDouble; };
    virtual const double& getW() const { return Errors::defaultDouble; };
    virtual Coord& getVelocity() { return Coord::defaultCoordNonConst; };
    virtual const Coord& getVelocity() const { return Coord::defaultCoord; };
    virtual Eos* getEos() const { return m_eos; };
    virtual const double& getEnergy() const { return m_energy; };
    virtual const double& getSoundSpeed() const { return m_soundSpeed; };
    virtual double getTemperature() const { return m_eos->computeTemperature(m_density, m_pressure); }; //KS//Modify with m_temperature directly?

    virtual void setAlpha(double alpha);
    virtual void setDensity(double density);
    virtual void setPressure(double pressure);
    virtual void setVelocity(const double& /*u*/, const double& /*v*/, const double& /*w*/) {};
    virtual void setVelocity(const Coord& /*vit*/) {};
    virtual void setU(const double& /*u*/) {};
    virtual void setV(const double& /*v*/) {};
    virtual void setW(const double& /*w*/) {};
    virtual void setEos(Eos* eos);
    virtual void setEnergy(double energy);
    virtual void setSoundSpeed(double soundSpeed);

    //Operators
    //---------
    virtual void changeSign();
    virtual void multiplyAndAdd(const Phase& slopesPhasesTemp, const double& coeff);
    virtual void divide(const double& coeff);

  protected:
    double m_alpha;           //!< phase volume fraction
    double m_density;         //!< phase specific mass
    double m_pressure;        //!< phase pressure
    double m_Y;               //!< phase mass fraction
    double m_temperature;     //!< phase temperature
    Eos* m_eos;               //!< pointer to phase equation of state
    double m_energy;          //!< phase internal energy
    double m_soundSpeed;      //!< phase speed of sound
};

#endif // PHASEUEQ_H
