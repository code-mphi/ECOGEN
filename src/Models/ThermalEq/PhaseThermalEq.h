//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

#ifndef PHASETHERMALEQ_H
#define PHASETHERMALEQ_H

//! \file      PhaseThermalEq.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      May 04 2018

#include "../Phase.h"
#include "../../Eos/Eos.h"
#include <fstream>

//! \class     PhaseThermalEq
//! \brief     Phase variables for ThermalEq system of equations (mechanical and thermal equilibrium)
class PhaseThermalEq : public Phase
{
  public:
    PhaseThermalEq();
    //! \brief     Phase constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!           ex:  <dataFluid alpha="0.5"/>
    //! \param     material           XML element to read for phase data
    //! \param     fileName           string name of readed XML file
    PhaseThermalEq(tinyxml2::XMLElement *material, Eos *eos, std::string fileName);
    virtual ~PhaseThermalEq();

    virtual void allocateAndCopyPhase(Phase **vecPhase);
    virtual void copyPhase(Phase &vecPhase);
    virtual void extendedCalculusPhase(const Coord &velocity);

    virtual void localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal) {};
    virtual void reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal) {};

    //Specific methods for data printing
    //----------------------------------
    virtual int getNumberScalars() const { return 2; };
    virtual int getNumberVectors() const { return 0; };
    virtual double returnScalar(const int &numVar) const;
    virtual Coord returnVector(const int &numVar) const { return 0; };
    virtual std::string returnNameScalar(const int &numVar) const;
    virtual std::string returnNameVector(const int &numVar) const { return 0; };

    //Specific method for reading from file
    //-------------------------------------
    virtual void setScalar(const int &numVar, const double &value);
    
    //Specific methods for parallel computing
    //---------------------------------------
    virtual int numberOfTransmittedVariables() const;
    virtual void fillBuffer(double *buffer, int &counter) const;
    virtual void getBuffer(double *buffer, int &counter, Eos **eos);

    //Specific methods for second order
    //---------------------------------
    virtual void computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double &distance);
    virtual void setToZero();
    virtual void extrapolate(const Phase &slope, const double &distance);
    virtual void limitSlopes(const Phase &slopeGauche, const Phase &slopeDroite, Limiter &globalLimiter, Limiter &volumeFractionLimiter);

    //Specific methods for parallele computing at second order
    //--------------------------------------------------------
	virtual int numberOfTransmittedSlopes() const;
	virtual void fillBufferSlopes(double *buffer, int &counter) const;
	virtual void getBufferSlopes(double *buffer, int &counter);

    //Verifications
    //-------------
    virtual void verifyPhase(const std::string &message = "") const;
    virtual void verifyAndCorrectPhase();

    //Accessors
    //---------
    virtual double getAlpha() const;
    virtual double getDensity() const;
    virtual double getPressure() const;
    virtual double getU() const { return 0.; };
    virtual double getV() const { return 0.; };
    virtual double getW() const { return 0.; };
    virtual Coord getVelocity() const { return 0; };
    virtual Eos* getEos() const;
    virtual double getEnergy() const;
    virtual double getSoundSpeed() const;
    virtual double getTotalEnergy() const;
    virtual double getTemperature() const;

    virtual void setAlpha(double alpha);
    virtual void setDensity(double density);
    virtual void setPressure(double pressure);
    virtual void setVelocity(const double &u, const double &v, const double &w) {};
    virtual void setVelocity(const Coord &vit) {};
    virtual void setU(const double &u) {};
    virtual void setV(const double &v) {};
    virtual void setW(const double &w) {};
    virtual void setEos(Eos *eos);
    virtual void setEnergy(double energie);
    virtual void setSoundSpeed(double soundSpeed);
    virtual void setTotalEnergy(double totalEnergy);

    //Operators
    //---------
    virtual void changeSign();
    virtual void multiplyAndAdd(const Phase &slopesPhasesTemp, const double &coeff);
    virtual void divide(const double &coeff);

  protected:
    double m_alpha;           //!< phase volume fraction
    double m_density;         //!< phase specific mass
    double m_pressure;        //!< phase pressure
    Eos *m_eos;               //!< pointer to phase equation of state
    double m_energie;         //!< phase internal energy
    double m_soundSpeed;      //!< phase speed of sound
    double m_totalEnergy;     //!< phase total energy
  private:
};

#endif // PHASETHERMALEQ_H
