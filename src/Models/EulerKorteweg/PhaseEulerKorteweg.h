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

#ifndef PHASEEULERKORTEWEG_H
#define PHASEEULERKORTEWEG_H

#include "../Phase.h"
#include "../../Eos/Eos.h"
#include <fstream>

//! \class     PhaseEulerKorteweg
//! \brief     Phase variables for Augmented Euler--Korteweg equations (single phase)
class PhaseEulerKorteweg : public Phase
{
  public:
    PhaseEulerKorteweg();
    //! \brief     Phase constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!            ex: <dataFluid density = "10.0">
    //!                  <velocity x = "1000." y = "1000." z = "0." / >
    //!                </dataFluid>
    //! \param     material           XML element to read for phase data
    //! \param     eos                EOS pointer to compute thermodynamic variables
    //! \param     fileName           string name of readed XML file
    PhaseEulerKorteweg(tinyxml2::XMLElement* material, Eos* eos, std::string fileName);
    virtual ~PhaseEulerKorteweg();

    virtual void allocateAndCopyPhase(Phase** vecPhase);
    virtual void copyPhase(Phase& vecPhase);
    virtual void extendedCalculusPhase(const Coord& /*velocity*/) {};

    virtual void localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal);
    virtual void reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal);

    //Specific methods for data printing
    //----------------------------------
    virtual int getNumberScalars() const { return 4; };
    virtual int getNumberVectors() const { return 2; };
    virtual double returnScalar(const int& numVar) const;
    virtual Coord returnVector(const int& numVar) const;
    virtual std::string returnNameScalar(const int& numVar) const;
    virtual std::string returnNameVector(const int& numVar) const;

    //Specific method for reading from file
    //-------------------------------------
    virtual void setScalar(const int& numVar, const double& value);
    virtual void setVector(const int& numVar, const Coord& value);
    
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
    virtual void limitSlopes(const Phase& slopeGauche, const Phase& slopeDroite, Limiter& globalLimiter, Limiter& /*volumeFractionLimiter*/);

    //Specific methods for parallele computing at second order
    //--------------------------------------------------------
    virtual int numberOfTransmittedSlopes() const;
    virtual void fillBufferSlopes(double* buffer, int& counter) const;
    virtual void getBufferSlopes(double* buffer, int& counter);

    //Verifications
    //-------------
    virtual void verifyPhase(const std::string& message = "") const;
    virtual void verifyAndCorrectPhase();
    virtual void verifyAndCorrectDensityMax();

    //Accessors
    //---------
    virtual const double& getAlpha() const { return Errors::defaultDouble; };
    virtual const double& getDensity() const { return m_density; };
    virtual const double& getOmega() const { return m_omega; };
    virtual const double& getEta() const { return m_eta; };
    virtual const double& getPressure() const { return m_pressure; };
    virtual const double& getU() const { return m_velocity.getX(); }; 
    virtual const double& getV() const { return m_velocity.getY(); };
    virtual const double& getW() const { return m_velocity.getZ(); };
    virtual Coord& getVelocity() { return m_velocity; };
    virtual const Coord& getVelocity() const { return m_velocity; };
    virtual const double& getVectorPX() const { return m_vectorP.getX(); };
    virtual const double& getVectorPY() const { return m_vectorP.getY(); };
    virtual const double& getVectorPZ() const { return m_vectorP.getZ(); };
    virtual Coord& getVectorP() { return m_vectorP; };
    virtual const Coord& getVectorP() const { return m_vectorP; };
    virtual Eos* getEos() const { return m_eos; };
    
    virtual void setDensity(double density);
    virtual void setOmega(const double& omega);
    virtual void setEta(const double& eta);
    virtual void setPressure(double pressure);
    virtual void setVelocity(const double& u, const double& v, const double& w);
    virtual void setVelocity(const Coord& vit);
    virtual void setU(const double& u);
    virtual void setV(const double& v);
    virtual void setW(const double& w);
    virtual void setVectorP(const double& Px, const double& Py, const double& Pz);
    virtual void setVectorP(const Coord& vecP);
    virtual void setVectorPX(const double& Px);
    virtual void setVectorPY(const double& Py);
    virtual void setVectorPZ(const double& Pz);
    virtual void setEos(Eos* eos);
    
    //Operators
    //---------
    virtual void changeSign();
    virtual void multiplyAndAdd(const Phase& slopesPhasesTemp, const double& coeff);
    virtual void divide(const double& coeff);

  protected:
    double m_density;          //!< Specific mass
    double m_omega;            //!< Time derivative of eta
    double m_eta;              //!< Analogue of density
    double m_pressure;         //!< Pressure, so far it is only for information
    Coord m_velocity;          //!< Velocity
    Coord m_vectorP;           //!< Gradient of eta
    Eos* m_eos;                //!< Pointer to equation of state
};

#endif // PHASEEULERKORTEWEG_H
