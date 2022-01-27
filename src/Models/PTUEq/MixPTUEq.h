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

#ifndef MIXPTUEQ_H
#define MIXPTUEQ_H

#include "../Mixture.h"

//! \class     MixPTUEq
//! \brief     Mixture variables for pressure-temperature-velocity (mechanical and thermal equilibrium) system of equations
class MixPTUEq : public Mixture
{
  public:
    MixPTUEq();
    //! \brief     Mixture constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!           ex: <mixture>
    //!                 <dataMix pressure = "1.e5" temperature = "300.">
    //!                 <velocity x = "0." y = "0." z = "0." />
    //!               </mixture>
    //! \param     state          XML element to read for mixture data
    //! \param     fileName       string name of readed XML file
    MixPTUEq(tinyxml2::XMLElement* state, std::string fileName);
    virtual ~MixPTUEq();

    virtual void allocateAndCopyMixture(Mixture** mixture);
    virtual void copyMixture(Mixture &mixture);
    virtual double computeDensity(const double* alphak, const double* rhok);
    virtual double computePressure(const double* alphak, const double* pk);
    virtual double computePressure(double* masses, const double& mixInternalEnerg, Phase** phases);
    virtual double computeTemperature(double* masses, const double& pressure, Phase** phases);
    virtual double computeInternalEnergy(const double* Yk, const double* ek);
    virtual double computeFrozenSoundSpeed(const double* Yk, const double* ck);
    
    //Specific thermodynamical evolutions
    virtual double computeTemperatureIsentrope(const double* Yk, const double& p0, const double& T0, const double& p, double* dTdp = 0);
    virtual double computeEnthalpyIsentrope(const double* Yk, const double& p0, const double& T0, const double& p, double* dhdp = 0);
    virtual double computeVolumeIsentrope(const double* Yk, const double& p0, const double& T0, const double& p, double* dvdp = 0);

    virtual void computeMixtureVariables(Phase** vecPhase);
    virtual void internalEnergyToTotalEnergy(std::vector<QuantitiesAddPhys*>& vecGPA);
    virtual void totalEnergyToInternalEnergy(std::vector<QuantitiesAddPhys*>& vecGPA);

    virtual void localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal);
    virtual void reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal);

    //Data printing
    virtual int getNumberScalars() const { return 3; };
    virtual int getNumberVectors() const { return 1; };
    virtual double returnScalar(const int& numVar) const;
    virtual Coord returnVector(const int& numVar) const;
    virtual std::string returnNameScalar(const int& numVar) const;
    virtual std::string returnNameVector(const int& numVar) const;

    //Data reading
    virtual void setScalar(const int& numVar, const double& value);
    virtual void setVector(const int& numVar, const Coord& value);

    //Parallel
    virtual int numberOfTransmittedVariables() const;
    virtual void fillBuffer(double* buffer, int& counter) const;
    virtual void getBuffer(double* buffer, int& counter);

    //Second order
    virtual void computeSlopesMixture(const Mixture &sLeft, const Mixture &sRight, const double& distance);
    virtual void setToZero();
    virtual void extrapolate(const Mixture &slope, const double& distance);
    virtual void limitSlopes(const Mixture &slopeGauche, const Mixture &slopeDroite, Limiter& globalLimiter);

    //Parallel second order
    virtual int numberOfTransmittedSlopes() const;
    virtual void fillBufferSlopes(double* buffer, int& counter) const;
    virtual void getBufferSlopes(double* buffer, int& counter);

    //Accessors
    virtual const double& getDensity() const { return m_density; };
    virtual const double& getPressure() const { return m_pressure; };
    virtual const double& getTemperature() const { return m_temperature; };
    virtual const double& getU() const { return m_velocity.getX(); };
    virtual const double& getV() const { return m_velocity.getY(); };
    virtual const double& getW() const { return m_velocity.getZ(); };
    virtual const Coord& getVelocity() const { return m_velocity; };
    virtual Coord& getVelocity() { return m_velocity; };
    virtual const double& getEnergy() const { return m_energy; };
    virtual const double& getTotalEnergy() const { return m_totalEnergy; };
    virtual const double& getMixSoundSpeed() const { return m_PTUEqSoundSpeed; };

    virtual void setPressure(const double& p);
    virtual void setVelocity(const double& u, const double& v, const double& w);
    virtual void setVelocity(const Coord& vit);
    virtual void setU(const double& u);
    virtual void setV(const double& v);
    virtual void setW(const double& w);
    virtual void setTotalEnergy(double& totalEnergy);

    //Operators
    virtual void changeSign();
    virtual void multiplyAndAdd(const Mixture &slopesMixtureTemp, const double& coeff);
    virtual void divide(const double& coeff);

  private:
    double m_density;              //!< mixture density
    double m_pressure;             //!< mixture pressure
    Coord m_velocity;              //!< mixture velocity
    double m_temperature;          //!< mixture temperature
    double m_energy;               //!< mixture internal specific energy
    double m_totalEnergy;          //!< mixture total specific energy
    double m_PTUEqSoundSpeed;      //!< mixture thermal equilibrium sound speed
};

#endif // MIXPTUEQ_H
