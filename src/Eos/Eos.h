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

#ifndef EOS_H
#define EOS_H

#include "../Maths/Coord.h"
#include "../Errors.h"
#include "../libTierces/tinyxml2.h"

//! \brief     Enumeration for the type of EOS (IG: ideal gas, SG: stiffened gas, NASG: Noble-Abel stiffened gas, VDW: Van Der Waals, Polynomial)
enum TypeEOS { IG, SG, NASG, VDW, Polynomial };

class Eos; //Predeclaration of class to include following .h

#include "../Tools.h"

//! \class     Eos
//! \brief     General class for Equation of State (EOS).
//! \details   This is a pure virtual class: can not be instantiated.
//! \details   This base class and its derived classes contain mainly computational functions to obtain thermodynamical values for each phase: pressure, temperature, internal energy...
//! \details   In this ECOGEN release these are the EOS available:
//! \details   ideal Gas (sub class EosIG), Stiffened (sub class EoSG), Noble-Abel Stiffened Gas (sub class EosNASG), Van der Waals (sub class EosVDW), polynomial (sub class EosPolynomial).
//! \details   Appropriate data for EOS leads to a large choice of gas, liquid or compressible solid.
class Eos
{
  public:
    Eos(int& number);
    virtual ~Eos();

    //Constant methods
    const std::string& getName() const { return m_name; };
    //! \brief See derived classes
    virtual TypeEOS getType() const { return TypeEOS::IG; };
    //! \brief  Return the number associated to the EOS
    //! \return  m_number
    const int& getNumber() const { return m_number; };

    //! \brief Read physical parameters (viscosity, thermal conductivity....)
    void readPhysicalParameter(tinyxml2::XMLNode *element);

    //! \brief    Compute the enthalpy of the phase
    //! \param    density         density (\f$\rho \f$)
    //! \param    pressure        pressure (p)
    //! \return   this -> h
    double computeEnthalpy(const double& density, const double& pressure) const;
    //! \brief    Compute the total enthalpy of the phase
    //! \param    density         density (\f$\rho \f$)
    //! \param    pressure        pressure (p)
    //! \param    velocity        velocity (m/s)
    //! \return   this -> hTotal
    //! \details  hTotal is computed as :\f$H_{total} =  \epsilon (\rho, p)  +  \frac{p}{\rho}  +  \frac{1}{2}u^2 \f$. 
    double computeTotalEnthalpy(const double& density, const double& pressure, const double& velX, const double& velY = 0., const double& velZ = 0.) const;
    //! \brief    Compute the total enthalpy of the phase
    //! \param    density         density (\f$\rho \f$)
    //! \param    pressure        pressure (p)
    //! \param    velocity        velocity vector (m/s)
    //! \return   this -> hTotal
    //! \details  hTotal is computed as :\f$H_{total} =  \epsilon (\rho, p)  +  \frac{p}{\rho}  +  \frac{1}{2}\bm{u}^2 \f$. 
    double computeTotalEnthalpy(const double& density, const double& pressure, const Coord& velocity) const;
    //! \brief  Return the dynamic viscosity of the fluid  
    //!return    \f$ \mu \f$ (Unit: Pa.s).
    const double& getMu() const { return m_mu; };
    //! \brief  get the thermal conductivity of the fluid
    //!return    \f$ \lambda \f$ (Unit: W/(m.K)).
    const double& getLambda() const { return m_lambda; };
    //! \brief  Assign the epsilonAlphaNull value for alphaNull option (alpha = 0 => epsilonAlphaNull != 0 ; alpha != 0 => epsilonAlphaNull = 0)
    //! \param     alphaNull       \f$\alpha = 0 \f$ (boolean value)
    void assignEpsilonForAlphaNull(bool alphaNull) const;

    //Gereral methods for all EOS

    //Virtual methods for child classes
    //! \brief See derived classes  
    virtual void assignParametersEos(std::string name, std::vector<double> parametersEos) = 0;
    //! \brief See derived classes 
    virtual double computeTemperature(const double& /*density*/, const double& /*pressure or energy*/) const { Errors::errorMessage("computeTemperature not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeEnergy(const double& /*density*/, const double& /*pressure or temperature*/) const { Errors::errorMessage("computeEnergy not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computePressure(const double& /*density*/, const double& /*energy or temperature*/) const { Errors::errorMessage("computePressure not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeDensity(const double& /*pressure*/, const double& /*temperature*/) const { Errors::errorMessage("computeDensity not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeSoundSpeed(const double& /*density*/, const double& /*pressure*/) const { Errors::errorMessage("computeSoundSpeed not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeInterfaceSoundSpeed(const double& /*density*/, const double& /*interfacePressure*/, const double& /*pressure*/) const { Errors::errorMessage("computeInterfaceSoundSpeed not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeAcousticImpedance(const double& /*density*/, const double& /*pressure*/) const { Errors::errorMessage("computeAcousticImpedance not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeDensityTimesInterfaceSoundSpeedSquare(const double& /*density*/, const double& /*interfacePressure*/, const double& /*pressure*/) const { Errors::errorMessage("computeDensityTimesInterfaceSoundSpeedSquare not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeEntropy(const double& /*temperature*/, const double& /*pressure*/) const { Errors::errorMessage("computeEntropy not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computePressureIsentropic(const double& /*initialPressure*/, const double& /*initialDensity*/, const double& /*finalDensity*/) const { Errors::errorMessage("computePressureIsentropic not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computePressureHugoniot(const double& /*initialPressure*/, const double& /*initialDensity*/, const double& /*finalDensity*/) const { Errors::errorMessage("computePressureHugoniot not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeDensityIsentropic(const double& /*initialPressure*/, const double& /*initialDensity*/, const double& /*finalPressure*/, double* /*drhodp*/ = 0) const { Errors::errorMessage("computeDensityIsentropic not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeDensityHugoniot(const double& /*initialPressure*/, const double& /*initialDensity*/, const double& /*finalPressure*/, double* /*drhodp*/ = 0) const { Errors::errorMessage("computeDensityHugoniot not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes
    virtual double computeDensityPfinal(const double& /*initialPressure*/, const double& /*initialDensity*/, const double& /*finalPressure*/, double* /*drhodp*/ = 0) const { Errors::errorMessage("computeDensityPfinal not yet programmed for EOS : " + m_name); return 0.; };

    //! \brief See derived classes 
    virtual double computeEnthalpyIsentropic(const double& /*initialPressure*/, const double& /*initialDensity*/, const double& /*finalPressure*/, double* /*dhdp*/ = 0) const { Errors::errorMessage("computeEnthalpyIsentropic not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeDensitySaturation(const double& /*pressure*/, const double& /*Tsat*/, const double& /*dTsatdP*/, double* /*drhodp*/ = 0) const { Errors::errorMessage("computeDensitySaturation not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double computeDensityEnergySaturation(const double& /*pressure*/, const double& /*rho*/, const double& /*drhodp*/, double* /*drhoedp*/ = 0) const { Errors::errorMessage("computeDensityEnergySaturation not available for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual void sendSpecialMixtureEos(double& /*gamPinfOverGamMinusOne*/, double& /*eRef*/, double& /*oneOverGamMinusOne*/, double& /*covolume*/) const { Errors::errorMessage("sendSpecialMixtureEos not yet programmed for EOS : " + m_name); };
    //! \brief See derived classes 
    virtual double vfpfh(const double& /*pressure*/, const double& /*enthalpy*/) const { Errors::errorMessage("vfpfh not yet programmed for EOS : " + m_name); return 0.; };

    //Partial derivatives
    //! \brief See derived classes 
    virtual double dvdpch(const double& /*pressure*/, const double& /*enthalpy*/) const { Errors::errorMessage("dvdpch not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double dvdhcp(const double& /*pressure*/) const { Errors::errorMessage("dvdhcp not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes
    virtual double drhodpcT(const double& /*pressure*/, const double& /*temperature*/) const { Errors::errorMessage("drhodpcT not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double dedrho(const double& /*density*/, const double& /*temperature*/) const { Errors::errorMessage("dedrho not yet programmed for EOS : " + m_name); return 0.; };
    //! \brief See derived classes 
    virtual double dedrhoSecond(const double& /*density*/, const double& /*temperature*/) const { Errors::errorMessage("dedrhoSecond not yet programmed for EOS : " + m_name); return 0.; };

    //Checking..
    virtual void verifyPressure(const double& /*pressure*/, const std::string& /*message*/ = " ") const { Errors::errorMessage("verifyPressure not yet programmed Eos : " + m_name); };
    //! \brief See derived classes 
    virtual void verifyAndModifyPressure(double& /*pressure*/) const { Errors::errorMessage("verifyAndModifyPressure not yet programmed Eos : "+ m_name); };
    //! \brief See derived classes 
    virtual void verifyAndCorrectDensityMax(const double& /*mass*/, double& /*alpha*/, double& /*density*/) const { Errors::errorMessage("verifyAndCorrectDensityMax not yet programmed Eos : "+ m_name); };
    //! \brief See derived classes 
    virtual void verifyAndCorrectDensityMax(double& /*density*/) const { Errors::errorMessage("verifyAndCorrectDensityMax not yet programmed Eos : "+ m_name); };

    //Get
    //! \brief See derived classes 
    virtual const double& getGamma() const { return Errors::defaultDouble; };
    //! \brief See derived classes 
    virtual const double& getPInf() const { return Errors::defaultDouble; };
    //! \brief See derived classes 
    virtual const double& getCv() const { return Errors::defaultDouble; };
    //! \brief See derived classes 
    virtual const double& getERef() const { return Errors::defaultDouble; };
    //! \brief See derived classes 
    virtual const double& getSRef() const { return Errors::defaultDouble; };

  protected:
    int m_number;       //!< Corresponding number of the equation of state
    std::string m_name; //!< Name of the equation of state

    double m_mu;        //!< Dynamic viscosity (kg/m/s or Pa.s)
    double m_lambda;    //!< Thermal conductivity (W/(m.K))
};

extern double epsilonAlphaNull;    //!< Epsilon value to avoid division by 0 when alpha = 0 is activated. If alpha = 0 desactivated, i.e. alpha != 0, then epsilonAlphaNull = 0.

#endif // EOS_H
