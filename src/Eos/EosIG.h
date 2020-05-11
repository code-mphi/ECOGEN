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

#ifndef EOSIG_H
#define EOSIG_H

//! \file      EosIG.h
//! \author    F. Petitpas, K. Schmidmayer, E. Daniel
//! \version   1.1
//! \date      June 5 2019

#include "Eos.h"

//! \class     EosIG
//! \brief     Class describing an ideal gas equation of state
class EosIG : public Eos
{
    public:
        EosIG();
        EosIG(std::vector<std::string> &nameParameterEos, int& number);
        virtual ~EosIG();

		//! \brief    assign the values of the attributes for EosIG from data defined in the code
		//! \param     name             string that contains the reduced name (sould be IG)
		//! \param     parametersEos    vector (size depending on the Eos, 4 for IG)
		//! \details  Are assigned the following attributes:  name, \f$ \gamma, \; c_v, \; e_{ref},  s_{ref}\f$. If the size of parameterEos \f$ \neq 4\f$  then the code aborts.
        virtual void assignParametersEos(std::string name, std::vector<double> parametersEos);
        //Constant methods (virtual because inherited from class Eos )
		//! \brief     Compute temperature
		//! \param     density             density (\f$\rho \f$)
		//! \param     pressure            pressure (p)
		//! \return    temperature  
		//! \details   with temperature : \f$  T(\rho, p)  =   \frac{p}{\rho c_v (\gamma-1)}\f$
        virtual double computeTemperature(const double &density,const double &pressure) const;
		//! \brief     Compute internal energy
		//! \param     density             density (\f$\rho\f$)
		//! \param     pressure            pressure (p)
		//! \return    internal energy 
		//! \details  with  internal energy : \f$  \epsilon (\rho, p)  = \frac{p}{\rho \ (\gamma-1)} \ + \ \epsilon_{ref}\f$
        virtual double computeEnergy(const double &density,const double &pressure) const;
		//! \brief     Compute pressure
		//! \param     density             density (\f$\rho\f$)
		//! \param     energy              internal energy (\f$\epsilon\f$)
		//! \return    pressure 
		//! \details  with  pressure : \f$  p(\rho, \epsilon)  = (\gamma-1)\rho(\epsilon  - \epsilon_{ref}) \f$
        virtual double computePressure(const double &density,const double &energy) const;
		//! \brief     Compute density
		//! \param     pressure            pressure (p)
		//! \param     temperature         temperature (T)
		//! \return    density  
		//! \details  with  density : \f$  \rho(p, T)  =  \frac{p}{T \ c_v \ (\gamma-1)}\ \f$
        virtual double computeDensity(const double &pressure, const double &temperature) const;
		//! \brief     Compute sound speed
		//! \param     density            density (\f$\rho\f$)
		//! \param     pressure           pressure (p)
		//! \return    soundSpeed
		//! \details  with  soundSpeed : \f$  c(p, \rho)  = \sqrt{  \frac{\gamma \ p}{\rho}} \f$
        virtual double computeSoundSpeed(const double &density,const double &pressure) const;
		//! \brief     Compute entropy
		//! \param     temperature        temperature (T)
		//! \param     pressure           pressure (p)
		//! \return    entropy
		//! \details  with s : \f$  s(T, p)  = c_v \ ln \left( \frac{T^\gamma}  { p^\frac{\gamma}{(\gamma-1)} } \right)+s_{ref} \f$
        virtual double computeEntropy(const double &temperature, const double &pressure) const;
		//! \brief     Compute phase pressure along an isentropic path
		//! \param     initialPressure    initial pressure (\f$ p_i \f$)
		//! \param     initialDensity     initial density (\f$ \rho_i \f$)
		//! \param     finalDensity       final density
		//! \return    finalPressure
		//! \details  with finalPressure :  \f$  p_f  = p_i  \left( \frac{\rho_f}{\rho_i} \right) ^\gamma   \f$
        virtual double computePressureIsentropic(const double &initialPressure, const double &initialDensity, const double &finalDensity) const; 
		//! \brief     Compute  pressure along the Hugoniot curve
		//! \param     initialPressure    initial pressure (\f$ p_i \f$)
		//! \param     initialDensity     initial density (\f$ \rho_i \f$)
		//! \param     finalDensity       final density (\f$ \rho_f \f$)
		//! \return    finalPressure  
		//! \details    with finalPressure : \f$ p_f = p_i \frac{  (\gamma+1)*\rho_f - (\gamma-1)*\rho_i  }{  (\gamma+1)*\rho_i-(\gamma-1)*\rho_f  }  \f$
        virtual double computePressureHugoniot(const double &initialPressure, const double &initialDensity, const double &finalDensity) const;
		//! \brief     Compute  density along an isentropic path AND its derivature versus the pressure at constant entropy
		//! \param     initialPressure    initial pressure \f$ (p_i) \f$
		//! \param     initialDensity     initial density   \f$ (\rho_i ) \f$
		//! \param     finalPressure      final pressure
		//! \return    finalDensity AND *drhodp (OPTIONAL)
		//! \detail    with finalDensity : \f$ \rho_f = \rho_i \left( \frac{p_f}{ p_i}\right) ^\frac{1}{\gamma} \f$ AND with   *drhodp :  \f$  \left.  \frac{\partial \rho}{\partial p} \right)_s   = \frac{\rho_f}{ \gamma \ p_f}  \f$
        virtual double computeDensityIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp=0) const;
		//! \brief     Compute  density along the Hugoniot curve path AND its derivature versus the pressure along Hugoniot curve
		//! \param     initialPressure    initial pressure (\f$ p_i \f$)
		//! \param     initialDensity     initial density (\f$ \rho_i \f$)
		//! \param     finalPressure      final pressure (\f$ p_f \f$)
		//! \return    finalDensity  AND  *drhodp (OPTIONAL)
		//! \detail    with finalDensity : \f$ \rho_f =\rho_i\frac{ (\gamma + 1)*p_f + (\gamma - 1)*p_i}{(\gamma - 1)p_f +(\gamma + 1)p_i) }\f$ AND with  *drhodp :  \f$  \left.  \frac{\partial \rho}{\partial p} \right)_Hug   =  \frac{4\gamma \rho_i \ p_i}{ \left( (\gamma - 1)p_f +(\gamma + 1)p_i \right) ^2}   \f$
        virtual double computeDensityHugoniot(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp = 0) const;
        //! \brief     Compute  enthalpy at the end of an isentropic path AND the enthalpy derivature versus the pressure
		//! \param     initialPressure	  Initial pressure (\f$ p_i \f$ at the beginning of the isentropic path)
		//! \param     initialDensity	  Initial density (\f$ \rho_i \f$ at the beginning of the isentropic path)
		//! \param     finalPressure	  final Pressure (\f$ p_f \f$ at the end of the isentropic path)
		//! \return    finalEnthalpy AND *dhdp (OPTIONAL)
		//! \detail with finalEnthalpy :  \f$ h(p_f,{\rho_f})= \frac{\gamma}{(\gamma-1)}\frac{p_f}{\rho_f} +e_{ref} \f$. The density \f$ \rho_f \f$ is fisrt computed along an isentropic path (see computeDensityIsentropic EosIG method ). AND with *dhdp : \f$ \frac{dh}{dp}=\frac{\gamma}{(\gamma-1)}\frac{ \rho_f - p_f  \frac{d\rho}{dp}  }{\rho_f^2} \f$. The term \f$\frac{d\rho}{dp}\f$ is the derivative obtained along an isentropic path.

        virtual double computeDensityPfinal(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp = 0) const;

        virtual double computeEnthalpyIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *dhdp = 0) const;
		//! \brief     Compute density on the saturation curve at the saturation temperature  AND  its derivative versus the pressure at constant density (OPTIONAL)
		//! \param     pressure          pressure (p)
		//! \param     Tsat        		 saturation temperature (\f$T_{sat}\f$)
		//! \param     dTsatdP  		 derivative of the saturation temperature versus pressure (\f$\frac{d T_{sat}}{d p}\f$)
		//! \return    rho 	AND  *drhodp (OPTIONAL)
		//! \detail  with  rho : \f$  \rho(p, T_{sat})  =  \frac{p}{T_{sat} \ c_v \ (\gamma-1)}\ \f$ AND \f$\frac{d \rho}{d p} = \left.  \frac{\partial \rho}{\partial p} \right)_{T_{sat}} + \left.  \frac{\partial \rho}{\partial T_{sat} } \right)_p \frac{d T_{sat}}{d p}  =\frac{(\gamma-1)c_v\left(T_{sat} - p\frac{d T_{sat}}{d p}\right) }{\left((\gamma-1)c_vT_{sat}\right)^2 }\f$
        virtual double computeDensitySaturation(const double &pressure, const double &Tsat, const double &dTsatdP, double *drhodp = 0) const;

		//! \brief     Compute  the volumic internal energy at the saturation AND  its derivative versus the pressure at constant density (OPTIONAL)
		//! \param     pressure			pressure (p)
		//! \param     density			density (\f$ \rho \f$)
		//! \param     drhodp			derivative of the density versus pressure (\f$\frac{d \rho}{d p}\f$ at the saturation state)
		//! \return    rhoe AND *drhoedp (OPTIONAL)
		//! \detail  with  rhoe : \f$  \rho\epsilon(p, \rho)  =  \frac{p}{(\gamma-1)} +\rho\epsilon_{ref}\f$ AND \f$ \frac{d \rho\epsilon}{d p} =\left.  \frac{\partial \rho\epsilon}{\partial p} \right)_{\rho} + \left.  \frac{\partial \rho\epsilon}{\partial\rho  } \right)_p \frac{d \rho}{d p}  = \frac{1}{(\gamma-1)} +\frac{d \rho}{d p}\epsilon_{ref} \f$
        virtual double computeDensityEnergySaturation(const double &pressure, const double &rho, const double &drhodp, double *drhoedp = 0) const;

		//! \brief    send specific values of the parameters useful for mixture EOS based on Ideal Gas and Stiffened Gas
		//! \param     gamPinfOverGamMinusOne	\f$  \frac{\gamma p_{\infty}}{(\gamma-1)} = 0 \;for\; Ideal\; Gas\f$ 		
		//! \param     eRef						\f$	  \epsilon_{ref} \f$ 
		//! \param     oneOverGamMinusOne		\f$  \frac{1}{(\gamma-1)}\f$ 	  
        virtual void sendSpecialMixtureEos(double &gamPinfOverGamMinusOne, double &eRef, double &oneOverGamMinusOne, double &covolume) const;
		//! \brief   compute the specific volume with the pressure and the enthalpy
		//! \param   pressure		pressure (p)  		
		//! \param   enthalpy       enthalpy (h)
		//! \return    specific volume (volume per mass unit) :  \f$ v= (h - \epsilon_{ref})\frac{(\gamma-1)}{\gamma p} \f$ from the definition of enthalpy \f$ h=\epsilon +pv \f$ 
		virtual double vfpfh(const double &pressure, const double &enthalpy) const;

		//Partial derivatives 
		//! \brief    compute the partial derivative of the specific volume versus pressure at constant enthalpy
		//! \param   pressure		pressure (p)  
		//! \param   enthalpy       enthalpy (h)
		//! \return  derivative : \f$ \left. \frac{\partial v}{\partial p} \right)_h = - \frac{(\gamma-1)(h-\epsilon_{ref})}{\gamma p^2} \f$	
        virtual double dvdpch(const double &pressure, const double &enthalpy) const;
		//! \brief   compute the partial derivative of the specific volume versus enthalpy at constant pressure
		//! \param   pressure		pressure (p)  
		//! \param   enthalpy       enthalpy (h)
		//! \return  derivative : \f$ \left. \frac{\partial v}{\partial h} \right)_p =  \frac{(\gamma-1)}{\gamma p} \f$
        virtual double dvdhcp(const double &pressure, const double &enthalpy) const;

    //Checking
		//! \brief  add a message error if the pressure is too low and under a limit value equal to \f$ 10^{-15} \f$ 
		//! \param   pressure		pressure (p)
		//! \param   message        ""
		//! \return  add a message in the  \f$ errors \f$ vector. The message is \f$ " :\;too\;low\;pressure \;in\; EosIG\;"\f$
        virtual void verifyPressure(const double &pressure, const std::string &message = "") const;
		//! \brief  modify the pressure if its value is under a limit value equal to \f$ 10^{-15} \f$ 
		//! \param   pressure		pressure (p)
		//! \return  pressure : \f$ if \; p < 10^{-15} \; then \;  p= 10^{-15} \f$
		virtual void verifyAndModifyPressure(double &pressure) const;

		//Mod

		//! \brief  get in a vector the values of the parameters of the Ideal Gas EOS
		//! \param   data		 
		//! \return  data[0]= \f$ \gamma \f$ ; data[1]= \f$ c_v \f$ ; data[2]= \f$ \epsilon_{ref} \f$
		   virtual void sendInfo(double *&data) const;

    //Get 
		//! \brief  get the adiabatic exponent of the fluid 
		//!return  m_gamma :  \f$ \gamma \f$
        virtual const double& getGamma() const { return m_gamma; };
		//! \brief  get the volume calorific energy of the fluid
		//!return     \f$ c_v \f$
        virtual const double& getCv() const { return m_cv; };
		//! \brief  get the energy of reference of the fluid
		//!return     \f$ \epsilon_{ref} \f$
        virtual const double& getERef() const { return m_eRef; };
		//! \brief  get the entropy of reference of the fluid
		//!return     \f$ \ s_{ref} \f$
        virtual const double& getSRef() const { return m_sRef; };
		//! \brief  get the type that is to say the  reduced name of the EOS in ECOGEN :
		//!return     \f$ \ "IG" \f$
        virtual std::string getType() const { return "IG"; };

    protected:
    private:
        double m_gamma;  //!< Adiabatic exponent of the fluid
        double m_cv;     //!< Volume calorific energy of the fluid
        double m_eRef;   //!< Energy of reference of the fluid
        double m_sRef;   //!< Entropy of the fluid
};

#endif // EOSIG_H
