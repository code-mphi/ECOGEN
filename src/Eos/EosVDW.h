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

#ifndef EOSVDW_H
#define EOSVDW_H

#include "Eos.h"

//! \class     EosVDW
//! \brief     Class describing a Van Der Waals equation of state
class EosVDW : public Eos
{
  public:
    EosVDW(std::vector<std::string>& nameParameterEos, int& number);
    virtual ~EosVDW();

    //! \brief     Assign the values of the attributes for EosVDW from data defined in the code
    //! \param     name             string that contains the reduced name (sould be VDW)
    //! \param     parametersEos    vector (size depending on the Eos, 6 for VDW)
    //! \details   Are assigned the following attributes:  name, \f$ \Tc, \; Pc, \; r \; \gamma \; p_0 \; \rho_0 \f$. If the size of parameterEos \f$ \neq 3\f$  then the code aborts.
    virtual void assignParametersEos(std::string name, std::vector<double> parametersEos);

    //Constant methods (virtual because inherited from class Eos )
    
    //! \brief     Compute internal energy
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (T)
    //! \return    internal energy 
    //! \details   with  internal energy : \f$ \epsilon (\rho, T) = - (r T log(\frac{1}{\rho} - b)) - \rho a \f$
    virtual double computeEnergy(const double& density, const double& /*temperature*/) const;

    //! \brief     Compute pressure
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (T)
    //! \return    pressure 
    //! \details   with  pressure : \f$ p (\rho, T) = \frac{\rho r T}{1 - \rho b} - \rho^2 a \f$ //KS//To update
    virtual double computePressure(const double& density, const double& /*temperature*/) const;
    
    //! \brief     Compute temperature
    //! \param     density             density (\f$\rho\f$)
    //! \param     energy              energy (\f$\epsilon\f$)
    //! \return    temperature 
    //! \details   with  temperature : \f$ T (\rho, \epsilon) = \frac{e + \rho a}{-r log(\frac{1}{\rho} - b})} \f$
    virtual double computeTemperature(const double& density, const double& /*energy*/) const;

    //Partial derivatives
    //! \brief     Compute partial derivative dedrho
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (T)
    //! \return    dedrho 
    //! \details   with  dedrho : \f$ \frac{\partial \epsilon}{\partial \rho} (\rho, T) = \frac{r T}{\rho (1 - \rho b)} - a \f$ //KS//To update
    virtual double dedrho(const double& density, const double& /*temperature*/) const;
    //! \brief     Compute partial derivative dedrhoSecond
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (T)
    //! \return    dedrhoSecond 
    //! \details   with  dedrhoSecond : \f$ \frac{\partial^2 \epsilon}{\partial \rho^2} (\rho, T) = \frac{r T (1 - 2 \rho b)}{\rho^2 (1 - \rho b)^2} \f$ //KS//To update
    virtual double dedrhoSecond(const double& density, const double& /*temperature*/) const;

    //Checking
    //! \brief   Do nothing for VDW
    virtual void verifyAndCorrectDensityMax(double& density) const;

    //Get 
    //! \brief  Get the type that is to say the reduced name of the EOS in ECOGEN
    //! \return \f$ \ "VDW" \f$
    virtual TypeEOS getType() const { return TypeEOS::VDW; };

  private:
    double m_r;      //!< Universal constant of the fluid
    double m_a;      //!< Constant parameter of VDW EOS
    double m_b;      //!< Constant parameter of VDW EOS
    double m_gamma;  //!< Polytropic constant
    double m_p0;     //!< Reference pressure
    double m_rho0;   //!< Reference density
};

#endif // EOSVDW_H
