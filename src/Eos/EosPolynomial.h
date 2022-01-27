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

#ifndef EOSPOLYNOMIAL_H
#define EOSPOLYNOMIAL_H

#include "Eos.h"

//! \class     EosPolynomial
//! \brief     Class describing a polynomial equation of state
class EosPolynomial : public Eos
{
  public:
    EosPolynomial(std::vector<std::string>& nameParameterEos, int& number);
    virtual ~EosPolynomial();

    //! \brief     Assign the values of the attributes for EosPolynomial from data defined in the code
    //! \param     name             string that contains the reduced name (should be Polynomial)
    //! \param     parametersEos    vector (size depending on the Eos, 3 for polynomial)
    //! \details   Are assigned the following attributes:  name, \f$ a, b, c\f$. If the size of parameterEos \f$ \neq 3\f$  then the code aborts.
    virtual void assignParametersEos(std::string name, std::vector<double> parametersEos);

    //Constant methods (virtual because inherited from class Eos)

    //! \brief     Compute pressure
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (T)
    //! \return    pressure
    //! \details   with  pressure : \f$ p (\epsilon = \frac{1}{\rho} - 1) = 4 a \epsilon^3 + 3 b \epsilon^2 + 2 c \epsilon \f$
    virtual double computePressure(const double& density, const double& /*temperature*/) const;

    //Partial derivatives
    //! \brief     Compute partial derivative dedrho
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (T)
    //! \return    dedrho 
    //! \details   with  dedrho : \f$ \frac{\partial \epsilon}{\partial \rho} (\rho) = \frac{p}{\rho^2} \f$
    virtual double dedrho(const double& density, const double& temperature) const;
    //! \brief     Compute partial derivative dedrhoSecond
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (T)
    //! \return    dedrhoSecond 
    //! \details   with  dedrhoSecond : \f$ \frac{\partial^2 \epsilon}{\partial \rho^2} (\rho) = \frac{1}{\rho^2} \frac{\partial p}{\partial \rho} - \frac{2}{\rho} \frac{\partial \epsilon}{\partial \rho} \f$
    virtual double dedrhoSecond(const double& density, const double& temperature) const;

    //Get 
    //! \brief  Get the type that is to say the reduced name of the EOS in ECOGEN
    //! \return \f$ \ "Polynomial" \f$
    virtual TypeEOS getType() const { return TypeEOS::Polynomial; };

  private:
    double m_a;   //!< Constant parameter of polynomial EOS
    double m_b;   //!< Constant parameter of polynomial EOS
    double m_c;   //!< Constant parameter of polynomial EOS
};

#endif // EOSPOLYNOMIAL_H
