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
    ~EosPolynomial() override;

    //! \brief     Assign the values of the attributes for EosPolynomial from data defined in the code
    //! \param     name             string that contains the reduced name (should be Polynomial)
    //! \param     parametersEos    vector (size depending on the Eos, 3 for polynomial)
    //! \details   Are assigned the following attributes:  name, \f$ a, b, c\f$. If the size of parameterEos \f$ \neq 3\f$  then the code aborts.
    void assignParametersEos(std::string name, std::vector<double> parametersEos) override;

    //Constant methods (virtual because inherited from class Eos)

    //! \brief     Compute internal energy
    //! \param     density             density (\f$\rho\f$)
    //! \param     pressure            pressure (\f$p\f$)
    //! \return    internal energy
    //! \details   with  internal energy : \f$ e (\rho, p)
    double computeEnergy(const double& density, const double& /*pressure*/) const override;

    //! \brief     Compute Entropy from density and internal energy
    //! \param     density             density (\f$\rho\f$)
    //! \param     energy              energy (\f$e\f$)
    //! \return    entropy
    double computeEntropyFromDensityAndEnergy(const double& density, const double& /*energy*/) const override;

    //! \brief     Compute Entropy from density and pressure
    //! \param     density             density (\f$\rho\f$)
    //! \param     pressure            pressure (\f$p\f$)
    //! \return    entropy
    double computeEntropyFromDensityAndPressure(const double& density, const double& /*pressure*/) const override;

    //! \brief     Compute pressure
    //! \param     density             density (\f$\rho\f$)
    //! \param     entropy             entropy (\f$s\f$)
    //! \return    pressure
    //! \details   with  pressure : \f$ p (\rho, s)
    double computePressureFromDensityAndEntropy(const double& density, const double& /*entropy*/) const override;

    //! \brief     Compute pressure
    //! \param     density             density (\f$\rho\f$)
    //! \param     energy              energy (\f$e\f$)
    //! \return    pressure
    //! \details   with  pressure : \f$ p (\rho, e)
    double computePressure(const double& density, const double& /*energy*/) const override;

    //! \brief     Compute temperature
    //! \param     density             density (\f$\rho\f$)
    //! \param     pressure            pressure (\f$p\f$)
    //! \return    temperature
    double computeTemperature(const double& density, const double& pressure) const override;

    //! \brief     Compute SquareSoundSpeed
    //! \param     density             density (\f$\rho\f$)
    //! \param     entropy             entropy (\f$s\f$)
    //! \return    Square Sound Speed
    double computeSquareSoundSpeed(const double& density, const double& entropy) const override;

    //! \brief     Compute the Cappilary coefficient (kappa) from density
    //! \param     density             density (\f$\rho\f$)
    //! \return    kappa
    double kappa(const double& density) const override;

    //! \brief     Compute the First derivative of the Cappilary coefficient (kappa) with respect to the density
    //! \param     density             density (\f$\rho\f$)
    //! \return    kappaPrime
    double kappaPrime(const double& density) const override;

    //! \brief     Compute the Second derivative of the Cappilary coefficient (kappa) with respect to the density
    //! \param     density             density (\f$\rho\f$)
    //! \return    kappaSecond
    double kappaSecond(const double& density) const override;

    //Checking
    //! \brief   Do nothing for Polynomial
    void verifyAndCorrectDensityMax(double& /*density*/) const override {};

    //Get
    //! \brief  Get the type that is to say the reduced name of the EOS in ECOGEN
    //! \return \f$ \ "Polynomial" \f$
    TypeEOS getType() const override { return TypeEOS::Polynomial; };

  private:
    double m_a; //!< Constant parameter of polynomial EOS
    double m_b; //!< Constant parameter of polynomial EOS
    double m_c; //!< Constant parameter of polynomial EOS
};

#endif // EOSPOLYNOMIAL_H
