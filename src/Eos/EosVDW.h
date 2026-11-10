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
    ~EosVDW() override;

    //! \brief     Assign the values of the attributes for EosVDW from data defined in the code
    //! \param     name             string that contains the reduced name (sould be VDW)
    //! \param     parametersEos    vector (size depending on the Eos, 6 for VDW)
    //! \details   Are assigned the following attributes:  name, \f$ \Tc, \; Pc, \; r \; \gamma \; p_0 \; \rho_0 \f$. If the size of parameterEos \f$ \neq 3\f$  then the code aborts.
    void assignParametersEos(std::string name, std::vector<double> parametersEos) override;

    //Constant methods (virtual because inherited from class Eos )

    //! \brief     Compute internal energy
    //! \param     density             density (\f$\rho\f$)
    //! \param     pressure            pressure (p)
    //! \return    internal energy
    //! \details   with  internal energy : \f$ e (\rho, p)
    double computeEnergy(const double& density, const double& pressure) const override;

    //! \brief     Compute internal energy from density and temperature
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (T)
    //! \return    internal energy
    //! \details   with  internal energy : \f$ e (\rho, T)
    double computeEnergyFromDensityAndTemperature(const double& density, const double& temperature) const override;

    //! \brief     Compute pressure
    //! \param     density             density (\f$\rho\f$)
    //! \param     energy              energy (e)
    //! \return    pressure
    //! \details   with  pressure : \f$ p (\rho, e)
    double computePressure(const double& density, const double& energy) const override;

    //! \brief     Compute density
    //! \param     pressure            pressure (\f$p\f$)
    //! \param     temperature         temperature (K)
    //! \return    density
    //! \details   with  density : \f$ \rho (p, T)
    double computeDensity(const double& pressure, const double& temperature) const override;

    //! \brief     Compute pressure
    //! \param     density             density (\f$\rho\f$)
    //! \param     entropy             entropy (s)
    //! \return    pressure
    //! \details   with  pressure : \f$ p (\rho, s)
    double computePressureFromDensityAndEntropy(const double& density, const double& entropy) const override;

    //! \brief     Compute pressure
    //! \param     density             density (\f$\rho\f$)
    //! \param     temperature         temperature (\f$T\f$)
    //! \return    pressure
    //! \details   with  pressure : \f$ p (\rho, T)
    double computePressureFromDensityAndTemperature(const double& density, const double& temperature) const override;

    //! \brief     Compute temperature
    //! \param     density             density (\f$\rho\f$)
    //! \param     pressure            pressure (\f$p\f$)
    //! \return    temperature
    //! \details   with  temperature : \f$ T (\rho, p) = \frac{ (p + a * \rho^2) * ( 1 / \rho - b) }{ r } \f$
    double computeTemperature(const double& density, const double& pressure) const override;

    //! \brief     Compute squareSoundSpeed
    //! \param     density             density (\f$\rho\f$)
    //! \param     entropy             entropy (s)
    //! \return    squareSoundSpeed
    double computeSquareSoundSpeed(const double& density, const double& energy) const override;

    //! \brief     Compute Entropy from density and internal energy
    //! \param     density             density
    //! \param     energy              energy
    //! \return    entropy
    double computeEntropyFromDensityAndEnergy(const double& density, const double& energy) const override;

    //! \brief     Compute Entropy from density and pressure
    //! \param     density             density
    //! \param     pressure            pressure
    //! \return    entropy
    double computeEntropyFromDensityAndPressure(const double& density, const double& pressure) const override;

    //! \brief     Compute Entropy from density and temperature
    //! \param     density             density
    //! \param     pressure            pressure
    //! \return    entropy
    double computeEntropyFromDensityAndTemperature(const double& density, const double& pressure) const override;

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
    //! \brief   Do nothing for VDW
    void verifyAndCorrectDensityMax(double& density) const override;

    //Get
    //! \brief  Get the type that is to say the reduced name of the EOS in ECOGEN
    //! \return \f$ \ "VDW" \f$
    TypeEOS getType() const override { return TypeEOS::VDW; };

  private:
    double m_r;          //!< Universal constant of the fluid
    double m_a;          //!< Constant parameter of VDW EOS
    double m_b;          //!< Constant parameter of VDW EOS
    double m_gamma;      //!< Not the polytropic constant !
    double m_cv;         //!< Capacité calorifique à volume constant
    double m_kappa;      //!< Capillarity coefficient
    bool m_IsKappaConst; //!< Is the capillarity coefficient constant ?
};

#endif // EOSVDW_H
