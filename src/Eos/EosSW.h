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

#ifndef EOSSW_H
#define EOSSW_H

#include "Eos.h"

//! \class     EosSW
//! \brief     Class describing an ideal gas equation of state
class EosSW : public Eos
{
  public:
    EosSW(std::vector<std::string>& nameParameterEos, int& number);
    ~EosSW() override;

    //! \brief     Assign the values of the attributes for EosSW from data defined in the code
    //! \param     name             string that contains the reduced name (sould be SW)
    //! \param     parametersEos    vector (size depending on the Eos, 2 for SW)
    //! \details   Assign 'name' and 'g' attributes. If the size of parameterEos \f$ \neq 1\f$  then the code aborts.
    void assignParametersEos(std::string name, std::vector<double> parametersEos) override;

    //! \brief     Compute pressure
    //! \param     height   phase height (h)
    //! \return    pressure
    //! \details   with  pressure : \f$  p(h)  = 0.5 g h^2 \f$
    double computePressure(const double& height) const override;

    //! \brief     Compute sound speed
    //! \param     height     phase height (h)
    //! \return    soundSpeed
    //! \details   with  soundSpeed : \f$  c(h)  = \sqrt{ g h } \f$
    double computeSoundSpeed(const double& height) const override;

    //Getters

  private:
    double m_gravity; //!< Gravity
};

#endif // EOSSW_H
