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

#ifndef QAPTURBULENCE_H
#define QAPTURBULENCE_H

#include "../QuantitiesAddPhys.h"

//! \class     QAPTurbulence
//! \brief     General class for turbulent quantities
class QAPTurbulence : public QuantitiesAddPhys
{
  public:
    QAPTurbulence(AddPhys* addPhys);
    ~QAPTurbulence() override;

    void computeQuantities(Cell* cell) override;

    //Accessors
    void setGrad(const Coord& grad, const int& num = -1) override;                         //1:U, 2:V, 3:W
    const Coord& getGrad(const int& num = -1) const override { return m_grads[num - 1]; }; //1:U, 2:V, 3:W

    //! \brief     Set the turbulent viscosity in the cell
    void setTurbulentViscosity(const double& turbulentViscosity) override;
    //! \brief     Return the turbulent viscosity
    double getTurbulentViscosity() const override;

  protected:
    std::vector<Coord> m_grads;  //!< Vector of the velocity gradients of the cell in x-, y- and z-directions (ex: m_grads[0] = (du/dx, du/dy, du/dz))
    double m_turbulentViscosity; //!< Turbulent viscosity

  private:
};

#endif // QAPTURBULENCE_H
