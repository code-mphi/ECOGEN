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

#ifndef BOUNDCONDINLETINJTEMP_H
#define BOUNDCONDINLETINJTEMP_H

#include "BoundCond.h"

class BoundCondInletInjTemp : public BoundCond
{
  public:
    BoundCondInletInjTemp(int numPhysique, tinyxml2::XMLElement* element, const int& numbPhases, Eos** eos, std::string fileName = "Unknown file");
    BoundCondInletInjTemp(const BoundCondInletInjTemp& Source, const int& lvl = 0); // Copy ctor (useful for AMR)
    ~BoundCondInletInjTemp() override;

    void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces) override;
    void solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax) override;

    int whoAmI() const override { return TypeBC::INLETINJTEMP; };
    void printInfo() override;

    // For AMR method
    void creerCellInterfaceChild() override; /*!< Create a child cell interface (uninitialized) */

  private:
    double m_m0;   //!< Inflow specific massflow
    double* m_Tk0; //!< Inflow temperature
    double* m_ak0; //!< Inflow volume fraction
};

#endif
