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

#ifndef BOUNDCONDMOVINGWALL_H
#define BOUNDCONDMOVINGWALL_H

#include "BoundCondWall.h"

class BoundCondMovingWall : public BoundCondWall
{
  public:
    BoundCondMovingWall(const BoundCondMovingWall& Source, const int& lvl = 0); //Copy ctor (useful for AMR)
    BoundCondMovingWall(int numPhysique, tinyxml2::XMLElement* element, std::string fileName);
    BoundCondMovingWall(int numPhysique);
    ~BoundCondMovingWall() override;

    void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces) override;
    void solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax) override;
    void solveRiemannTransportBoundary(Cell& /*cellLeft*/) const override;

    int whoAmI() const override { return TypeBC::MOVINGWALL; };

    //Accessors
    double getWallVelocity() override;

    //For AMR method
    void creerCellInterfaceChild() override; /*!< Create a child cell interface (not initialized) */

  protected:
    double m_amplitude;       //!< Amplitude of the moving wall following oscillatory behaviour [m/s]
                              //!< If frequency is null and offset is pi/2 it corresponds to a constant piston velocity
    double m_frequency;       //!< Frequency of the moving wall following oscillatory behaviour [Hz]
    double m_frequencyOffset; //!< Frequency offset of the sinus [Hz]
};

#endif // BOUNDCONDMOVINGWALL_H
