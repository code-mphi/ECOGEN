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

#ifndef BOUNDCONDWALL_H
#define BOUNDCONDWALL_H

#include "BoundCond.h"


class BoundCondWall : public BoundCond
{
public:
  BoundCondWall(const BoundCondWall& Source, const int& lvl = 0); //Copy ctor (useful for AMR)
  BoundCondWall(int numPhysique, tinyxml2::XMLElement* element, std::string fileName);
  BoundCondWall(int numPhysique);
  virtual ~BoundCondWall();

  virtual void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces);
  virtual void solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax);
  virtual void solveRiemannTransportBoundary(Cell& /*cellLeft*/) const;

  virtual int whoAmI() const { return WALL; };
  virtual int whoAmIHeat() const { return m_heatCondition; }
  virtual bool isMRFWall() const { return m_isMRFWall; }

  virtual double getBoundaryHeatQuantity() const { return m_imposedHeatQuantity; }
  virtual Coord& getWallRotationalVelocityMRF() { return m_omegaWall; }

  //For AMR method
  virtual void creerCellInterfaceChild();  /*!< Create a child cell interface (not initialized) */

protected:
  TypeBCHeat m_heatCondition;   //!< Specific heat boundary condition, could be imposed temperature or flux density (default is adiabatic)
  double m_imposedHeatQuantity; //!< Imposed heat quantity on the wall. Depending on input (m_heatCondition) could be imposed temperature or flux density. This option requires conductivity additionnal physics to work
  bool m_isMRFWall;             //!< Flag to check if a special boundary treatment is required with MRF
  Coord m_omegaWall;            //!< Rotational wall velocity for MRF computation with viscous additionnal physic
private:
};

#endif // BOUNDCONDWALL_H
