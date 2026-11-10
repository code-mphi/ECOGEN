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

#include "BoundCondPeriodic.h"

//****************************************************************************

BoundCondPeriodic::BoundCondPeriodic() {}

//****************************************************************************

BoundCondPeriodic::BoundCondPeriodic(const BoundCondPeriodic& Source, const int& lvl) : BoundCond(Source, lvl) {}

//****************************************************************************

BoundCondPeriodic::BoundCondPeriodic(int numPhysique) : BoundCond(numPhysique) {}

//****************************************************************************

BoundCondPeriodic::~BoundCondPeriodic() {}

//****************************************************************************

void BoundCondPeriodic::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondPeriodic(*(this)));
}

//****************************************************************************

void BoundCondPeriodic::initialize(Cell* cellLeft, Cell* cellRight)
{
  m_cellLeft  = cellLeft;
  m_cellRight = cellRight;
}

//***********************************************************************

void BoundCondPeriodic::solveRiemann(double& dtMax,
                                     Limiter& /*globalLimiter*/,
                                     Limiter& /*interfaceLimiter*/,
                                     Limiter& /*globalVolumeFractionLimiter*/,
                                     Limiter& /*interfaceVolumeFractionLimiter*/,
                                     Prim /*type*/)
{
  //Vector and tensor projections on geometric reference frame of the face
  m_cellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
  m_cellRight->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());

  //Riemann problem
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  double dxRight(m_cellRight->getElement()->getLCFL());
  dxLeft  = dxLeft * std::pow(2., (double)m_lvl);
  dxRight = dxRight * std::pow(2., (double)m_lvl);
  model->solveRiemannIntern(*m_cellLeft, *m_cellRight, dxLeft, dxRight, dtMax);
  //Handling of transport functions (m_Sm known: need to be called after Riemann solver)
  if (numberTransports > 0) {
    model->solveRiemannTransportIntern(*m_cellLeft, *m_cellRight);
  }

  //Flux projection on absolute reference frame
  model->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
  m_cellLeft->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
  m_cellRight->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
}

//****************************************************************************
