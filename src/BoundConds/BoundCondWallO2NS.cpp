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

#include "BoundCondWallO2NS.h"

//****************************************************************************

BoundCondWallO2NS::BoundCondWallO2NS(const BoundCondWallO2NS& Source, const int& lvl) : BoundCondWall(Source, lvl)
{}

//****************************************************************************

BoundCondWallO2NS::BoundCondWallO2NS(int numPhysique, tinyxml2::XMLElement* element, std::string fileName) : BoundCondWall(numPhysique, element, fileName)
{}

//****************************************************************************

BoundCondWallO2NS::BoundCondWallO2NS(int numPhysique) : BoundCondWall(numPhysique)
{}

//****************************************************************************

BoundCondWallO2NS::~BoundCondWallO2NS()
{}

//****************************************************************************

void BoundCondWallO2NS::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondWallO2NS(*(this)));
}

//***********************************************************************

void BoundCondWallO2NS::solveRiemann(double& dtMax, Limiter& /*globalLimiter*/, Limiter& /*interfaceLimiter*/, Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/, Prim type)
{
  // Initialize buffer cells used in Riemann problem
  bufferCellLeft->copyVec(m_cellLeft->getPhases(type), m_cellLeft->getMixture(type), m_cellLeft->getTransports(type));

  // For NS extrapolation a scalar product between the distance 
  // Coord rij and the gradient is already done during CellO2NS::computeLocalSlopes 
  // and stored in slopesPhase/MixtureLocal1
  // Therefore, to stay compliant with cartesian 2nd order, distance is equal to 1 here
  double distanceLeft = 1.0; 

  // Left side of interface extrapolation
  m_cellLeft->computeLocalSlopes(*this); // Build slopesPhasesLocal1 = theta_i * rij . grad(Wi)
  for (int k = 0; k < numberPhases; k++) {
    bufferCellLeft->getPhase(k)->extrapolate(*slopesPhasesLocal1[k], distanceLeft); // Build Wij_lim = Wi + theta_i * rij . grad(Wi)
    bufferCellLeft->getPhase(k)->verifyAndCorrectPhase();
    bufferCellLeft->getPhase(k)->verifyAndCorrectDensityMax();
  }
  bufferCellLeft->getMixture()->extrapolate(*slopesMixtureLocal1, distanceLeft);
	for (int k = 0; k < numberTransports; k++) {
		bufferCellLeft->getTransport(k).extrapolate(slopesTransportLocal1[k], distanceLeft);
	}

  // Compute extended variables since projection change some of them (Phases, Mixture, AddPhys)
  bufferCellLeft->fulfillState();
  
  // Vector and tensor projections into reference frame attached to the face
  bufferCellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());

  // Riemann problem
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  this->solveRiemannBoundary(*bufferCellLeft, dxLeft, dtMax);
  // Handling of transport functions (m_Sm known: need to be called after Riemann solver)
  if (numberTransports > 0) { this->solveRiemannTransportBoundary(*bufferCellLeft); }

  // Flux projection into the absolute reference frame
  model->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
}

//***********************************************************************