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

#include "APEViscosity.h"

//***********************************************************************

APEViscosity::APEViscosity() {}

//***********************************************************************

APEViscosity::APEViscosity(int& numberQPA, Eos** eos)
{
  m_mu     = eos[0]->getMu();
  m_numQPA = numberQPA++;
}

//***********************************************************************

APEViscosity::~APEViscosity() {}

//***********************************************************************

void APEViscosity::addQuantityAddPhys(Cell* cell) { cell->getVecQuantitiesAddPhys().push_back(new QAPViscosity(this)); }

//***********************************************************************

void APEViscosity::solveFluxAddPhys(CellInterface* cellInterface)
{
  // Copy velocities and gradients of left and right cells
  m_velocityLeft  = cellInterface->getCellLeft()->getPhase(0)->getVelocity();
  m_velocityRight = cellInterface->getCellRight()->getPhase(0)->getVelocity();

  m_gradULeft  = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(1);
  m_gradURight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(1);
  m_gradVLeft  = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(2);
  m_gradVRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(2);
  m_gradWLeft  = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(3);
  m_gradWRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(3);

  m_normal   = cellInterface->getFace()->getNormal();
  m_tangent  = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Projection on orientation axes attached to the edge of velocities and gradients
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_velocityRight.localProjection(m_normal, m_tangent, m_binormal);

  m_tensorLeft.setTensorByLines(m_gradULeft, m_gradVLeft, m_gradWLeft);
  m_tensorRight.setTensorByLines(m_gradURight, m_gradVRight, m_gradWRight);
  m_tensorLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_tensorRight.localProjection(m_normal, m_tangent, m_binormal);
  m_tensorLeft.tensorToCoords(m_gradULeft, m_gradVLeft, m_gradWLeft);
  m_tensorRight.tensorToCoords(m_gradURight, m_gradVRight, m_gradWRight);

  this->solveFluxViscosityInner(m_velocityLeft, m_velocityRight, m_gradULeft, m_gradURight, m_gradVLeft, m_gradVRight, m_gradWLeft, m_gradWRight, m_mu);

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APEViscosity::solveFluxAddPhysBoundary(CellInterface* cellInterface)
{
  //DEV// BC Injection, Tank, Outflow to do

  // Copy velocities and gradients of left and right cells
  m_velocityLeft = cellInterface->getCellLeft()->getPhase(0)->getVelocity();

  m_gradULeft = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(1);
  m_gradVLeft = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(2);
  m_gradWLeft = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(3);

  m_normal   = cellInterface->getFace()->getNormal();
  m_tangent  = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Projection on orientation axes attached to the edge of velocities and gradients
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);

  m_tensorLeft.setTensorByLines(m_gradULeft, m_gradVLeft, m_gradWLeft);
  m_tensorLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_tensorLeft.tensorToCoords(m_gradULeft, m_gradVLeft, m_gradWLeft);

  // Distances cells/cell interfaces for weighting on the flux
  double distLeft = cellInterface->getCellLeft()->distance(cellInterface);

  switch (cellInterface->whoAmI()) {
  case TypeBC::NONREFLECTING:
  case TypeBC::OUTLETPRESSURE:
  case TypeBC::INLETINJSTAGSTATE:
  case TypeBC::INLETTANK:
  case TypeBC::INLETINJTEMP:
    this->solveFluxViscosityNonReflecting(m_velocityLeft, m_gradULeft, m_gradVLeft, m_gradWLeft, m_mu);
    break;
  case TypeBC::WALL:
    if (!cellInterface->isMRFWall()) {
      this->solveFluxViscosityWall(m_velocityLeft, m_mu, distLeft);
    }
    else {
      Coord velocityWall(0.);
      velocityWall = cellInterface->getWallRotationalVelocityMRF().cross(cellInterface->getFace()->getPos());
      velocityWall.localProjection(m_normal, m_tangent, m_binormal);
      velocityWall.setX(0.); // Fluid not crossing boundary
      this->solveFluxViscosityWall(m_velocityLeft, m_mu, distLeft, velocityWall);
    }
    break;
  case TypeBC::SYMMETRY:
    this->solveFluxViscositySymmetry(m_gradULeft, m_gradVLeft, m_gradWLeft, m_mu);
    break;
  case TypeBC::PERIODIC:
    m_velocityRight = cellInterface->getCellRight()->getPhase(0)->getVelocity();

    m_gradURight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(1);
    m_gradVRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(2);
    m_gradWRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(3);

    m_tensorRight.setTensorByLines(m_gradURight, m_gradVRight, m_gradWRight);
    m_tensorRight.localProjection(m_normal, m_tangent, m_binormal);
    m_tensorRight.tensorToCoords(m_gradURight, m_gradVRight, m_gradWRight);

    this->solveFluxViscosityInner(
      m_velocityLeft, m_velocityRight, m_gradULeft, m_gradURight, m_gradVLeft, m_gradVRight, m_gradWLeft, m_gradWRight, m_mu);
    break;
  default:
    this->solveFluxViscosityOther();
  }

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APEViscosity::solveFluxViscosityInner(const Coord& velocityLeft,
                                           const Coord& velocityRight,
                                           const Coord& gradULeft,
                                           const Coord& gradURight,
                                           const Coord& gradVLeft,
                                           const Coord& gradVRight,
                                           const Coord& gradWLeft,
                                           const Coord& gradWRight,
                                           const double& mu) const
{
  //Extraction of data
  double uL(0.), vL(0.), wL(0.), uR(0.), vR(0.), wR(0.);
  double dudxL(0.), dudyL(0.), dudzL(0.), dudxR(0.), dudyR(0.), dudzR(0.);
  double dvdxL(0.), dvdyL(0.), dvdxR(0.), dvdyR(0.);
  double dwdxL(0.), dwdzL(0.), dwdxR(0.), dwdzR(0.);
  uL = velocityLeft.getX();
  vL = velocityLeft.getY();
  wL = velocityLeft.getZ();
  uR = velocityRight.getX();
  vR = velocityRight.getY();
  wR = velocityRight.getZ();

  dudxL = gradULeft.getX();
  dudyL = gradULeft.getY();
  dudzL = gradULeft.getZ();
  dudxR = gradURight.getX();
  dudyR = gradURight.getY();
  dudzR = gradURight.getZ();

  dvdxL = gradVLeft.getX();
  dvdyL = gradVLeft.getY();
  dvdxR = gradVRight.getX();
  dvdyR = gradVRight.getY();

  dwdxL = gradWLeft.getX();
  dwdzL = gradWLeft.getZ();
  dwdxR = gradWRight.getX();
  dwdzR = gradWRight.getZ();

  //Data of the cell interface
  double u(0.), v(0.), w(0.);
  double dudx(0.), dudy(0.), dudz(0.);
  double dvdx(0.), dvdy(0.);
  double dwdx(0.), dwdz(0.);

  u = (uL + uR) / 2.;
  v = (vL + vR) / 2.;
  w = (wL + wR) / 2.;

  dudx = (dudxL + dudxR) / 2.;
  dudy = (dudyL + dudyR) / 2.;
  dudz = (dudzL + dudzR) / 2.;

  dvdx = (dvdxL + dvdxR) / 2.;
  dvdy = (dvdyL + dvdyR) / 2.;

  dwdx = (dwdxL + dwdxR) / 2.;
  dwdz = (dwdzL + dwdzR) / 2.;

  //Writing of viscous terms on each equation of fluxBuffEuler
  //AF//TODO// Precompute mu/3.
  static_cast<FluxEuler*>(fluxBuff)->m_mass = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setX(-mu / 3. * (4. * dudx - 2. * (dvdy + dwdz)));
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setY(-mu * (dvdx + dudy));
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setZ(-mu * (dwdx + dudz));
  static_cast<FluxEuler*>(fluxBuff)->m_energ = -mu * (1. / 3. * u * (4 * dudx - 2. * (dvdy + dwdz)) + (dvdx + dudy) * v + (dwdx + dudz) * w);
}

//***********************************************************************

void APEViscosity::solveFluxViscosityNonReflecting(
  const Coord& velocityLeft, const Coord& gradULeft, const Coord& gradVLeft, const Coord& gradWLeft, const double& muLeft) const
{
  this->solveFluxViscosityInner(velocityLeft, velocityLeft, gradULeft, gradULeft, gradVLeft, gradVLeft, gradWLeft, gradWLeft, muLeft);
}

//***********************************************************************

void APEViscosity::solveFluxViscosityWall(const Coord& velocityLeft, const double& muLeft, const double& distLeft, Coord& velocityWall) const
{
  //Computed the gradients locally because of the particular condition at the wall
  double dudx = velocityWall.getX() - velocityLeft.getX() / distLeft;
  double dvdx = velocityWall.getY() - velocityLeft.getY() / distLeft;
  double dwdx = velocityWall.getZ() - velocityLeft.getZ() / distLeft;

  //Writing of viscous terms on each equation of fluxBuffEuler
  static_cast<FluxEuler*>(fluxBuff)->m_mass = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setX(-muLeft / 3. * 4. * dudx);
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setY(-muLeft * dvdx);
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setZ(-muLeft * dwdx);
  static_cast<FluxEuler*>(fluxBuff)->m_energ = 0.;
}

//***********************************************************************

void APEViscosity::solveFluxViscositySymmetry(const Coord& gradULeft, const Coord& gradVLeft, const Coord& gradWLeft, const double& mu) const
{
  //Extraction of data
  //Note that dudy, dudz, dvdx and dwdx are nulls
  double dudx(0.), dvdy(0.), dwdz(0.);
  dudx = gradULeft.getX();
  dvdy = gradVLeft.getY();
  dwdz = gradWLeft.getZ();

  //Writing of viscous terms on each equation of fluxBuffEuler
  static_cast<FluxEuler*>(fluxBuff)->m_mass = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setX(-mu / 3. * (4. * dudx - 2. * (dvdy + dwdz)));
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxEuler*>(fluxBuff)->m_energ = 0.;
}

//***********************************************************************

void APEViscosity::solveFluxViscosityOther() const
{
  //Not manage at the moment, just an example
  //std::cout << "Viscous boundary not managed" << std::endl;

  // To avoid bug when not managed
  static_cast<FluxEuler*>(fluxBuff)->m_mass     = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_momentum = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_energ    = 0.;
}

//***********************************************************************

void APEViscosity::communicationsAddPhys(const int& dim, const int& lvl)
{
  parallel.communicationsVector(Variable::QPA, dim, lvl, m_numQPA, 1); //m_gradU
  parallel.communicationsVector(Variable::QPA, dim, lvl, m_numQPA, 2); //m_gradV
  parallel.communicationsVector(Variable::QPA, dim, lvl, m_numQPA, 3); //m_gradW
}

//***********************************************************************

void APEViscosity::computeTurbulenceParameters() const
{
  double global_sum_u2v2w2(0.);

  MPI_Allreduce(&sum_u2v2w2, &global_sum_u2v2w2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  long l_nCells = static_cast<long>(nCells);
  long global_nCells(0);
  MPI_Allreduce(&l_nCells, &global_nCells, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  const double mean_u2v2w2 = global_sum_u2v2w2 / static_cast<double>(global_nCells);
  const double urms_glob   = std::sqrt(mean_u2v2w2);        //Root mean square speed
  const double Mt_num      = urms_glob / initialSoundSpeed; // Turbulent Mach

  if (rankCpu == 0) {
    std::cout << "-------------------" << std::endl;
    std::cout << "Initial field check" << std::endl;
    std::cout << "  turbulent Mach         = " << Mt_num << std::endl;
    std::cout << "  conductivity           = " << initialConductivity << std::endl;
    std::cout << "  viscosity              = " << initialViscosity << std::endl;
    std::cout << "  root mean square speed = " << urms_glob / std::sqrt(3.) << std::endl;
    std::cout << "  sound speed            = " << initialSoundSpeed << std::endl;
    std::cout << "-------------------" << std::endl;
  }
}

//***********************************************************************
