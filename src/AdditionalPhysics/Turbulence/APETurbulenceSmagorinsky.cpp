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

#include "APETurbulenceSmagorinsky.h"

using namespace tinyxml2;

//***********************************************************************

APETurbulenceSmagorinsky::APETurbulenceSmagorinsky() {}

//***********************************************************************

APETurbulenceSmagorinsky::APETurbulenceSmagorinsky(int& numberQPA, Eos** eos, XMLElement* elementSmagorinsky, std::string fileName)
{
  XMLError error;
  //Cs
  error = elementSmagorinsky->QueryDoubleAttribute("C_s", &m_Cs);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("C_s", fileName, __FILE__, __LINE__);
  //CI
  error = elementSmagorinsky->QueryDoubleAttribute("C_I", &m_Ci);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("C_I", fileName, __FILE__, __LINE__);
  //Pr
  error = elementSmagorinsky->QueryDoubleAttribute("Pr_T", &m_PrT);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("Pr_T", fileName, __FILE__, __LINE__);
  //Ducros
  error = elementSmagorinsky->QueryBoolAttribute("ducros", &m_Ducros);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("m_Ducros", fileName, __FILE__, __LINE__);

  m_gamma  = eos[0]->getGamma();
  m_cv     = eos[0]->getCv();
  m_mu     = eos[0]->getMu();
  m_lambda = eos[0]->getLambda();
  m_numQPA = numberQPA++;
}

//***********************************************************************

APETurbulenceSmagorinsky::~APETurbulenceSmagorinsky() {}

//***********************************************************************

void APETurbulenceSmagorinsky::addQuantityAddPhys(Cell* cell) { cell->getVecQuantitiesAddPhys().push_back(new QAPTurbulence(this)); }

//***********************************************************************

void APETurbulenceSmagorinsky::solveFluxAddPhys(CellInterface* cellInterface)
{
  // Copy velocities, densities, cell volume and gradients of left and right cells
  //add volume calculation and pass as argument in the method
  m_velocityLeft  = cellInterface->getCellLeft()->getPhase(0)->getVelocity();
  m_velocityRight = cellInterface->getCellRight()->getPhase(0)->getVelocity();

  double rhoLeft  = cellInterface->getCellLeft()->getPhase(0)->getDensity();
  double rhoRight = cellInterface->getCellRight()->getPhase(0)->getDensity();

  double cellVolumeLeft = cellInterface->getCellLeft()->getElement()->getVolume();

  m_gradULeft  = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(1);
  m_gradURight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(1);
  m_gradVLeft  = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(2);
  m_gradVRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(2);
  m_gradWLeft  = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(3);
  m_gradWRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(3);

  m_gradTLeft  = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(4);
  m_gradTRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(4);

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

  m_gradTLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradTRight.localProjection(m_normal, m_tangent, m_binormal);

  this->solveFluxTurbulenceInner(cellInterface,
                                 m_velocityLeft,
                                 m_velocityRight,
                                 m_gradULeft,
                                 m_gradURight,
                                 m_gradVLeft,
                                 m_gradVRight,
                                 m_gradWLeft,
                                 m_gradWRight,
                                 m_gradTLeft,
                                 m_gradTRight,
                                 m_gamma,
                                 m_cv,
                                 m_mu,
                                 m_lambda,
                                 rhoLeft,
                                 rhoRight,
                                 cellVolumeLeft,
                                 m_tensorLeft,
                                 m_tensorRight);

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APETurbulenceSmagorinsky::solveFluxAddPhysBoundary(CellInterface* cellInterface)
{
  // Copy cell volume, density, gradients and velocity of left and right cells
  double cellVolumeLeft = cellInterface->getCellLeft()->getElement()->getVolume();

  double rhoRight(0.);
  double rhoLeft = cellInterface->getCellLeft()->getPhase(0)->getDensity();

  m_gradULeft = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(1);
  m_gradVLeft = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(2);
  m_gradWLeft = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(3);

  m_gradTLeft = cellInterface->getCellLeft()->getQPA(m_numQPA)->getGrad(4);

  m_velocityLeft = cellInterface->getCellLeft()->getPhase(0)->getVelocity();

  // Projection on orientation axes attached to the edge of velocities and gradients
  m_normal   = cellInterface->getFace()->getNormal();
  m_tangent  = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);

  m_tensorLeft.setTensorByLines(m_gradULeft, m_gradVLeft, m_gradWLeft);
  m_tensorLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_tensorLeft.tensorToCoords(m_gradULeft, m_gradVLeft, m_gradWLeft);

  m_gradTLeft.localProjection(m_normal, m_tangent, m_binormal);

  // Distances cells/cell interfaces for weighting on the flux
  // double distLeft = cellInterface->getCellLeft()->distance(cellInterface); //ET// for wall BCs
  switch (cellInterface->whoAmI()) {
  case TypeBC::NONREFLECTING:
  case TypeBC::OUTLETPRESSURE:
  case TypeBC::INLETINJSTAGSTATE:
  case TypeBC::INLETTANK:
  case TypeBC::INLETINJTEMP:
    this->solveFluxTurbulenceNonReflecting(
      cellInterface, m_velocityLeft, m_gradULeft, m_gradVLeft, m_gradWLeft, m_gradTLeft, rhoLeft, cellVolumeLeft, m_tensorLeft);
    break;
  // case TypeBC::WALL:
  //   if (!cellInterface->isMRFWall()) {
  //     this->solveFluxTurbulenceWall(m_velocityLeft, rhoLeft, distLeft);
  //   }
  //   else {
  //     Coord velocityWall(0.);
  //     velocityWall = cellInterface->getWallRotationalVelocityMRF().cross(cellInterface->getFace()->getPos());
  //     velocityWall.localProjection(m_normal, m_tangent, m_binormal);
  //     velocityWall.setX(0.); // Fluid not crossing boundary
  //     this->solveFluxTurbulenceWall(m_velocityLeft, rhoLeft, distLeft, velocityWall);
  //   }
  //   break;
  case TypeBC::SYMMETRY:
    this->solveFluxTurbulenceSymmetry(cellInterface, m_gradULeft, m_gradVLeft, m_gradWLeft, rhoLeft, cellVolumeLeft, m_tensorLeft, m_mu);
    break;
  case TypeBC::PERIODIC:
    rhoRight = cellInterface->getCellRight()->getPhase(0)->getDensity();

    m_velocityRight = cellInterface->getCellRight()->getPhase(0)->getVelocity();
    m_velocityRight.localProjection(m_normal, m_tangent, m_binormal);

    m_gradURight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(1);
    m_gradVRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(2);
    m_gradWRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(3);

    m_gradTRight = cellInterface->getCellRight()->getQPA(m_numQPA)->getGrad(4);

    m_tensorRight.setTensorByLines(m_gradURight, m_gradVRight, m_gradWRight);
    m_tensorRight.localProjection(m_normal, m_tangent, m_binormal);
    m_tensorRight.tensorToCoords(m_gradURight, m_gradVRight, m_gradWRight);

    m_gradTRight.localProjection(m_normal, m_tangent, m_binormal);

    this->solveFluxTurbulenceInner(cellInterface,
                                   m_velocityLeft,
                                   m_velocityRight,
                                   m_gradULeft,
                                   m_gradURight,
                                   m_gradVLeft,
                                   m_gradVRight,
                                   m_gradWLeft,
                                   m_gradWRight,
                                   m_gradTLeft,
                                   m_gradTRight,
                                   m_gamma,
                                   m_cv,
                                   m_mu,
                                   m_lambda,
                                   rhoLeft,
                                   rhoRight,
                                   cellVolumeLeft,
                                   m_tensorLeft,
                                   m_tensorRight);
    break;
  default:
    this->solveFluxTurbulenceOther();
  }

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APETurbulenceSmagorinsky::solveFluxTurbulenceInner(CellInterface* cellInterface,
                                                        const Coord& velocityLeft,
                                                        const Coord& velocityRight,
                                                        const Coord& gradULeft,
                                                        const Coord& gradURight,
                                                        const Coord& gradVLeft,
                                                        const Coord& gradVRight,
                                                        const Coord& gradWLeft,
                                                        const Coord& gradWRight,
                                                        const Coord& gradTLeft,
                                                        const Coord& gradTRight,
                                                        const double& gamma,
                                                        const double& cv,
                                                        const double& mu,
                                                        const double& lambda,
                                                        const double& rhoLeft,
                                                        const double& rhoRight,
                                                        const double& cellVolumeLeft,
                                                        const Tensor& strainRateLeft,
                                                        const Tensor& strainRateRight) const
{
  //Extraction of data
  double uL(0.), vL(0.), wL(0.), uR(0.), vR(0.), wR(0.);
  double dudxL(0.), dudyL(0.), dudzL(0.), dudxR(0.), dudyR(0.), dudzR(0.);
  double dvdxL(0.), dvdyL(0.), dvdxR(0.), dvdyR(0.);
  double dwdxL(0.), dwdzL(0.), dwdxR(0.), dwdzR(0.);
  double dTdxL(0.), dTdxR(0.);
  double cbrDelta(0.);                // Cell size (cube root delta)
  double magSLeft(0.), magSRight(0.); // Magnitude stress tensor S left or right
  double magSgsS(0.);                 // Magnitude of sub-grid stress tensor S

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

  dTdxL = gradTLeft.getX();
  dTdxR = gradTRight.getX();

  //Data of the cell interface
  double u(0.), v(0.), w(0.);
  double rho(0.);
  double dudx(0.), dudy(0.), dudz(0.);
  double dvdx(0.), dvdy(0.);
  double dwdx(0.), dwdz(0.);
  double dTdx(0.);

  u = (uL + uR) * 0.5;
  v = (vL + vR) * 0.5;
  w = (wL + wR) * 0.5;

  rho = (rhoLeft + rhoRight) * 0.5;

  dudx = (dudxL + dudxR) * 0.5;
  dudy = (dudyL + dudyR) * 0.5;
  dudz = (dudzL + dudzR) * 0.5;

  dvdx = (dvdxL + dvdxR) * 0.5;
  dvdy = (dvdyL + dvdyR) * 0.5;

  dwdx = (dwdxL + dwdxR) * 0.5;
  dwdz = (dwdzL + dwdzR) * 0.5;

  dTdx = (dTdxL + dTdxR) * 0.5;

  //Model variables
  tensorStrainRate  = strainRateLeft;
  tensorStrainRate += strainRateRight;
  tensorStrainRate *= 0.5;
  tensorStrainRate.transpose(tensorStrainRateTranspose);
  tensorStrainRate += tensorStrainRateTranspose;
  tensorStrainRate *= 0.5; //symmetric part of the gradU
  magSgsS           = std::sqrt(2.) * tensorStrainRate.frobeniusNorm();
  cbrDelta          = std::cbrt(cellVolumeLeft);

  double Cs2delta2 = (m_Cs * m_Cs) * (cbrDelta * cbrDelta);
  double Ci2Delta2 = m_Ci * 2. * (cbrDelta * cbrDelta);

  double mu_T     = Cs2delta2 * rho * magSgsS;
  double lambda_T = gamma * cv * mu_T / m_PrT;

  //Preparing quantities for the output of turbulent viscosity field
  const double sxxL = strainRateLeft.getXX();
  const double syyL = strainRateLeft.getYY();
  const double szzL = strainRateLeft.getZZ();
  const double sxyL = 0.5 * (strainRateLeft.getXY() + strainRateLeft.getYX());
  const double sxzL = 0.5 * (strainRateLeft.getXZ() + strainRateLeft.getZX());
  const double syzL = 0.5 * (strainRateLeft.getYZ() + strainRateLeft.getZY());
  magSLeft          = std::sqrt(2. * (sxxL * sxxL + syyL * syyL + szzL * szzL + 2. * (sxyL * sxyL + sxzL * sxzL + syzL * syzL)));

  const double sxxR = strainRateRight.getXX();
  const double syyR = strainRateRight.getYY();
  const double szzR = strainRateRight.getZZ();
  const double sxyR = 0.5 * (strainRateRight.getXY() + strainRateRight.getYX());
  const double sxzR = 0.5 * (strainRateRight.getXZ() + strainRateRight.getZX());
  const double syzR = 0.5 * (strainRateRight.getYZ() + strainRateRight.getZY());
  magSRight         = std::sqrt(2. * (sxxR * sxxR + syyR * syyR + szzR * szzR + 2. * (sxyR * sxyR + sxzR * sxzR + syzR * syzR)));

  const double Skk = dudx + dvdy + dwdz;
  const double Sxy = 0.5 * (dudy + dvdx);
  const double Sxz = 0.5 * (dudz + dwdx);

  double taukk = Ci2Delta2 * rho * (magSgsS * magSgsS);

  //Ducros sensor to turn off the SGS stresses in the zone with strong shocks
  double switch_LES(1.);
  if (m_Ducros) {
    double dvdzL(0.), dvdzR(0.), dwdyL(0.), dwdyR(0.);
    dvdzL = gradVLeft.getZ();
    dvdzR = gradVRight.getZ();
    dwdyL = gradWLeft.getY();
    dwdyR = gradWRight.getY();

    double dvdz(0.), dwdy(0.);
    dvdz = (dvdzL + dvdzR) * 0.5;
    dwdy = (dwdyL + dwdyR) * 0.5;

    double w1      = dwdy - dvdz;
    double w2      = dudz - dwdx;
    double w3      = dvdx - dudy;
    double omegasq = (w1 * w1) + (w2 * w2) + (w3 * w3);

    double sensor_ducros = (Skk * Skk) / ((Skk * Skk) + omegasq + 1e-30); // Skk is equal to theta in the nomenclature of the sensor

    // double switch_LES = std::max(0.0, 1.0 - sensor_ducros); //linear option
    switch_LES = std::pow(1.0 - sensor_ducros, 9); //sharp cutoff switch
  }

  mu_T     *= switch_LES;
  taukk    *= switch_LES;
  lambda_T *= switch_LES;

  // Output of nu_T
  double nuTLeft  = Cs2delta2 * magSLeft * switch_LES;
  double nuTRight = Cs2delta2 * magSRight * switch_LES;

  if (cellInterface->getCellLeft()) {
    auto qpaL = cellInterface->getCellLeft()->getQPA(m_numQPA);
    if (qpaL) qpaL->setTurbulentViscosity(nuTLeft);
  }
  if (cellInterface->getCellRight()) {
    auto qpaR = cellInterface->getCellRight()->getQPA(m_numQPA);
    if (qpaR) qpaR->setTurbulentViscosity(nuTRight);
  }

  const double tau_xx_phys = -2. * mu * (dudx - Skk / 3.);
  const double tau_xy_phys = -2. * mu * Sxy;
  const double tau_xz_phys = -2. * mu * Sxz;

  const double tau_xx_sgs = -2. * mu_T * (dudx - Skk / 3.) + taukk / 3.;
  const double tau_xy_sgs = -2. * mu_T * Sxy;
  const double tau_xz_sgs = -2. * mu_T * Sxz;

  double Jx = 0.5 * (u * tau_xx_sgs + v * tau_xy_sgs + w * tau_xz_sgs);

  lambda_T += lambda; //Physical conductivity flux

  double Qx = (-lambda_T * dTdx);

  //Writing of SGS terms on each equation of fluxBuffEuler
  static_cast<FluxEuler*>(fluxBuff)->m_mass = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setX(tau_xx_phys + tau_xx_sgs);
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setY(tau_xy_phys + tau_xy_sgs);
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setZ(tau_xz_phys + tau_xz_sgs);
  static_cast<FluxEuler*>(fluxBuff)->m_energ = Qx + Jx;
}

//***********************************************************************

void APETurbulenceSmagorinsky::solveFluxTurbulenceNonReflecting(CellInterface* cellInterface,
                                                                const Coord& velocityLeft,
                                                                const Coord& gradULeft,
                                                                const Coord& gradVLeft,
                                                                const Coord& gradWLeft,
                                                                const Coord& gradTLeft,
                                                                const double& rhoLeft,
                                                                const double& cellVolumeLeft,
                                                                const Tensor& strainRateLeft) const
{
  this->solveFluxTurbulenceInner(cellInterface,
                                 velocityLeft,
                                 velocityLeft,
                                 gradULeft,
                                 gradULeft,
                                 gradVLeft,
                                 gradVLeft,
                                 gradWLeft,
                                 gradWLeft,
                                 gradTLeft,
                                 gradTLeft,
                                 m_gamma,
                                 m_cv,
                                 m_mu,
                                 m_lambda,
                                 rhoLeft,
                                 rhoLeft,
                                 cellVolumeLeft,
                                 strainRateLeft,
                                 strainRateLeft);
}

//***********************************************************************

// void APETurbulenceSmagorinsky::solveFluxTurbulenceWall(const Coord& velocityLeft, const double& rhoLeft, const double& distLeft, Coord& velocityWall) const
// {
//   //ET//Wall treatment to come
//   static_cast<FluxEuler*>(fluxBuff)->m_mass     = 0.;
//   static_cast<FluxEuler*>(fluxBuff)->m_momentum = 0.;
//   static_cast<FluxEuler*>(fluxBuff)->m_energ    = 0.;
// }

//***********************************************************************

void APETurbulenceSmagorinsky::solveFluxTurbulenceSymmetry(CellInterface* cellInterface,
                                                           const Coord& gradULeft,
                                                           const Coord& gradVLeft,
                                                           const Coord& gradWLeft,
                                                           const double& rhoLeft,
                                                           const double& cellVolumeLeft,
                                                           const Tensor& strainRateLeft,
                                                           const double& mu) const
{
  //Extraction of data
  //Note that dudy, dudz, dvdx and dwdx are nulls
  double dudx(0.), dvdy(0.), dwdz(0.);
  double cbrDelta(0.); // Cell size (cube root delta)
  double magSgsS(0.);  // Magnitude of sub-grid stress tensor S
  double magSLeft(0.); // Magnitude stress tensor S left
  Tensor strainRateLeftTransposed(0., 0., 0., 0., 0., 0., 0., 0., 0.), symLeftS(0., 0., 0., 0., 0., 0., 0., 0., 0.);
  double rho = rhoLeft;

  dudx = gradULeft.getX();
  dvdy = gradVLeft.getY();
  dwdz = gradWLeft.getZ();

  //Model variables
  tensorStrainRate = strainRateLeft;
  tensorStrainRate.transpose(tensorStrainRateTranspose);
  tensorStrainRate += tensorStrainRateTranspose;
  tensorStrainRate *= 0.5; //symmetric part of the gradU
  magSgsS           = std::sqrt(2.) * tensorStrainRate.frobeniusNorm();
  cbrDelta          = std::cbrt(cellVolumeLeft);

  //Output turbulent viscosity field
  strainRateLeft.transpose(strainRateLeftTransposed);
  symLeftS        = strainRateLeft;
  symLeftS       += strainRateLeftTransposed;
  symLeftS       *= 0.5;
  magSLeft        = std::sqrt(2.) * symLeftS.frobeniusNorm();
  double nuTLeft  = (m_Cs * m_Cs) * (cbrDelta * cbrDelta) * magSLeft;

  cellInterface->getCellLeft()->getQPA(m_numQPA)->setTurbulentViscosity(nuTLeft);

  double mu_T = (m_Cs * m_Cs) * (cbrDelta * cbrDelta) * rho * magSgsS;

  const double Skk = dudx + dvdy + dwdz;

  double taukk        = m_Ci * 2. * rho * (cbrDelta * cbrDelta) * (magSgsS * magSgsS);
  const double tau_xx = -2. * mu_T * (dudx - Skk / 3.) + taukk / 3.;

  //Ducros sensor to turn off the SGS stresses in the zone with strong shocks
  double switch_LES(1.);
  if (m_Ducros) {
    double dwdy(0.), dvdz(0.);
    dwdy = gradWLeft.getY();
    dvdz = gradVLeft.getZ();

    double theta = dudx + dvdy + dwdz;

    double w1      = dwdy - dvdz;
    double omegasq = (w1 * w1);

    double sensor_ducros = (theta * theta) / ((theta * theta) + (omegasq) + 1e-30);

    // switch_LES = std::max(0.0, 1.0 - sensor_ducros); //linear option
    switch_LES = std::pow(1.0 - sensor_ducros, 9); //sharp cutoff switch
  }

  mu_T  *= switch_LES;
  taukk *= switch_LES;

  mu_T += mu;

  //Writing of SGS terms on each equation of fluxBuffEuler
  static_cast<FluxEuler*>(fluxBuff)->m_mass = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setX(tau_xx);
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxEuler*>(fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxEuler*>(fluxBuff)->m_energ = 0.;
}

//***********************************************************************

void APETurbulenceSmagorinsky::solveFluxTurbulenceOther() const
{
  //Not manage at the moment, just an example
  //std::cout << "Turbulent boundary not managed" << std::endl;

  // To avoid bug when not managed
  static_cast<FluxEuler*>(fluxBuff)->m_mass     = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_momentum = 0.;
  static_cast<FluxEuler*>(fluxBuff)->m_energ    = 0.;
}

//***********************************************************************

void APETurbulenceSmagorinsky::communicationsAddPhys(const int& dim, const int& lvl)
{
  parallel.communicationsVector(Variable::QPA, dim, lvl, m_numQPA, 1); //m_gradU
  parallel.communicationsVector(Variable::QPA, dim, lvl, m_numQPA, 2); //m_gradV
  parallel.communicationsVector(Variable::QPA, dim, lvl, m_numQPA, 3); //m_gradW
  parallel.communicationsVector(Variable::QPA, dim, lvl, m_numQPA, 4); //m_gradT
}

//***********************************************************************

void APETurbulenceSmagorinsky::computeTurbulenceParameters() const
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
