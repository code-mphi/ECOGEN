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

#include "GDEntireDomainWithParticularities.h"
#include <random>

bool algorithmSaadOn = false;

//***********************************************************

GDEntireDomainWithParticularities::GDEntireDomainWithParticularities(
  std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, const int& physicalEntity) :
  GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{}

//***********************************************************

GDEntireDomainWithParticularities::~GDEntireDomainWithParticularities() {}

//***********************************************************

bool GDEntireDomainWithParticularities::belong(Coord& /*posElement*/, const int& /*lvl*/) const
{
  //1. Laplace pressure initialization
  //----------------------------------
  return true; //always belong to entire domain

  //2. Respecting special coordinates
  //---------------------------------
  //bool result(false);
  //if (posElement.getY() - posElement.getX() >= -2.-1.e-8) { result = true; }
  //return result;

  //3. Respecting special coordinates (with AMR test)
  //-------------------------------------------------
  // bool result(false);
  // if (lvl > 0) {
  //   if (posElement.getX() < 0.02 / std::pow(2., (double)(lvl))) { result = true; }
  // }
  // else {
  //   if (posElement.getX() < 0.02) { result = true; }
  // }
  // return result;

  //4. Random velocity perturbations
  //--------------------------------
  // return true; //always belong to entire domain

  //5. Rayleigh-Taylor instability
  //------------------------------
  // return true; //always belong to entire domain

  //6. Blast-wave equation
  //----------------------
  // return true; //always belong to entire domain

  //7. Ultrasound equation
  //----------------------
  // return true; //always belong to entire domain

  //8. Exponential decay for wave equation
  //--------------------------------------
  // return true; //always belong to entire domain

  //9. Sine wave viscosity decay test
  //---------------------------------
  // return true; //always belong to entire domain

  //10. 3D isotropic compressible turbulence decay
  //----------------------------------------------
  // return true; //always belong to entire domain
}

//******************************************************************

void GDEntireDomainWithParticularities::fillIn(Cell* cell) const
{
  //As basic fillIn: Test if the cell belongs to the geometrical domain
  bool belongs(true);
  if (cell->getElement() != 0) {
    Coord coord(cell->getPosition());
    if (!this->belong(coord, cell->getLvl())) {
      belongs = false;
    }
    //Test if the cell belongs to physical mesh entity (for unstructured meshes)
    if (cell->getElement()->getAppartenancePhysique() > 0 && m_physicalEntity > 0) {
      if (cell->getElement()->getAppartenancePhysique() != m_physicalEntity) {
        belongs = false;
      }
    }
  }

  if (belongs) {
    for (int k = 0; k < numberPhases; k++) {
      cell->copyPhase(k, m_vecPhases[k]);
    }
    cell->copyMixture(m_mixture);
    for (int k = 0; k < numberTransports; k++) {
      cell->setTransport(m_vecTransports[k].getValue(), k);
    }
    if (m_physicalEntity == -1) {
      cell->setWall(true);
    }
    else {
      cell->setWall(false);
    }

    //Particularities
    //1. Laplace pressure initialization
    //----------------------------------
    // if (cell->getElement() != 0) {
    //   double pressure(0.);
    //   Coord posElement(cell->getPosition());
    //   double radius;
    //   //radius = posElement.getX(); //1D
    //   radius = std::pow(std::pow(posElement.getX(), 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   //radius = std::pow(std::pow(posElement.getX() - 0.75e-3, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   // radius = std::pow(std::pow(posElement.getX() - 2.e-4, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   //radius = std::pow(std::pow(posElement.getX(), 2.) + std::pow(posElement.getY(), 2.) + std::pow(posElement.getZ(), 2.), 0.5); //3D
    //   //radius = std::pow(std::pow(posElement.getX() - 153.6e-3, 2.) + std::pow(posElement.getY(), 2.) + std::pow(posElement.getZ(), 2.), 0.5); //3D
    //   //pressure = 1.e5 + 1.e-3 / radius * (1.e4 - 1.e5);
    //   //pressure = 1.e5 + 1.e-3 / radius * (4.e3 - 1.e5);
    //   //pressure = 1.e5 + 1.e-3 / radius * (1.e3 - 1.e5);
    //   //pressure = 353.e5 + 1.e-4 / radius * (1.e5 - 353.e5);
    //   //pressure = 50.6625e5 + 1.e-4 / radius * (3.55e3 - 50.6625e5);
    //   // pressure = 1.01325e5 + 2.72e-3 / radius * (2700. - 1.01325e5);
    //   pressure = 1.01325e5 + 2.72e-3 / radius * (2700. - 1.01325e5);
    //   double radius2 = std::pow(std::pow(posElement.getX() - 3.7e-3, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.51e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;
    //   for (int k = 0; k < numberPhases; k++) { cell->getPhase(k)->setPressure(pressure); }
    //   cell->getMixture()->setPressure(pressure);
    // }

    //  BEGO - 2 SIMULATION
    //  Feb-Mar 2025
    //  ----------------------------------------------------------------------------------------------------------------------------------------------

    //  1.b Laplace pressure initialization
    //  -----------------------------------
    // if (cell->getElement() != 0) {
    //   double pressure(0.);
    //   Coord posElement(cell->getPosition());
    //   double radius;

    //   CASE WITH NO DISTANCE BUBBLE
    //   case 1 - no small bubble, offset big bubble = 0.0004
    //   radius = std::pow(std::pow(posElement.getX() - 0.2009, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //
    //   case 2 - no small bubble, offset big bubble = 0
    //   radius = std::pow(std::pow(posElement.getX() - 0.2005, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D

    //   CASE WITH DISTANCE BUBBLE
    //   case 3 - small bubble, offset big bubble = 0.0004
    //   radius = std::pow(std::pow(posElement.getX() - 0.2009, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D

    //   Add effect of the big bubble
    //   pressure = 1.01325e5 + 2.7e-3 / radius * (2700. - 1.01325e5);

    //   SENSITIVITY ON DISTANCE BUBBLE STAND OFF
    //   stand-off 1.10
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.20415, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 1.50
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.20435, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 1.375
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.2042875, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 1.425
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.2043125, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 1.75
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.204475, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 2.125
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.2046625, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 2.00
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.2046, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 2.50
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.20485, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 3.75
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.205475, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    //   stand-off 5.00
    //   double radius2 = std::pow(std::pow(posElement.getX() - 0.2061, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    //   pressure = (pressure + pressure2) / 2.;

    // ---------------------- MATCH EXPERIMENT (BASE CASES) ---------------------------------------------------------
    // Small cavity - Ellipse parameters (Case Experiment standoff 1.375, with r1 = 0.4)
    // stand-off 1.375
    // double radius2 = std::pow(std::pow(posElement.getX() - 0.20415, 2.) + std::pow(posElement.getY(), 2.), 0.5);     //2D
    // double pressure2 = 1.01325e5 + 0.4e-3 / radius2 * (2700. - 1.01325e5);
    // pressure = (pressure + pressure2) / 2.;

    // Small cavity - Ellipse parameters (Case Experiment standoff 1.375, with r1 = 0.6)
    // stand-off 1.375
    // double radius2 = std::pow(std::pow(posElement.getX() - 0.204425, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    // double pressure2 = 1.01325e5 + 0.6e-3 / radius2 * (2700. - 1.01325e5);
    // pressure = (pressure + pressure2) / 2.;

    // ----------------------ASPECT RATIO (BASE CASES) ------------------------------------------------------------------
    // Small cavity - Ellipse parameters (Case Experiment with standoff 1.375, with r1 = 0.5. r2 =0.5, aspect ratio = 1)
    // double radius2 = std::pow(std::pow(posElement.getX() - 0.2042875, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    // double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    // pressure = (pressure + pressure2) / 2.;

    // Small cavity - Ellipse parameters (Case Experiment with standoff 1.375, with r1 = 0.5. r2 = 0.575, aspect ratio = 1.15)
    // double radius2 = std::pow(std::pow(posElement.getX() - 0.2042875, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    // double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    // pressure = (pressure + pressure2) / 2.;

    // Small cavity - Ellipse parameters (Case Experiment with standoff 1.375, with r1 = 0.5. r2 = 0.725, aspect ratio = 1.45)
    // double radius2 = std::pow(std::pow(posElement.getX() - 0.2042875, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    // double pressure2 = 1.01325e5 + 0.5e-3 / radius2 * (2700. - 1.01325e5);
    // pressure = (pressure + pressure2) / 2.;

    // --------------------------------------------------------------------------------------------------------------

    //    for (int k = 0; k < numberPhases; k++) { cell->getPhase(k)->setPressure(pressure); }
    //    cell->getMixture()->setPressure(pressure);
    // }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    // END BEGO - 2 SIMULATION

    //2. Respecting special coordinates
    //---------------------------------
    //Nothing special here

    //3. Respecting special coordinates (with AMR test)
    //-------------------------------------------------
    //if (cell->getElement() != 0) {
    //  double pressure(0.);
    //  Coord posElement(cell->getPosition());
    //  if (posElement.getX() < 0.025 && posElement.getY() < 1.265) {
    //    //pressure = 1e1 + (1.55e6 - 1e1) / 1.265 * posElement.getY();
    //    pressure = 1e1;
    //    for (int k = 0; k < numberPhases; k++) { cell->getPhase(k)->setPressure(pressure); }
    //    cell->getMixture()->setPressure(pressure);
    //  }
    //}

    //4. Random velocity perturbations: O(1e−4 u_s)
    //---------------------------------------------
    // if (cell->getElement() != 0) {
    //   Coord perturbedVelocity(cell->getMixture()->getVelocity());
    //   perturbedVelocity.setX(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getX());
    //   perturbedVelocity.setY(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getY());
    //   perturbedVelocity.setZ(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getZ());
    //   cell->getMixture()->setVelocity(perturbedVelocity);
    // }

    //5. Rayleigh-Taylor instability
    //------------------------------
    // Sinus shape parameters
    // const double pi = std::acos(-1); // Pi constant

    // // * R-T same kinematic viscosity
    // double lambda = 0.2;             // Width of the domain
    // double h = 0.7;                  // Height of the interface

    // // * R-T same dynamic viscosity
    // double lambda = 2.5;             // Width of the domain
    // double h = 12.5;                 // Height of the interface

    // double k = 2 * pi / lambda;      // Wave-number
    // int nx = 1000;                   // Nb of points to plot interface function
    // double dx = lambda / (nx - 1.);  // Points spacing for the plot of the interface fn
    // double amp = 0.05 / k;           // Amplitude of the sinus function

    // double rhoHeavy = cell->getPhase(0)->getDensity(),
    //   rhoLight = cell->getPhase(1)->getDensity();
    // std::vector<double> alpha(2, 0.);

    // // Hydrostatic pressure
    // double pref = 1.e5, pinterface = pref, pressure = 0.;

    // // * R-T same kinematic viscosity
    // double g = 9.81, ly = 1.2;

    // // * R-T same dynamic viscosity
    // double g = 1., ly = 25.;

    // std::vector<double> interfaceX, interfaceY; // Interface fn coordinates
    // for (int i = 0; i < nx; i++) {
    //   interfaceX.push_back(i * dx);
    //   interfaceY.push_back(amp * std::sin(2. * pi * interfaceX[i] / lambda + pi / 2.) + h);
    // }

    // int index = 0;

    // if (cell->getElement() != 0) {
    //   Coord posElement(cell->getPosition());
    //   double minVal = 1.;
    //   for (unsigned int i = 0; i < interfaceX.size(); i++) { // Find nearest index corresponding to x-position of interface fn
    //     if (std::abs(posElement.getX() - interfaceX[i]) < minVal) {
    //       minVal = std::abs(posElement.getX() - interfaceX[i]);
    //       index = i;
    //     }
    //   }

    //   // Check location to interface and initialize domain accordingly
    //   if (posElement.getY() > interfaceY[index]) {
    //   alpha[0] = 1.;
    //   alpha[1] = 0.;
    //   pressure = pref + rhoHeavy * g * (ly - posElement.getY());
    //   }
    //   else {
    //     alpha[0] = 0.;
    //     alpha[1] = 1.;
    //     pinterface = pref + rhoHeavy * g * (ly - interfaceY[index]);
    //     pressure = pinterface + rhoLight * g * (interfaceY[index] - posElement.getY());
    //   }

    //   for (int k = 0; k < numberPhases; k++) {
    //     cell->getPhase(k)->setAlpha(alpha[k]);
    //     cell->getPhase(k)->setPressure(pressure);
    //   }
    //   cell->getMixture()->setPressure(pressure);
    // }

    //6. Blast-wave equation
    //----------------------
    //p(t) = p0 + 2 p_s exp(−αt) * cos(ωt + π/3)
    //p(x) = p0 + 2 p_s exp(-αr/c) * cos(ωr/c + π/3)
    // if (cell->getElement() != 0) {
    //   double pressure(0.), velocity(0.), pk(0.);
    //   double beta(1.48e6), omega(1.21e6); //beta here is the alpha variable of the equation
    //   // double beta(10), omega(0.0009);
    //   double posX(cell->getPosition().getX());
    //   double p0(1.01325e5), pS(1.e6), soundSpeed(1625.), density(998.);
    //   // double p0(1.01325e5), pS(94.4e6), soundSpeed(1534.6), density(998.);
    //   double shockFront(2.e-3);

    //   //Right-facing shock
    //   if (posX < shockFront) {
    //     double r = posX - shockFront;
    //     pressure = p0 + 2. * pS * exp(beta * r / soundSpeed) * std::cos(- omega * r / soundSpeed + M_PI / 3.);
    //     // pressure = p0 + pS * exp(beta * r / omega) * (1. - (r / omega));
    //     velocity = (pressure - p0) / (density * soundSpeed);

    //     for (int k = 0; k < numberPhases; k++) {
    //       pk = pressure;
    //       cell->getPhase(k)->getEos()->verifyAndModifyPressure(pk);
    //       cell->getPhase(k)->setPressure(pk);
    //     }
    //     cell->getMixture()->setPressure(pressure);
    //     cell->getMixture()->setU(velocity);
    //   }

    //   // //Left-facing shock
    //   // if (posX > shockFront) {
    //   //   double r = shockFront - posX;
    //   //   pressure = p0 + 2. * pS * exp(beta * r / soundSpeed) * std::cos(- omega * r / soundSpeed + M_PI / 3.);
    //   //   velocity = - (pressure - p0) / (density * soundSpeed);

    //   //   for (int k = 0; k < numberPhases; k++) {
    //   //     pk = pressure;
    //   //     cell->getPhase(k)->getEos()->verifyAndModifyPressure(pk);
    //   //     cell->getPhase(k)->setPressure(pk);
    //   //   }
    //   //   cell->getMixture()->setPressure(pressure);
    //   //   cell->getMixture()->setU(velocity);
    //   // }
    // }

    //7. Ultrasound equation
    //----------------------
    // if (cell->getElement() != 0) {
    //   double pressure(0.), velocity(0.), pk(0.);
    //   double omega(2. * M_PI * 5.e6);
    //   double posX(cell->getPosition().getX());
    //   double p0(1.01325e5), pS(2.e6), soundSpeed(1500.), density(1000.);
    //   double ultrasoundFront(300.e-6);

    //   //Right-facing ultrasound
    //   if (posX < ultrasoundFront) {
    //     double r = posX - ultrasoundFront;
    //     pressure = p0 + pS * std::sin(- omega * r / soundSpeed);
    //     velocity = (pressure - p0) / (density * soundSpeed);

    //     for (int k = 0; k < numberPhases; k++) {
    //       pk = pressure;
    //       cell->getPhase(k)->getEos()->verifyAndModifyPressure(pk);
    //       cell->getPhase(k)->setPressure(pk);
    //     }
    //     cell->getMixture()->setPressure(pressure);
    //     cell->getMixture()->setU(velocity);
    //   }
    // }

    //8. Exponential decay for wave equation
    //--------------------------------------
    //p(x,t) = f(x) * g(x + ct) + f(x) * g(x - ct)
    //f(x) = 1 if cartesian; f(x) = 1 / x if spherical
    //g(x +- ct) = ( x +- ct ) / 2 * p0 * exp( -pi * ( x +- ct ) ^ 2 ) -- expected analytic solution
    //g(x) = x / 2 * p0 * exp( -pi * x ^ 2 ) -- implemented initial condition
    // if (cell->getElement() != 0) {
    //   double pressure(0.), pk(0.);
    //   double posX(cell->getPosition().getX());
    //   double p0(101325);

    //   pressure = p0 * (1. + 1.e-4 * exp(-M_PI * pow(posX, 2.)));

    //   for (int k = 0; k < numberPhases; k++) {
    //     pk = pressure;
    //     cell->getPhase(k)->getEos()->verifyAndModifyPressure(pk);
    //     cell->getPhase(k)->setPressure(pk);
    //   }
    //   cell->getMixture()->setPressure(pressure);
    // }

    //9. Sine wave viscosity decay test
    //---------------------------------
    //v(x,t) = v0 * exp(- nu * k ^ 2 * t) * sin(k * x) -- expected analytic solution
    //v(x,t) = v0 * sin(k * x) -- implemented initial condition
    // if (cell->getElement() != 0) {
    //   const double pi    = std::acos(-1.);
    //   const double Lx    = 1.;
    //   const int    mode  = 1;
    //   const double kappax    = 2. * pi * mode / Lx;
    //   const double V0    = 18.8;

    //   const double x     = cell->getPosition().getX();
    //   const double v_val = V0 * std::sin(kappax * x);

    //   for (int i = 0; i < numberPhases; ++i) {
    //     cell->getPhase(i)->setV(v_val);
    //   }
    //   cell->getMixture()->setV(v_val);
    // }

    //10. 3D isotropic compressible turbulence decay
    //----------------------------------------------
    //generate the synthetic initial conditions for a compressible homogeneous isotropic
    //turbulence decay simulation. The script generates a solenoidal velocity field that
    //is constructed to have the following energy E(k) ~ k^4 * exp(-2 * (k/k0)^2)
    // if (cell->getElement() != 0) {
    //   algorithmSaadOn = true;

    //   // constants
    //   const double pi     = 3.1415926535897932;
    //   const double gamma  = 1.4;
    //   const double cp     = 1173;         // J/kg/K
    //   const double Pr     = 0.71;
    //   const double T0     = 1200.;        // K
    //   const double p0     = 1.01300e5;    // Pa
    //   const double k0     = 4.;           // most energetic wavenumber
    //   const double R      = 335.143;
    //   const double rho0   = p0 / (R * T0);

    //   // prescribed quantities
    //   double Mt(0.6), Re_psi(100.);

    //   // Step 1 - thermodynamic parameters
    //   double initialRootMeanSquareSpeed(0.), psi0(0.);
    //   psi0    = 2. / k0;
    //   initialSoundSpeed = std::sqrt(gamma * R * T0);
    //   initialRootMeanSquareSpeed   = (Mt * initialSoundSpeed) / std::sqrt(3.);
    //   initialViscosity     = (rho0 * psi0 * initialRootMeanSquareSpeed) / Re_psi;
    //   initialConductivity = initialViscosity * cp / Pr;
    //   //check that the mu and lambda parameters in the IG_cHIT_decay.xml are the same
    //   //as the ones calculated here (debug will show in console)

    //   // Step 2 - target energy spectrum
    //   // E(k) = a * k^4 * exp(-2*(k/k0)^2)
    //   double a = (32. / std::sqrt(2. * pi)) * (initialRootMeanSquareSpeed * initialRootMeanSquareSpeed) / std::pow(k0, 5.);

    //   const int    M  = 8000;     // number of modes; about 8000 for 64^3, up to 20000 for 512^3
    //   const double L  = 2. * pi;   // box length
    //   const int    Nx = 64;        // number of cells along direction x
    //   const int    Ny = 64;
    //   const int    Nz = 64;

    //   double Etot                   = 1.5 * initialRootMeanSquareSpeed * initialRootMeanSquareSpeed;
    //   static bool modes_initialized = false;
    //   static double kappa[M];
    //   static double q    [M];
    //   static double khat_x[M], khat_y[M], khat_z[M];
    //   static double sigma_x[M], sigma_y[M], sigma_z[M];
    //   static double phase[M];

    //   if (!modes_initialized) {

    //     const double dx = L / static_cast<double>(Nx);
    //     const double dy = L / static_cast<double>(Ny);
    //     const double dz = L / static_cast<double>(Nz);

    //     const double kmin = 2. * pi / L;     // fundamental mode
    //     const double kmax_x = pi  / dx;
    //     const double kmax_y = pi  / dy;
    //     const double kmax_z = pi  / dz;
    //     const double kmax   = std::min({kmax_x, kmax_y, kmax_z}); // Nyquist limit

    //     std::mt19937 rng(42); //seed for reproducibility, change if needed (fun fact: rng(42) is lucky for 64^3 with 8k modes ^^)
    //     std::uniform_real_distribution<double> U01(0., 1.);
    //     std::uniform_real_distribution<double> Um11(-1., 1.);
    //     std::uniform_real_distribution<double> Uphase(-0.5 * pi, 0.5 * pi);

    //     const double log_ratio = std::log(kmax / kmin);

    //     // precompute shell edges (M+1 values)
    //     // e[m] = kmin * (kmax/kmin)^(m/M)
    //     auto edge = [&](int m) {
    //       return kmin * std::exp((double)m / (double)M * log_ratio);
    //     };

    //     double S0 = 0.;

    //     for (int m = 0; m < M; ++m) {

    //       // 3.4 - log-spaced bin center (geometric mean of shell edges)
    //       const double e_left  = edge(m);
    //       const double e_right = edge(m + 1);
    //       kappa[m]             = std::sqrt(e_left * e_right);   // geometric midpoint

    //       // mode-dependent bin width for amplitude computation
    //       const double dk_m    = e_right - e_left;

    //       // 3.5 - wave amplitude q_m = sqrt(E(kappa_m) * dk_m)
    //       const double ratio = kappa[m] / k0;
    //       const double Em    = a * std::pow(kappa[m], 4.) * std::exp(-2. * ratio * ratio);
    //       q[m]               = std::sqrt(Em * dk_m);
    //       S0                += Em * dk_m;

    //       // 3.6 - unit wave vector k̂_m on unit sphere
    //       const double phi   = 2. * pi * U01(rng);
    //       const double t     = Um11(rng);
    //       const double theta = std::acos(t);
    //       const double s     = std::sin(theta);
    //       khat_x[m] = s * std::cos(phi);
    //       khat_y[m] = s * std::sin(phi);
    //       khat_z[m] = std::cos(theta);

    //       // 3.7 - solenoidal direction vector sigma_m
    //       const double alfa  = 2. * pi * U01(rng);
    //       sigma_x[m] =  std::cos(phi)*std::cos(theta)*std::cos(alfa) - std::sin(phi)*std::sin(alfa);
    //       sigma_y[m] =  std::sin(phi)*std::cos(theta)*std::cos(alfa) + std::cos(phi)*std::sin(alfa);
    //       sigma_z[m] = -std::sin(theta)*std::cos(alfa);
    //       // sigma is already a unit vector by construction — no normalization needed

    //       // 3.8 - random phase psi_m in [-pi/2, pi/2]
    //       phase[m] = Uphase(rng);
    //     }

    //     // Rescaling: correct for quadrature error in discrete spectrum
    //     // This is more important with log spacing because dk_m varies widely.
    //     // The factor enforces int E(k)dk = (3/2)*urms^2.
    //     const double sqrtfac = std::sqrt(Etot / S0);
    //     for (int m = 0; m < M; ++m) q[m] *= sqrtfac;

    //     modes_initialized = true;
    //   }

    //   // Get cell-center coordinates
    //   const double x = cell->getPosition().getX();
    //   const double y = cell->getPosition().getY();
    //   const double z = cell->getPosition().getZ();

    //   double u(0.), v(0.), w(0.);

    //   // 3.9 - Saad velocity series
    //   // u(x) = 2 * SUM_m q_m * cos(kappa_m * k̂_m · x + psi_m) * sigma_m
    //   for (int m = 0; m < M; ++m) {
    //     const double kdotx = kappa[m] * (khat_x[m]*x + khat_y[m]*y + khat_z[m]*z);
    //     const double coeff = 2. * q[m] * std::cos(kdotx + phase[m]);
    //     u += coeff * sigma_x[m];
    //     v += coeff * sigma_y[m];
    //     w += coeff * sigma_z[m];
    //   }

    //   sum_u2v2w2 += u*u + v*v + w*w;
    //   nCells++;

    //   // Step 4 - set values
    //   for (int i = 0; i < numberPhases; ++i) {
    //     cell->getPhase(i)->setU(u);
    //     cell->getPhase(i)->setV(v);
    //     cell->getPhase(i)->setW(w);
    //   }
    //   cell->getMixture()->setU(u);
    //   cell->getMixture()->setV(v);
    //   cell->getMixture()->setW(w);

    // } // end of Saad algorithm
  }
}
//***********************************************************
