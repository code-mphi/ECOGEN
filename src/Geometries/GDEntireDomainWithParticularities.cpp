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

//***********************************************************

GDEntireDomainWithParticularities::GDEntireDomainWithParticularities(std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, const int& physicalEntity) :
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
}

//******************************************************************

void GDEntireDomainWithParticularities::fillIn(Cell* cell) const
{
  //As basic fillIn: Test if the cell belongs to the geometrical domain
  bool belongs(true);
  if (cell->getElement() != 0) {
    Coord coord(cell->getPosition());
    if (!this->belong(coord, cell->getLvl())) { belongs = false; }
    //Test if the cell belongs to physical mesh entity (for unstructured meshes)
    if (cell->getElement()->getAppartenancePhysique() > 0 && m_physicalEntity > 0) {
      if (cell->getElement()->getAppartenancePhysique() != m_physicalEntity) { belongs = false; }
    }
  }

  if (belongs) {
    for (int k = 0; k < numberPhases; k++) { cell->copyPhase(k, m_vecPhases[k]); }
    cell->copyMixture(m_mixture);
    for (int k = 0; k < numberTransports; k++) { cell->setTransport(m_vecTransports[k].getValue(), k); }
    if(m_physicalEntity == -1){ cell->setWall(true); }
    else{ cell->setWall(false); }

    //Particularities
    //1. Laplace pressure initialization
    //----------------------------------
    // if (cell->getElement() != 0) {
    //   double pressure(0.);
    //   Coord posElement(cell->getPosition());
    //   double radius;
    //   //radius = posElement.getX(); //1D
    //   //radius = std::pow(std::pow(posElement.getX(), 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   //radius = std::pow(std::pow(posElement.getX() - 0.75e-3, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   radius = std::pow(std::pow(posElement.getX() - 2.e-4, 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //   //radius = std::pow(std::pow(posElement.getX(), 2.) + std::pow(posElement.getY(), 2.) + std::pow(posElement.getZ(), 2.), 0.5); //3D
    //   //radius = std::pow(std::pow(posElement.getX() - 153.6e-3, 2.) + std::pow(posElement.getY(), 2.) + std::pow(posElement.getZ(), 2.), 0.5); //3D
    //   //pressure = 1.e5 + 1.e-3 / radius * (1.e4 - 1.e5);
    //   //pressure = 1.e5 + 1.e-3 / radius * (4.e3 - 1.e5);
    //   //pressure = 1.e5 + 1.e-3 / radius * (1.e3 - 1.e5);
    //   //pressure = 353.e5 + 1.e-4 / radius * (1.e5 - 353.e5);
    //   pressure = 50.6625e5 + 1.e-4 / radius * (3.55e3 - 50.6625e5);
    //   for (int k = 0; k < numberPhases; k++) { cell->getPhase(k)->setPressure(pressure); }
    //   cell->getMixture()->setPressure(pressure);
    // }

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
    // Coord perturbedVelocity(cell->getMixture()->getVelocity());
    // perturbedVelocity.setX(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getX());
    // perturbedVelocity.setY(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getY());
    // perturbedVelocity.setZ(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getZ());
    // cell->getMixture()->setVelocity(perturbedVelocity);
    // }

    //5. Rayleigh-Taylor instability
    //------------------------------
    // // Sinus shape parameters
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
    //   double posX(cell->getPosition().getX());
    //   double p0(1.01325e5), pS(35e6), soundSpeed(1625.), density(1000.);
    //   double shockFront(7.5e-3);

    //   if (posX < shockFront) {
    //     double r = posX - shockFront;
    //     pressure = p0 + 2 * pS * exp(beta * r / soundSpeed) * std::cos(- omega * r / soundSpeed + M_PI / 3.);
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
  }
}

//***********************************************************
