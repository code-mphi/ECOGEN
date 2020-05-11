//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

//! \file      GDEntireDomainWithParticularities.cpp
//! \author    K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "GDEntireDomainWithParticularities.h"

//***********************************************************

GDEntireDomainWithParticularities::GDEntireDomainWithParticularities(std::string name, std::vector<Phase*> vecPhases, Mixture *mixture, std::vector<Transport> vecTransports, const int &physicalEntity) :
  GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{}

//***********************************************************

GDEntireDomainWithParticularities::~GDEntireDomainWithParticularities() {}

//***********************************************************

bool GDEntireDomainWithParticularities::belong(Coord &posElement, const int &lvl) const
{
  //1. Laplace pressure initialization
  //----------------------------------
  //return true; //always belong to entire domain

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
  return true; //always belong to entire domain
}

//******************************************************************

void GDEntireDomainWithParticularities::fillIn(Cell *cell, const int &numberPhases, const int &numberTransports) const
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

    //Particularities
    //1. Laplace pressure initialization
    //----------------------------------
    // if (cell->getElement() != 0) {
    //  double pressure(0.);
    //  Coord posElement(cell->getPosition());
    //  double radius;
    //  //radius = posElement.getX(); //1D
    //  radius = std::pow(std::pow(posElement.getX(), 2.) + std::pow(posElement.getY(), 2.), 0.5); //2D
    //  //radius = std::pow(std::pow(posElement.getX(), 2.) + std::pow(posElement.getY(), 2.) + std::pow(posElement.getZ(), 2.), 0.5); //3D
    //  //radius = std::pow(std::pow(posElement.getX() - 153.6e-3, 2.) + std::pow(posElement.getY(), 2.) + std::pow(posElement.getZ(), 2.), 0.5); //3D
    //  pressure = 1.e5 + 1.e-3 / radius * (1.e4 - 1.e5);
    //  //pressure = 1.e5 + 1.e-3 / radius * (4.e3 - 1.e5);
    //  //pressure = 1.e5 + 1.e-3 / radius * (1.e3 - 1.e5);
    //  //pressure = 50.6625e5 + 100.e-6 / radius * (3.55e3 - 50.6625e5);
    //  for (int k = 0; k < numberPhases; k++) { cell->getPhase(k)->setPressure(pressure); }
    //  cell->getMixture()->setPressure(pressure);
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
    Coord perturbedVelocity(cell->getMixture()->getVelocity());
    perturbedVelocity.setX(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getX());
    perturbedVelocity.setY(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getY());
    perturbedVelocity.setZ(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getZ());
    cell->getMixture()->setVelocity(perturbedVelocity);
    
  }
}

//***********************************************************