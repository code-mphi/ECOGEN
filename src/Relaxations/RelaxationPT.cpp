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

//! \file      RelaxationPT.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      October 16 2018

#include "RelaxationPT.h"

//***********************************************************************

RelaxationPT::RelaxationPT(){}

//***********************************************************************

RelaxationPT::~RelaxationPT(){}

//***********************************************************************

void RelaxationPT::stiffRelaxation(Cell *cell, const int &numberPhases, Prim type) const
{
  //Is the pressure-Temperature relaxation procedure necessary?
  //If alpha = 0 is activated, a test is done to know if the relaxation procedure is necessary or not
  //Else, i.e. alpha = 0 is desactivated (alpha != 0), the relaxation procedure is always done (relax = true)
  bool relax(true);
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    for (int k = 0; k < numberPhases; k++) {
      if (cell->getPhase(k, type)->getAlpha() >(1. - 1.e-5)) relax = false;
    }
  }

  if (relax) {
    //Restrictions //FP//TODO// to improve
    if (numberPhases > 2) Errors::errorMessage("More than two phases not permitted in RelaxationPT::stiffRelaxation");
    for (int k = 0; k < numberPhases; k++) {
      if (cell->getPhase(k, type)->getEos()->getType() != "IG" && cell->getPhase(k, type)->getEos()->getType() != "SG") { Errors::errorMessage("Only IG and SG permitted in RelaxationPT::stiffRelaxation"); }
    }
       
    //Relaxaed pressure for 2 phases (SG or IG EOS)
    double pStar = analyticalPressure(cell, numberPhases, type);
    //Temperature for SG or IG
    double TStar = analyticalTemperature(pStar, cell, numberPhases, type);

    for (int k = 0; k < numberPhases; k++) {
      TB->rhokS[k] = TB->eos[k]->computeDensity(pStar,TStar);
      TB->akS[k] = cell->getPhase(k, type)->getAlpha()*cell->getPhase(k, type)->getDensity()/ TB->rhokS[k];
      //phase->verifyPhase();
    }

    //Cell update
    for (int k = 0; k < numberPhases; k++) {
      cell->getPhase(k, type)->setAlpha(TB->akS[k]);
      cell->getPhase(k, type)->setDensity(TB->rhokS[k]);
      cell->getPhase(k, type)->setPressure(pStar);
    }
    cell->getMixture(type)->setPressure(pStar);
  }
}

//***********************************************************************

double RelaxationPT::analyticalPressure(Cell *cell, const int &numberPhases, Prim type) const
{
  //Restrictions
  if (numberPhases > 2) Errors::errorMessage("More than two phases not permitted in RelaxationPT::analyticalPressure");
  for (int k = 0; k < numberPhases; k++) {
    if (cell->getPhase(k)->getEos()->getType() != "IG" && cell->getPhase(k)->getEos()->getType() != "SG") { Errors::errorMessage("Only IG and SG permitted in RelaxationPT::analyticalPressure" + cell->getPhase(k)->getEos()->getType()); }
  }

  double e0(cell->getMixture(type)->getEnergy());
  double rho0(cell->getMixture(type)->getDensity());

  //Formulae of pressure for 2 phases goverened by SG EOS (Le Martelot, 2013, thesis)
  double gamma1 = cell->getPhase(0)->getEos()->getGamma();
  double pInf1 = cell->getPhase(0)->getEos()->getPInf();
  double cv1 = cell->getPhase(0)->getEos()->getCv();
  double e01 = cell->getPhase(0)->getEos()->getERef();
  double Y1 = cell->getPhase(0,type)->getY();
  double gamma2 = cell->getPhase(1)->getEos()->getGamma();
  double pInf2 = cell->getPhase(1)->getEos()->getPInf();
  double cv2 = cell->getPhase(1)->getEos()->getCv();
  double e02 = cell->getPhase(1)->getEos()->getERef();
  double Y2 = cell->getPhase(1, type)->getY();

  double q = Y1 * e01 + Y2 * e02;
  double cvMel = Y1 * cv1 + Y2 * cv2;
  double A1 = Y1 * (gamma1 - 1.)*cv1 / cvMel * (rho0*(e0 - q) - pInf1);
  double A2 = Y2 * (gamma2 - 1.)*cv2 / cvMel * (rho0*(e0 - q) - pInf2);
  double pressure = 0.5*(A1 + A2 - (pInf1 + pInf2)) + sqrt(0.25*(A2 - A1 - (pInf2 - pInf1))*(A2 - A1 - (pInf2 - pInf1)) + A1 * A2);
  return pressure;

}

//***********************************************************************

double RelaxationPT::analyticalTemperature(double pressure, Cell *cell, const int &numberPhases, Prim type) const
{
  //Restrictions
  for (int k = 0; k < numberPhases; k++) {
    if (cell->getPhase(k)->getEos()->getType() != "IG" && cell->getPhase(k)->getEos()->getType() != "SG") { Errors::errorMessage("Only IG and SG permitted in RelaxationPT::analyticalPressure" + cell->getPhase(k)->getEos()->getType()); }
  }

  double rho0(cell->getMixture(type)->getDensity());

  //Formulae for phases goverened by SG EOS (Le Martelot, 2013, phd thesis)
  double gammak, pInfk, cvk, Yk;
  double temperature(0.);
  for (int k = 0; k < numberPhases; k++) {
    gammak = cell->getPhase(k)->getEos()->getGamma();
    pInfk = cell->getPhase(k)->getEos()->getPInf();
    cvk = cell->getPhase(k)->getEos()->getCv();
    Yk = cell->getPhase(k, type)->getY();
    temperature += Yk * (gammak - 1.)*cvk / (pressure + pInfk);
  }
  temperature = 1. / (temperature*rho0);

  return temperature;
}

//***********************************************************************