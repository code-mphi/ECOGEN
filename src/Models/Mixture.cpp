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

#include "Mixture.h"

int numberScalarsMixture;

//***************************************************************************

Mixture::Mixture(){}

//***************************************************************************

Mixture::~Mixture(){}

//***************************************************************************

void Mixture::printMixture(std::ofstream &fileStream) const
{
  //Scalar variables
  for (int var = 1; var <= this->getNumberScalars(); var++) {
    fileStream << this->returnScalar(var) << " ";
  }
  //Vector variables
  for (int var = 1; var <= this->getNumberVectors(); var++) {
    fileStream << this->returnVector(var).norm() << " ";
  } 
}

//***************************************************************************

double Mixture::computeTsat(const Eos* eosLiq, const Eos* eosVap, const double& pressure, double* dTsat)
{
  //Restrictions //FP//TODO// to improve
  if (eosLiq->getType() != TypeEOS::IG && eosLiq->getType() != TypeEOS::SG) { Errors::errorMessage("Only IG and SG permitted in thermal equilibrium model: MixPTUEq::computeTsat"); }
  if (eosVap->getType() != TypeEOS::IG && eosVap->getType() != TypeEOS::SG) { Errors::errorMessage("Only IG and SG permitted in thermal equilibrium model: MixPTUEq::computeTsat"); }

  double gammaL = eosLiq->getGamma();
  double pInfL = eosLiq->getPInf();
  double cvL = eosLiq->getCv();
  double e0L = eosLiq->getERef();
  double s0L = eosLiq->getSRef();

  double gammaV = eosVap->getGamma();
  double pInfV = eosVap->getPInf();
  double cvV = eosVap->getCv();
  double e0V = eosVap->getERef();
  double s0V = eosVap->getSRef();

  double A, B, C, D;
  A = (gammaL*cvL - gammaV*cvV + s0V - s0L) / (gammaV*cvV - cvV);
  B = (e0L - e0V) / (gammaV*cvV - cvV);
  C = (gammaV*cvV - gammaL*cvL) / (gammaV*cvV - cvV);
  D = (gammaL*cvL - cvL) / (gammaV*cvV - cvV);

  //iterative process to catch saturation temperature
  int iteration(0);
  double Tsat(0.1*B / C);
  double f(0.), df(1.);
  do {
    Tsat -= f / df; iteration++;
    if (iteration > 50) {
      errors.push_back(Errors("number iterations trop grand dans recherche Tsat", __FILE__, __LINE__));
      break;
    }
    f = A + B / Tsat + C*log(Tsat) - log(pressure + pInfV) + D*log(pressure + pInfL);
    df = C / Tsat - B / (Tsat*Tsat);
  } while (std::fabs(f)>1e-10);

  double dfdp = -1. / (pressure + pInfV) + D / (pressure + pInfL);
  if (dTsat != 0) *dTsat = -dfdp / df;
  return Tsat;
}

//***************************************************************************