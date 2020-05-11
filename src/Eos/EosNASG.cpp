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

//! \file      EosNASG.cpp
//! \author    J. Caze
//! \version   1.1
//! \date      October 11 2019

#include <cmath>
#include <algorithm>
#include "EosNASG.h"

using namespace std;

//***********************************************************************

EosNASG::EosNASG(){}

//***********************************************************************

EosNASG::EosNASG(std::vector<std::string> &nameParameterEos, int &number) :
    Eos(number)
{
  nameParameterEos.push_back("gamma");
  nameParameterEos.push_back("pInf");
  nameParameterEos.push_back("b");
  nameParameterEos.push_back("cv");
  nameParameterEos.push_back("energyRef");
  nameParameterEos.push_back("entropyRef");
}

//***********************************************************************

EosNASG::~EosNASG(){}

//***********************************************************************

//Mod
void EosNASG::sendInfo(double *&data) const
{
	int number = 5;
	data = new double[number];

	data[0] = m_gamma;
	data[1] = m_pInf;
  data[2] = m_b;
	data[3] = m_cv;
	data[4] = m_eRef;
}

//***********************************************************************

void EosNASG::assignParametersEos(std::string name, std::vector<double> parametersEos)
{
    m_name = name;
    assert(parametersEos.size() == 6);
    m_gamma = parametersEos[0];
    m_pInf  = parametersEos[1];
    m_b     = parametersEos[2];
    m_cv    = parametersEos[3];
    m_eRef  = parametersEos[4];
    m_sRef  = parametersEos[5];
}


//Constant methods
//****************
double EosNASG::computeTemperature(const double &density, const double &pressure) const
{
  return (pressure+m_pInf)*(1.-density*m_b)/(m_gamma-1.)/m_cv/std::max(density, epsilonAlphaNull);
}

//***********************************************************************

double EosNASG::computeEnergy(const double &density, const double &pressure) const
{
  return (pressure+m_gamma*m_pInf)*(1.-density*m_b)/(m_gamma-1.)/std::max(density, epsilonAlphaNull) + m_eRef;
}

//***********************************************************************

double EosNASG::computePressure(const double &density, const double &energy) const
{
  return (density*(m_gamma-1.)*(energy-m_eRef)/std::max((1.-density*m_b), epsilonAlphaNull) - m_gamma*m_pInf);
}

//***********************************************************************

double EosNASG::computeDensity(const double &pressure, const double &temperature) const
{
  return (pressure + m_pInf)/std::max((m_gamma-1.)*m_cv*temperature + m_b*(pressure + m_pInf), epsilonAlphaNull);
}

//***********************************************************************

double EosNASG::computeSoundSpeed(const double &density, const double &pressure) const
{
  return sqrt(m_gamma*(pressure+m_pInf)/std::max(density*(1.-density*m_b), epsilonAlphaNull));
}

//***********************************************************************

double EosNASG::computeEntropy(const double &temperature, const double &pressure) const
{
  return m_cv*log(std::pow(temperature, m_gamma) / std::max(std::pow(pressure + m_pInf, m_gamma - 1.), epsilonAlphaNull)) + m_sRef;
}

//***********************************************************************

double EosNASG::computePressureIsentropic(const double &initialPressure, const double &initialDensity, const double &finalDensity) const
{
  return (initialPressure+m_pInf)*std::pow(finalDensity/std::max(initialDensity, epsilonAlphaNull),m_gamma)*std::pow((1.-initialDensity*m_b)/std::max(1.-finalDensity*m_b,epsilonAlphaNull),m_gamma)-m_pInf;
}

//***********************************************************************

double EosNASG::computePressureHugoniot(const double &initialPressure, const double &initialDensity, const double &finalDensity) const
{
  return ((initialDensity*(2.*m_gamma*m_pInf+initialPressure*(m_gamma-1.))-finalDensity*(2.*m_gamma*m_pInf+initialPressure*(m_gamma+1.))+2.*m_b*initialDensity*initialPressure)/std::max(finalDensity*(m_gamma-1.)-initialDensity*(m_gamma+1.)+2.*m_b*initialDensity,epsilonAlphaNull));
}

//***********************************************************************

double EosNASG::computeDensityIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp) const
{
  double finalDensity(initialDensity*std::pow(finalPressure + m_pInf, 1./m_gamma)/std::max(std::pow(initialPressure+m_pInf, 1./m_gamma)*(1.-initialDensity*m_b)+initialDensity*m_b*std::pow(finalPressure+m_pInf, 1./m_gamma), epsilonAlphaNull));
  if (drhodp != NULL) *drhodp = finalDensity*(1-finalDensity*m_b)/std::max((m_gamma*(finalPressure+m_pInf)), epsilonAlphaNull);
  return finalDensity;
}

//***********************************************************************

double EosNASG::computeDensityHugoniot(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp) const
{
  double num((m_gamma+1.)*(finalPressure+m_pInf)+ (m_gamma - 1.)*(initialPressure + m_pInf));
  double denom((m_gamma - 1.)*(finalPressure + m_pInf) + (m_gamma + 1.)*(initialPressure + m_pInf)+2.*m_b*initialDensity*(finalPressure-initialPressure));
  double finalDensity(initialDensity*num/std::max(denom, epsilonAlphaNull));
  if (drhodp != NULL) *drhodp = initialDensity*4.*(m_gamma)*(initialPressure+m_pInf)*(1.-m_b*initialDensity)/std::max((denom*denom), epsilonAlphaNull);
  return finalDensity;
}

//***********************************************************************

double EosNASG::computeDensityPfinal(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp) const
{
  double num((m_gamma)*(finalPressure + m_pInf));
  double denom(num + initialPressure - finalPressure + m_b*initialDensity*(finalPressure-initialPressure));
  double finalDensity(initialDensity*num/std::max(denom, epsilonAlphaNull));
  if (drhodp != NULL) *drhodp = initialDensity*m_gamma*(initialPressure + m_pInf - m_b*initialDensity*(initialPressure+m_pInf)) / std::max((denom*denom), epsilonAlphaNull);
  return finalDensity;
}

//***********************************************************************

double EosNASG::computeEnthalpyIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *dhdp) const
{
  double finalRho, drho;
  finalRho = this->computeDensityIsentropic(initialPressure, initialDensity, finalPressure, &drho);
  double finalEnthalpy(m_gamma*(finalPressure+m_pInf)*(1.-finalRho*m_b) / (m_gamma - 1.) / std::max(finalRho, epsilonAlphaNull) + m_b*finalPressure + m_eRef);
  if (dhdp != NULL) *dhdp = m_gamma / (m_gamma - 1.)*((finalRho - (finalPressure+m_pInf)*drho) / std::max((finalRho*finalRho), epsilonAlphaNull)-m_b)+m_b;
  return finalEnthalpy;
}

//***********************************************************************

double EosNASG::computeDensitySaturation(const double &pressure, const double &Tsat, const double &dTsatdP, double *drhodp) const
{
  double rho;
  if (drhodp != NULL) {
    *drhodp = (m_gamma - 1.)*m_cv*Tsat - (pressure + m_pInf)*(m_gamma - 1.)*m_cv*dTsatdP;
    *drhodp /= std::max((((m_gamma - 1.)*m_cv*Tsat + m_b*(pressure + m_pInf))*((m_gamma - 1.)*m_cv*Tsat + m_b*(pressure + m_pInf))), epsilonAlphaNull);
  }
  rho = (pressure + m_pInf)/std::max(((m_gamma - 1.)*m_cv*Tsat + m_b*(pressure + m_pInf)), epsilonAlphaNull);
  return rho;
}

//***********************************************************************

double EosNASG::computeDensityEnergySaturation(const double &pressure, const double &rho, const double &drhodp, double *drhoedp) const
{
  double rhoe;
  if (drhoedp != NULL) { *drhoedp = (1.-rho*m_b)/(m_gamma-1.) + drhodp*(m_eRef - m_b*(pressure + m_gamma*m_pInf)/(m_gamma-1.)); }
  rhoe = (pressure + m_gamma*m_pInf)*(1.-rho*m_b) / (m_gamma - 1.) + rho*m_eRef;
  return rhoe;
}

//***********************************************************************

void EosNASG::sendSpecialMixtureEos(double &gamPinfOverGamMinusOne, double &eRef, double &oneOverGamMinusOne, double &covolume) const
{
  gamPinfOverGamMinusOne = m_gamma*m_pInf/(m_gamma-1.);
  eRef = m_eRef;
  oneOverGamMinusOne = 1. / (m_gamma - 1.);
  covolume = m_b;
}

//***********************************************************************

double EosNASG::vfpfh(const double &pressure, const double &enthalpy) const
{
  return ((m_gamma - 1.)*(enthalpy - m_b*pressure - m_eRef) / std::max((m_gamma*(pressure+m_pInf)), epsilonAlphaNull) + m_b);
}

//***********************************************************************

double EosNASG::dvdpch(const double &pressure, const double &enthalpy) const
{
  return (-(m_gamma-1.)/m_gamma/((pressure + m_pInf)*(pressure + m_pInf))*(m_b*(pressure+m_pInf) + enthalpy - m_b*pressure - m_eRef));
}

//***********************************************************************

double EosNASG::dvdhcp(const double &pressure, const double &enthalpy) const
{
  return (m_gamma - 1.) / m_gamma / std::max((pressure + m_pInf), epsilonAlphaNull);
}

//***********************************************************************

void EosNASG::verifyPressure(const double &pressure, const std::string &message) const
{
  if (pressure <= -(1. - 1.e-15)*m_pInf + 1.e-15) errors.push_back(Errors(message + " : too low pressure in EosNASG"));
}

//***********************************************************************

void EosNASG::verifyAndModifyPressure(double &pressure) const
{
  if (pressure <= -(1. - 1.e-15)*m_pInf + 1.e-15) pressure = -(1. - 1.e-15)*m_pInf + 1.e-15;
}

//***********************************************************************
