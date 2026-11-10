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

#include "EosVDW.h"

//***********************************************************************

EosVDW::EosVDW(std::vector<std::string>& nameParameterEos, int& number) : Eos(number)
{
  nameParameterEos.push_back("Tc");    // Critical temperature of the fluid (Kelvin)
  nameParameterEos.push_back("Pc");    // Critical pressure of the fluid (Pascal)
  nameParameterEos.push_back("M");     // Molar mass (g/mol)
  nameParameterEos.push_back("c_v");   // Heat capacity at constant volume (J / Kg / K) at the critical temperature Tc
  nameParameterEos.push_back("kappa"); // Capillarity coefficient
}

//***********************************************************************

EosVDW::~EosVDW() {}

//***********************************************************************

void EosVDW::assignParametersEos(std::string name, std::vector<double> parametersEos)
{
  assert(parametersEos.size() == 5);

  m_name    = name;
  double Tc = parametersEos[0];                 // Critical temperature of the fluid (Kelvin)
  double Pc = parametersEos[1];                 // Critical pressure of the fluid (Pascal)
  m_r       = 8.314 / parametersEos[2] * 1000.; // Specific gaz constant : R/M (J / K / Kg)
  m_a       = 27. * m_r * m_r * Tc * Tc / (64. * Pc);
  m_b       = m_r * Tc / (8. * Pc);
  m_cv      = parametersEos[3];
  m_gamma   = 1. + m_r / m_cv;
  m_kappa   = parametersEos[4];

  if (m_kappa == 0.) {
    m_IsKappaConst = false;
  }
  else {
    m_IsKappaConst = true;
  }
}

//***********************************************************************

double EosVDW::computeEnergy(const double& density, const double& pressure) const
{
  return (pressure + m_a * density * density) * (1. / density - m_b) / (m_gamma - 1.) - m_a * density;
}

//***********************************************************************

double EosVDW::computeEnergyFromDensityAndTemperature(const double& density, const double& temperature) const
{
  return m_cv * temperature - m_a * density;
}

//***********************************************************************

double EosVDW::computePressure(const double& density, const double& energy) const
{
  return (m_gamma - 1.) * (energy + m_a * density) / (1. / density - m_b) - m_a * density * density;
}

//***********************************************************************

double EosVDW::computeDensity(const double& pressure, const double& temperature) const
{
  // Cubic equation of van der Waals in terms of the compressibility factor Z = P*V/(r*T) => the coefficients are dimensionless
  // The following algorithm assumes that the liquid/vapor state is not in the saturation dome.
  double A = m_a * pressure / (m_r * m_r * temperature * temperature);
  double B = pressure * m_b / (m_r * temperature);
  double C = 1. + B;
  // Solve the depressed cubic
  double p    = A - C * C / 3.;
  double q    = C * (9. * A - 2. * C * C) / 27. - A * B;
  double disc = p * p * p / 27. + q * q / 4.;

  if (disc < 0.) {
    // Viete's formula for multiple roots
    double theta = std::acos(3. * q / (2. * p) * std::sqrt(3. / (-p))) / 3.;
    double Z_g   = C / 3. + 2. * std::sqrt(-p / 3.) * std::cos(theta);
    double Z_l   = C / 3. + 2. * std::sqrt(-p / 3.) * std::cos(theta - 4. * M_PI / 3.);
    double V_g   = Z_g * m_r * temperature / pressure;
    double V_l   = Z_l * m_r * temperature / pressure;
    // Diff of Gibbs potential
    double area = m_r * temperature * std::log((V_g - m_b) / (V_l - m_b)) + m_a * (1. / V_g - 1. / V_l) - pressure * (V_g - V_l);

    if (area < 0.) {
      return 1. / V_l;
    }
    else {
      return 1. / V_g;
    }
  }
  else {
    // Hyperbolic solution for one root (better than Cardano's formula)
    if (p < 0.) {
      double theta = std::acosh(-3. * std::abs(q) / (2. * p) * std::sqrt(3. / -p)) / 3.;
      double Z_r   = C / 3. - 2. * std::abs(q) / q * std::sqrt(-p / 3.) * std::cosh(theta);
      return pressure / (Z_r * m_r * temperature);
    }
    else {
      double theta = std::asinh(3. * q / (2. * p) * std::sqrt(3. / p)) / 3.;
      double Z_r   = C / 3. - 2. * std::sqrt(p / 3.) * std::sinh(theta); // if pressure > 0
      return pressure / (Z_r * m_r * temperature);
    }
  }
}

//***********************************************************************

double EosVDW::computePressureFromDensityAndEntropy(const double& density, const double& entropy) const
{
  return m_r * std::exp(entropy / m_cv) * std::pow(1. / density - m_b, -m_gamma) - m_a * density * density;
}

//***********************************************************************

double EosVDW::computePressureFromDensityAndTemperature(const double& density, const double& temperature) const
{
  return m_r * temperature / (1. / density - m_b) - m_a * density * density;
}

//***********************************************************************

double EosVDW::computeTemperature(const double& density, const double& pressure) const
{
  return (pressure + m_a * density * density) * (1. / density - m_b) / m_r;
}

//***********************************************************************

double EosVDW::computeEntropyFromDensityAndEnergy(const double& density, const double& energy) const
{
  return m_cv * std::log((energy + m_a * density) / m_cv) + m_r * std::log(1. / density - m_b);
}

//***********************************************************************

double EosVDW::computeEntropyFromDensityAndPressure(const double& density, const double& pressure) const
{
  return m_cv * std::log((pressure + m_a * density * density) * (1. / density - m_b) / m_r) + m_r * std::log(1. / density - m_b);
}

//***********************************************************************

double EosVDW::computeEntropyFromDensityAndTemperature(const double& density, const double& temperature) const
{
  return m_cv * std::log(temperature) + m_r * std::log(1. / density - m_b);
}

//***********************************************************************

double EosVDW::computeSquareSoundSpeed(const double& density, const double& energy) const
{
  return m_gamma * (m_gamma - 1.) * (energy + m_a * density) / (1. - density * m_b) / (1. - density * m_b) - 2. * m_a * density;
}

//***********************************************************************

void EosVDW::verifyAndCorrectDensityMax(double& density) const
{
  //AF//TODO// Precompute 1./m_b
  double density_max = 1. / m_b;
  if (density > density_max - 1.e-10) {
    density = density_max - 1.e-10;
  }
}

//***********************************************************************

double EosVDW::kappa(const double& density) const
{
  if (m_IsKappaConst == true) {
    return m_kappa;
  }
  else {
    return 0.25 / density;
  }
}

//****************************************************************************

double EosVDW::kappaPrime(const double& density) const
{
  if (m_IsKappaConst == true) {
    return 0.;
  }
  else {
    double v = 1. / density;
    return -0.25 * v * v;
  }
}

//****************************************************************************

double EosVDW::kappaSecond(const double& density) const
{
  if (m_IsKappaConst == true) {
    return 0.;
  }
  else {
    double v = 1. / density;
    return 0.5 * v * v * v;
  }
}

//***********************************************************************
