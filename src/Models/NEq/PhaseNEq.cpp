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

#include "PhaseNEq.h"

using namespace tinyxml2;

//***************************************************************************

PhaseNEq::PhaseNEq() :
  m_alpha(1.), m_density(0.), m_pressure(0.), m_velocity(0.), m_Y(1.), m_temperature(0.), m_eos(0), m_energy(0.), m_totalEnergy(0.), m_soundSpeed(0.)
{}

//***************************************************************************

PhaseNEq::PhaseNEq(XMLElement* material, Eos* eos, std::string fileName) :
  m_alpha(1.), m_density(0.), m_pressure(0.), m_velocity(0.), m_Y(1.), m_temperature(0.), m_eos(eos), m_energy(0.), m_totalEnergy(0.), m_soundSpeed(0.)
{
  XMLElement* sousElement(material->FirstChildElement("dataFluid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //alpha
  error = sousElement->QueryDoubleAttribute("alpha", &m_alpha);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alpha", fileName, __FILE__, __LINE__);

  //Thermodynamic data reading
  int presenceDensity(0), presencePressure(0), presenceTemperature(0);
  if (sousElement->QueryDoubleAttribute("density", &m_density) == XML_NO_ERROR) presenceDensity = 1;
  if (sousElement->QueryDoubleAttribute("pressure", &m_pressure) == XML_NO_ERROR) presencePressure = 1;
  if (sousElement->QueryDoubleAttribute("temperature", &m_temperature) == XML_NO_ERROR) presenceTemperature = 1;

  //Attribute error gestion
  if (presenceDensity + presencePressure + presenceTemperature != 2)
    throw ErrorXMLAttribut("only two of following is required : density, pressure, temperature", fileName, __FILE__, __LINE__);

  //Thermodynamic reconstruction if needed
  if (presencePressure && presenceTemperature) m_density = m_eos->computeDensity(m_pressure, m_temperature);
  if (presenceDensity && presencePressure) m_temperature = m_eos->computeTemperature(m_density, m_pressure);
  if (presenceDensity && presenceTemperature)
    throw ErrorXMLAttribut("impossible to initialize UEq phase with density and temperature", fileName, __FILE__, __LINE__);

  m_energy     = m_eos->computeEnergy(m_density, m_pressure);
  m_soundSpeed = m_eos->computeSoundSpeed(m_density, m_pressure);

  //velocity
  XMLElement* xmlVelocity(sousElement->FirstChildElement("velocity"));
  if (xmlVelocity == NULL) throw ErrorXMLElement("velocity", fileName, __FILE__, __LINE__);
  double velocityX(0.), velocityY(0.), velocityZ(0.);
  error = xmlVelocity->QueryDoubleAttribute("x", &velocityX);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
  error = xmlVelocity->QueryDoubleAttribute("y", &velocityY);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
  error = xmlVelocity->QueryDoubleAttribute("z", &velocityZ);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName, __FILE__, __LINE__);
  m_velocity.setXYZ(velocityX, velocityY, velocityZ);
}

//***************************************************************************

PhaseNEq::~PhaseNEq() {}

//***************************************************************************

void PhaseNEq::allocateAndCopyPhase(Phase** vecPhase) { *vecPhase = new PhaseNEq(*this); }

//***************************************************************************

void PhaseNEq::copyPhase(Phase& phase)
{
  m_alpha       = phase.getAlpha();
  m_density     = phase.getDensity();
  m_pressure    = phase.getPressure();
  m_velocity    = phase.getVelocity();
  m_Y           = phase.getMassFraction();
  m_temperature = phase.getTemperature();
  m_eos         = phase.getEos();
  m_energy      = phase.getEnergy();
  m_soundSpeed  = phase.getSoundSpeed();
  m_totalEnergy = phase.getTotalEnergy();
}

//***************************************************************************

void PhaseNEq::extendedCalculusPhase(const Coord& /*velocity*/)
{
  m_energy      = m_eos->computeEnergy(m_density, m_pressure);
  m_soundSpeed  = m_eos->computeSoundSpeed(m_density, m_pressure);
  m_totalEnergy = m_energy + 0.5 * m_velocity.squaredNorm();
  m_temperature = m_eos->computeTemperature(m_density, m_pressure);
}

//****************************************************************************

void PhaseNEq::computeMassFraction(const double& density) { m_Y = m_alpha * m_density / std::max(density, epsilonAlphaNull); }

//***************************************************************************

void PhaseNEq::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void PhaseNEq::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseNEq::returnScalar(const int& numVar) const
{
  switch (numVar) {
  case 1:
    if (m_alpha < 1.e-20) {
      return 0.;
    }
    else {
      return m_alpha;
    }
    break;
  case 2:
    if (m_density < 1.e-20) {
      return 1.e-20;
    }
    else {
      return m_density;
    }
    break;
  case 3:
    return m_pressure;
    break;
  case 4:
    if (m_temperature < 1.e-20) {
      return 1.e-20;
    }
    else {
      return m_temperature;
    }
    break;
  case 5:
    if (m_Y < 1.e-20) {
      return 0.;
    }
    else {
      return m_Y;
    }
    break;
  default:
    return 0.;
    break;
  }
}

//***************************************************************************

Coord PhaseNEq::returnVector(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return m_velocity;
    break;
  default:
    return 0;
    break;
  }
}

//***************************************************************************

std::string PhaseNEq::returnNameScalar(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return "Alpha";
    break;
  case 2:
    return "Density";
    break;
  case 3:
    return "Pressure";
    break;
  case 4:
    return "Temperature";
    break;
  case 5:
    return "Mass fraction";
    break;
  default:
    return "NoName";
    break;
  }
}

//***************************************************************************

std::string PhaseNEq::returnNameVector(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return "Velocity";
    break;
  default:
    return "NoName";
    break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhaseNEq::setScalar(const int& numVar, const double& value)
{
  switch (numVar) {
  case 1:
    m_alpha = value;
    break;
  case 2:
    m_density = value;
    break;
  case 3:
    m_pressure = value;
    break;
  case 4:
    m_temperature = value;
    break;
  case 5:
    m_Y = value;
    break;
  default:
    Errors::errorMessage("numVar not found in PhaseNEq::setScalar");
    break;
  }
}

//****************************************************************************

void PhaseNEq::setVector(const int& numVar, const Coord& value)
{
  switch (numVar) {
  case 1:
    m_velocity = value;
    break;
  default:
    Errors::errorMessage("numVar not found in PhaseNEq::setVector");
    break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseNEq::numberOfTransmittedVariables() const
{
  //6 variables (alpha,density,pressure,velocity(3)) + number EOS
  return 7;
}

//***************************************************************************

void PhaseNEq::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_alpha;
  buffer[++counter] = m_density;
  buffer[++counter] = m_pressure;
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();
  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseNEq::fillBuffer(std::vector<double>& dataToSend) const
{
  dataToSend.push_back(m_alpha);
  dataToSend.push_back(m_density);
  dataToSend.push_back(m_pressure);
  dataToSend.push_back(m_velocity.getX());
  dataToSend.push_back(m_velocity.getY());
  dataToSend.push_back(m_velocity.getZ());
  dataToSend.push_back(static_cast<double>(m_eos->getNumber()));
}

//***************************************************************************

void PhaseNEq::getBuffer(double* buffer, int& counter, Eos** eos)
{
  m_alpha    = buffer[++counter];
  m_density  = buffer[++counter];
  m_pressure = buffer[++counter];
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//***************************************************************************

void PhaseNEq::getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos)
{
  m_alpha    = dataToReceive[counter++];
  m_density  = dataToReceive[counter++];
  m_pressure = dataToReceive[counter++];
  m_velocity.setX(dataToReceive[counter++]);
  m_velocity.setY(dataToReceive[counter++]);
  m_velocity.setZ(dataToReceive[counter++]);
  m_eos = eos[static_cast<int>(std::round(dataToReceive[counter++]))];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseNEq::computeSlopesPhase(const Phase& sLeft, const Phase& sRight, const double& distance)
{
  m_alpha    = (sRight.getAlpha() - sLeft.getAlpha()) / distance;
  m_density  = (sRight.getDensity() - sLeft.getDensity()) / distance;
  m_pressure = (sRight.getPressure() - sLeft.getPressure()) / distance;
  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);
}

//***************************************************************************

void PhaseNEq::setToZero()
{
  m_alpha    = 0.;
  m_density  = 0.;
  m_pressure = 0.;
  m_velocity.setX(0.);
  m_velocity.setY(0.);
  m_velocity.setZ(0.);
}

//***************************************************************************

void PhaseNEq::setToMax()
{
  m_alpha    = 1.e15;
  m_density  = 1.e15;
  m_pressure = 1.e15;
  m_velocity.setX(1.e15);
  m_velocity.setY(1.e15);
  m_velocity.setZ(1.e15);
}

//***************************************************************************

void PhaseNEq::extrapolate(const Phase& slope, const double& distance)
{
  m_alpha    += slope.getAlpha() * distance;
  m_density  += slope.getDensity() * distance;
  m_pressure += slope.getPressure() * distance;
  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);
}

//***************************************************************************

void PhaseNEq::limitSlopes(const Phase& slopeGauche, const Phase& slopeDroite, Limiter& globalLimiter, Limiter& volumeFractionLimiter)
{
  m_alpha    = volumeFractionLimiter.limiteSlope(slopeGauche.getAlpha(), slopeDroite.getAlpha());
  m_density  = globalLimiter.limiteSlope(slopeGauche.getDensity(), slopeDroite.getDensity());
  m_pressure = globalLimiter.limiteSlope(slopeGauche.getPressure(), slopeDroite.getPressure());
  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));
}

//****************************************************************************

void PhaseNEq::setMin(const Phase& phase1, const Phase& phase2)
{
  m_alpha    = std::min(phase1.getAlpha(), phase2.getAlpha());
  m_density  = std::min(phase1.getDensity(), phase2.getDensity());
  m_pressure = std::min(phase1.getPressure(), phase2.getPressure());

  m_velocity.setX(std::min(phase1.getVelocity().getX(), phase2.getVelocity().getX()));
  m_velocity.setY(std::min(phase1.getVelocity().getY(), phase2.getVelocity().getY()));
  m_velocity.setZ(std::min(phase1.getVelocity().getZ(), phase2.getVelocity().getZ()));
}

//****************************************************************************

void PhaseNEq::setMax(const Phase& phase1, const Phase& phase2)
{
  m_alpha    = std::max(phase1.getAlpha(), phase2.getAlpha());
  m_density  = std::max(phase1.getDensity(), phase2.getDensity());
  m_pressure = std::max(phase1.getPressure(), phase2.getPressure());

  m_velocity.setX(std::max(phase1.getVelocity().getX(), phase2.getVelocity().getX()));
  m_velocity.setY(std::max(phase1.getVelocity().getY(), phase2.getVelocity().getY()));
  m_velocity.setZ(std::max(phase1.getVelocity().getZ(), phase2.getVelocity().getZ()));
}

//****************************************************************************

void PhaseNEq::computeGradientLimiter(
  const Limiter& globalLimiter, const Phase& phase, const Phase& phaseMin, const Phase& phaseMax, const Phase& slope)
{
  m_alpha    = std::min(m_alpha, globalLimiter.computeGradientLimiter(phase.getAlpha(), phaseMin.getAlpha(), phaseMax.getAlpha(), slope.getAlpha()));
  m_density  = std::min(m_density,
                       globalLimiter.computeGradientLimiter(phase.getDensity(), phaseMin.getDensity(), phaseMax.getDensity(), slope.getDensity()));
  m_pressure = std::min(
    m_pressure, globalLimiter.computeGradientLimiter(phase.getPressure(), phaseMin.getPressure(), phaseMax.getPressure(), slope.getPressure()));

  m_velocity.setX(std::min(m_velocity.getX(),
                           globalLimiter.computeGradientLimiter(
                             phase.getVelocity().getX(), phaseMin.getVelocity().getX(), phaseMax.getVelocity().getX(), slope.getVelocity().getX())));
  m_velocity.setY(std::min(m_velocity.getY(),
                           globalLimiter.computeGradientLimiter(
                             phase.getVelocity().getY(), phaseMin.getVelocity().getY(), phaseMax.getVelocity().getY(), slope.getVelocity().getY())));
  m_velocity.setZ(std::min(m_velocity.getZ(),
                           globalLimiter.computeGradientLimiter(
                             phase.getVelocity().getZ(), phaseMin.getVelocity().getZ(), phaseMax.getVelocity().getZ(), slope.getVelocity().getZ())));
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseNEq::numberOfTransmittedSlopes() const { return 6; }

//***************************************************************************

void PhaseNEq::fillBufferSlopes(double* buffer, int& counter) const
{
  buffer[++counter] = m_alpha;
  buffer[++counter] = m_density;
  buffer[++counter] = m_pressure;
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();
}

//***************************************************************************

void PhaseNEq::getBufferSlopes(double* buffer, int& counter)
{
  m_alpha    = buffer[++counter];
  m_density  = buffer[++counter];
  m_pressure = buffer[++counter];
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseNEq::verifyPhase(const std::string& message) const
{
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    if (m_alpha < 0.) errors.push_back(Errors(message + "too small alpha in verifyPhase"));
    if (m_alpha > 1.) errors.push_back(Errors(message + "too big alpha in verifyPhase"));
    if (m_density < 0.) errors.push_back(Errors(message + "too small density in verifyPhase"));
    m_eos->verifyPressure(m_pressure, message);
  }
  else { // alpha = 0 is desactivated (alpha != 0)
    if (m_alpha <= 1e-15) errors.push_back(Errors(message + "too small alpha in verifyPhase"));
    if (m_alpha >= 1. - 1e-15) errors.push_back(Errors(message + "too big alpha in verifyPhase"));
    if (m_density <= 1.e-15) errors.push_back(Errors(message + "too small density in verifyPhase"));
  }
}

//***************************************************************************

void PhaseNEq::verifyAndCorrectPhase()
{
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    if (m_alpha < 0.) m_alpha = 0.;
    if (m_alpha > 1.) m_alpha = 1.;
    if (m_density <= 1.e-15) m_density = 1.e-15;
  }
  else { // alpha = 0 is desactivated (alpha != 0)
    if (m_alpha < 1e-15) m_alpha = 1e-14;
    if (m_alpha > 1. - 1e-15) m_alpha = 1. - 1e-14;
    if (m_density < 1.e-15) m_density = 1.e-14;
  }
  m_eos->verifyAndModifyPressure(m_pressure);
}

//***************************************************************************

void PhaseNEq::verifyAndCorrectDensityMax(const double& mass) { m_eos->verifyAndCorrectDensityMax(mass, m_alpha, m_density); }

//***************************************************************************

void PhaseNEq::verifyAndCorrectDensityMax() { m_eos->verifyAndCorrectDensityMax(m_density); }

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseNEq::setAlpha(double alpha) { m_alpha = alpha; }

//***************************************************************************

void PhaseNEq::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseNEq::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseNEq::setVelocity(const double& u, const double& v, const double& w)
{
  m_velocity.setX(u);
  m_velocity.setY(v);
  m_velocity.setZ(w);
}

//***************************************************************************

void PhaseNEq::setVelocity(const Coord& vit) { m_velocity = vit; }

//***************************************************************************

void PhaseNEq::setEos(Eos* eos) { m_eos = eos; }

//***************************************************************************

void PhaseNEq::setEnergy(double energy) { m_energy = energy; }

//***************************************************************************

void PhaseNEq::setSoundSpeed(double soundSpeed) { m_soundSpeed = soundSpeed; }

//****************************************************************************
//****************************** OPERATORS ***********************************
//****************************************************************************

void PhaseNEq::changeSign()
{
  m_alpha    = -m_alpha;
  m_density  = -m_density;
  m_pressure = -m_pressure;
  m_velocity = m_velocity * -1.;
}

//***************************************************************************

void PhaseNEq::multiplyAndAdd(const Phase& slopesPhasesTemp, const double& coeff)
{
  m_alpha    += slopesPhasesTemp.getAlpha() * coeff;
  m_density  += slopesPhasesTemp.getDensity() * coeff;
  m_pressure += slopesPhasesTemp.getPressure() * coeff;
  m_velocity += slopesPhasesTemp.getVelocity() * coeff;
}

//***************************************************************************

void PhaseNEq::divide(const double& coeff)
{
  m_alpha    /= coeff;
  m_density  /= coeff;
  m_pressure /= coeff;
  m_velocity /= coeff;
}

//***************************************************************************
