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

#include "PhaseExtendedNavierStokesKorteweg1D.h"
#include "../../Eos/Eos.h"

using namespace tinyxml2;

//***************************************************************************
double lambdaEK, betaEK, tauEK, epsEK;

PhaseExtendedNavierStokesKorteweg1D::PhaseExtendedNavierStokesKorteweg1D() :
  m_density(0.),
  m_omega(0.),
  m_densitySurEta(0.),
  m_pressure(0.),
  m_eos(0),
  m_energy(0.),
  m_entropy(0.),
  m_temperature(0.),
  m_totalEnergy(0.),
  m_squareSoundSpeed(0.)
{
  m_velocity.setXYZ(0., 0., 0.);
  m_vectorP.setXYZ(0., 0., 0.);
}

//***************************************************************************

PhaseExtendedNavierStokesKorteweg1D::PhaseExtendedNavierStokesKorteweg1D(XMLElement* material, Eos* eos, std::string fileName) :
  m_density(0.),
  m_omega(0.),
  m_densitySurEta(0.),
  m_pressure(0.),
  m_velocity(0),
  m_vectorP(0),
  m_eos(eos),
  m_energy(0.),
  m_entropy(0.),
  m_temperature(0.),
  m_totalEnergy(0.),
  m_squareSoundSpeed(0.)
{
  XMLElement* sousElement(material->FirstChildElement("dataFluid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);

  //Attributes reading
  //------------------
  XMLError error;
  //Density
  //error = sousElement->QueryDoubleAttribute("density", &m_density);
  //if (error != XML_NO_ERROR) throw ErrorXMLAttribut("density", fileName, __FILE__, __LINE__);

  //Thermodynamic data reading
  int presenceDensity(0), presencePressure(0), presenceTemperature(0);
  if (sousElement->QueryDoubleAttribute("density", &m_density) == XML_NO_ERROR) presenceDensity = 1;
  if (sousElement->QueryDoubleAttribute("pressure", &m_pressure) == XML_NO_ERROR) presencePressure = 1;
  if (sousElement->QueryDoubleAttribute("temperature", &m_temperature) == XML_NO_ERROR) presenceTemperature = 1;

  //Attribute error gestion
  if (presenceDensity + presencePressure + presenceTemperature != 2)
    throw ErrorXMLAttribut("only two of following is required : density, pressure, temperature", fileName, __FILE__, __LINE__);

  //Thermodynamic reconstruction
  if (presencePressure && presenceTemperature) {
    m_density = m_eos->computeDensity(m_pressure, m_temperature);
    m_entropy = m_eos->computeEntropyFromDensityAndPressure(m_density, m_pressure);
  }
  if (presenceDensity && presencePressure) {
    m_entropy = m_eos->computeEntropyFromDensityAndPressure(m_density, m_pressure);
  }
  if (presenceDensity && presenceTemperature) {
    m_pressure = m_eos->computePressureFromDensityAndTemperature(m_density, m_temperature);
    m_entropy  = m_eos->computeEntropyFromDensityAndTemperature(m_density, m_temperature);
  }

  //Velocity
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

PhaseExtendedNavierStokesKorteweg1D::~PhaseExtendedNavierStokesKorteweg1D() {}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::allocateAndCopyPhase(Phase** vecPhase) { *vecPhase = new PhaseExtendedNavierStokesKorteweg1D(*this); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::copyPhase(Phase& phase)
{
  m_density          = phase.getDensity();
  m_omega            = phase.getOmega();
  m_densitySurEta    = phase.getDensitySurEta();
  m_pressure         = phase.getPressure();
  m_velocity         = phase.getVelocity();
  m_vectorP          = phase.getVectorP();
  m_eos              = phase.getEos();
  m_energy           = phase.getEnergy();
  m_entropy          = phase.getEntropy();
  m_totalEnergy      = phase.getTotalEnergy();
  m_squareSoundSpeed = phase.getSquareSoundSpeed();
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::extendedCalculusPhase(const Coord& /*velocity*/)
{
  m_pressure         = m_eos->computePressureFromDensityAndEntropy(m_density, m_entropy);
  m_energy           = m_eos->computeEnergy(m_density, m_pressure);
  m_squareSoundSpeed = m_eos->computeSquareSoundSpeed(m_density, m_energy);

  m_totalEnergy = m_energy + 0.5 * m_velocity.squaredNorm() + 0.5 * betaEK * m_omega * m_omega +
                  0.5 * lambdaEK * (1. - 1. / m_densitySurEta) * (1. - 1. / m_densitySurEta) +
                  0.5 * epsEK * epsEK * m_eos->kappa(m_density) / m_density * m_vectorP.squaredNorm();
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
  m_vectorP.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
  m_vectorP.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseExtendedNavierStokesKorteweg1D::returnScalar(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return m_density;
    //return 1. / m_density - 1.; //QD//For polynomial EOS
    break;
  case 2:
    return m_omega;
    break;
  case 3:
    return m_densitySurEta;
    break;
  case 4:
    return m_pressure;
    break;
  case 5:
    return m_entropy;
    break;
  case 6:
    return m_pressure + 0.5 * epsEK * epsEK * (m_density * m_eos->kappaPrime(m_density) + m_eos->kappa(m_density)) * m_vectorP.squaredNorm() +
           lambdaEK * m_density / m_densitySurEta * (1. - 1. / m_densitySurEta);
    break;
  case 7:
    return m_totalEnergy;
    break;
  case 8:
    return m_squareSoundSpeed;
    break;
  case 9:
    return m_velocity.getX();
    break;
  case 10:
    return m_vectorP.getX();
    break;
  case 11:
    return m_eos->computeTemperature(m_density, m_pressure);
    break;
  default:
    return 0.;
    break;
  }
}

//***************************************************************************

Coord PhaseExtendedNavierStokesKorteweg1D::returnVector(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return m_velocity;
    break;
  case 2:
    return m_vectorP;
    break;
  default:
    return 0;
    break;
  }
}

//***************************************************************************

std::string PhaseExtendedNavierStokesKorteweg1D::returnNameScalar(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return "Density";
    //return "Epsilon"; //QD//For polynomial EOS
    break;
  case 2:
    return "Psi";
    break;
  case 3:
    return "Density / Eta";
    break;
  case 4:
    return "Pressure";
    break;
  case 5:
    return "Entropy";
    break;
  case 6:
    return "Total Pressure";
    break;
  case 7:
    return "Total Energy";
    break;
  case 8:
    return "Square sound speed";
    break;
  case 9:
    return "Velocity_x";
    break;
  case 10:
    return "Phi_x";
    break;
  case 11:
    return "Temperature";
    break;
  default:
    return "NoName";
    break;
  }
}

//***************************************************************************

std::string PhaseExtendedNavierStokesKorteweg1D::returnNameVector(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return "Velocity";
    break;
  case 2:
    return "VectorP";
    break;
  default:
    return "NoName";
    break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setScalar(const int& numVar, const double& value)
{
  switch (numVar) {
  case 1:
    m_density = value;
    break;
  case 2:
    m_omega = value;
    break;
  case 3:
    m_densitySurEta = value;
    break;
  case 4:
    m_pressure = value;
    break;
  case 5:
    m_entropy = value;
    break;
  default:
    Errors::errorMessage("numVar not found in PhaseExtendedNavierStokesKorteweg1D::setScalar");
    break;
  }
}

//****************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setVector(const int& numVar, const Coord& value)
{
  switch (numVar) {
  case 1:
    m_velocity = value;
    break;
  case 2:
    m_vectorP = value;
    break;
  default:
    Errors::errorMessage("numVar not found in PhaseExtendedNavierStokesKorteweg1D::setVector");
    break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseExtendedNavierStokesKorteweg1D::numberOfTransmittedVariables() const
{
  //11 variables (4 scalar + 2*3 vector + number EOS)
  return 11;
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_density;
  buffer[++counter] = m_omega;
  buffer[++counter] = m_densitySurEta;
  buffer[++counter] = m_entropy;

  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();

  buffer[++counter] = m_vectorP.getX();
  buffer[++counter] = m_vectorP.getY();
  buffer[++counter] = m_vectorP.getZ();

  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::fillBuffer(std::vector<double>& dataToSend) const
{
  dataToSend.push_back(m_density);
  dataToSend.push_back(m_omega);
  dataToSend.push_back(m_densitySurEta);
  dataToSend.push_back(m_entropy);

  dataToSend.push_back(m_velocity.getX());
  dataToSend.push_back(m_velocity.getY());
  dataToSend.push_back(m_velocity.getZ());

  dataToSend.push_back(m_vectorP.getX());
  dataToSend.push_back(m_vectorP.getY());
  dataToSend.push_back(m_vectorP.getZ());

  dataToSend.push_back(static_cast<double>(m_eos->getNumber()));
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::getBuffer(double* buffer, int& counter, Eos** eos)
{
  m_density       = buffer[++counter];
  m_omega         = buffer[++counter];
  m_densitySurEta = buffer[++counter];
  m_entropy       = buffer[++counter];

  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);

  m_vectorP.setX(buffer[++counter]);
  m_vectorP.setY(buffer[++counter]);
  m_vectorP.setZ(buffer[++counter]);

  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos)
{
  m_density       = dataToReceive[counter++];
  m_omega         = dataToReceive[counter++];
  m_densitySurEta = dataToReceive[counter++];
  m_entropy       = dataToReceive[counter++];

  m_velocity.setX(dataToReceive[counter++]);
  m_velocity.setY(dataToReceive[counter++]);
  m_velocity.setZ(dataToReceive[counter++]);

  m_vectorP.setX(dataToReceive[counter++]);
  m_vectorP.setY(dataToReceive[counter++]);
  m_vectorP.setZ(dataToReceive[counter++]);

  m_eos = eos[static_cast<int>(dataToReceive[counter++])];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::computeSlopesPhase(const Phase& sLeft, const Phase& sRight, const double& distance)
{
  m_density       = (sRight.getDensity() - sLeft.getDensity()) / distance;
  m_omega         = (sRight.getOmega() - sLeft.getOmega()) / distance;
  m_densitySurEta = (sRight.getDensitySurEta() - sLeft.getDensitySurEta()) / distance;
  m_entropy       = (sRight.getEntropy() - sLeft.getEntropy()) / distance;

  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);

  m_vectorP.setX((sRight.getVectorP().getX() - sLeft.getVectorP().getX()) / distance);
  m_vectorP.setY((sRight.getVectorP().getY() - sLeft.getVectorP().getY()) / distance);
  m_vectorP.setZ((sRight.getVectorP().getZ() - sLeft.getVectorP().getZ()) / distance);
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setToZero()
{
  m_density       = 0.;
  m_omega         = 0.;
  m_densitySurEta = 0.;
  m_entropy       = 0.;
  m_velocity.setX(0.);
  m_velocity.setY(0.);
  m_velocity.setZ(0.);
  m_vectorP.setX(0.);
  m_vectorP.setY(0.);
  m_vectorP.setZ(0.);
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::extrapolate(const Phase& slope, const double& distance)
{
  m_density       += slope.getDensity() * distance;
  m_omega         += slope.getOmega() * distance;
  m_densitySurEta += slope.getDensitySurEta() * distance;
  m_entropy       += slope.getEntropy() * distance;

  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);

  m_vectorP.setX(m_vectorP.getX() + slope.getVectorP().getX() * distance);
  m_vectorP.setY(m_vectorP.getY() + slope.getVectorP().getY() * distance);
  m_vectorP.setZ(m_vectorP.getZ() + slope.getVectorP().getZ() * distance);
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::limitSlopes(const Phase& slopeGauche,
                                                      const Phase& slopeDroite,
                                                      Limiter& globalLimiter,
                                                      Limiter& /*volumeFractionLimiter*/)
{
  m_density       = globalLimiter.limiteSlope(slopeGauche.getDensity(), slopeDroite.getDensity());
  m_omega         = globalLimiter.limiteSlope(slopeGauche.getOmega(), slopeDroite.getOmega());
  m_densitySurEta = globalLimiter.limiteSlope(slopeGauche.getDensitySurEta(), slopeDroite.getDensitySurEta());
  m_entropy       = globalLimiter.limiteSlope(slopeGauche.getEntropy(), slopeDroite.getEntropy());

  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));

  m_vectorP.setX(globalLimiter.limiteSlope(slopeGauche.getVectorP().getX(), slopeDroite.getVectorP().getX()));
  m_vectorP.setY(globalLimiter.limiteSlope(slopeGauche.getVectorP().getY(), slopeDroite.getVectorP().getY()));
  m_vectorP.setZ(globalLimiter.limiteSlope(slopeGauche.getVectorP().getZ(), slopeDroite.getVectorP().getZ()));
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseExtendedNavierStokesKorteweg1D::numberOfTransmittedSlopes() const { return 10; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::fillBufferSlopes(double* buffer, int& counter) const
{
  buffer[++counter] = m_density;
  buffer[++counter] = m_omega;
  buffer[++counter] = m_densitySurEta;
  buffer[++counter] = m_entropy;

  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();

  buffer[++counter] = m_vectorP.getX();
  buffer[++counter] = m_vectorP.getY();
  buffer[++counter] = m_vectorP.getZ();
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::getBufferSlopes(double* buffer, int& counter)
{
  m_density       = buffer[++counter];
  m_omega         = buffer[++counter];
  m_densitySurEta = buffer[++counter];
  m_entropy       = buffer[++counter];

  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);

  m_vectorP.setX(buffer[++counter]);
  m_vectorP.setY(buffer[++counter]);
  m_vectorP.setZ(buffer[++counter]);
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::verifyPhase(const std::string& message) const
{
  if (m_density <= 1.e-10) errors.push_back(Errors(message + "too small density in verifyPhase"));
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::verifyAndCorrectPhase()
{
  if (m_density < 1.e-10) m_density = 1.e-10;
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::verifyAndCorrectDensityMax() { m_eos->verifyAndCorrectDensityMax(m_density); }

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setOmega(const double& omega) { m_omega = omega; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setDensitySurEta(const double& densitySurEta) { m_densitySurEta = densitySurEta; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setVelocity(const double& u, const double& v, const double& w) { m_velocity.setXYZ(u, v, w); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setVelocity(const Coord& vit) { m_velocity = vit; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setU(const double& u) { m_velocity.setX(u); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setV(const double& v) { m_velocity.setY(v); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setW(const double& w) { m_velocity.setZ(w); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setVectorP(const double& Px, const double& Py, const double& Pz) { m_vectorP.setXYZ(Px, Py, Pz); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setVectorP(const Coord& vecP) { m_vectorP = vecP; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setVectorPX(const double& Px) { m_vectorP.setX(Px); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setVectorPY(const double& Py) { m_vectorP.setY(Py); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setVectorPZ(const double& Pz) { m_vectorP.setZ(Pz); }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setEos(Eos* eos) { m_eos = eos; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setEntropy(double entropy) { m_entropy = entropy; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setEnergy(double energy) { m_energy = energy; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setTotalEnergy(double totalEnergy) { m_totalEnergy = totalEnergy; }

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::setSquareSoundSpeed(double squareSoundSpeed) { m_squareSoundSpeed = squareSoundSpeed; }

//****************************************************************************
//***************************** OPERATORS ************************************
//****************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::changeSign()
{
  m_density       = -m_density;
  m_omega         = -m_omega;
  m_densitySurEta = -m_densitySurEta;
  m_entropy       = -m_entropy;
  m_velocity      = m_velocity * -1.;
  m_vectorP       = m_vectorP * -1;
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::multiplyAndAdd(const Phase& slopesPhasesTemp, const double& coeff)
{
  m_density       += slopesPhasesTemp.getDensity() * coeff;
  m_omega         += slopesPhasesTemp.getOmega() * coeff;
  m_densitySurEta += slopesPhasesTemp.getDensitySurEta() * coeff;
  m_entropy       += slopesPhasesTemp.getEntropy() * coeff;
  m_velocity      += slopesPhasesTemp.getVelocity() * coeff;
  m_vectorP       += slopesPhasesTemp.getVectorP() * coeff;
}

//***************************************************************************

void PhaseExtendedNavierStokesKorteweg1D::divide(const double& coeff)
{
  m_density       /= coeff;
  m_omega         /= coeff;
  m_densitySurEta /= coeff;
  m_entropy       /= coeff;
  m_velocity      /= coeff;
  m_vectorP       /= coeff;
}

//***************************************************************************
