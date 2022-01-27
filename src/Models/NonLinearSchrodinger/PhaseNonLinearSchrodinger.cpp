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

#include "PhaseNonLinearSchrodinger.h"
#include "../../Eos/Eos.h"

using namespace tinyxml2;

//***************************************************************************

PhaseNonLinearSchrodinger::PhaseNonLinearSchrodinger() : PhaseEulerKorteweg(){}

//***************************************************************************

PhaseNonLinearSchrodinger::PhaseNonLinearSchrodinger(XMLElement* material, Eos* eos, std::string fileName) :
  PhaseEulerKorteweg(material, eos, fileName) {}

//***************************************************************************

PhaseNonLinearSchrodinger::~PhaseNonLinearSchrodinger(){}

//***************************************************************************

void PhaseNonLinearSchrodinger::allocateAndCopyPhase(Phase** vecPhase)
{
  *vecPhase = new PhaseNonLinearSchrodinger(*this);
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseNonLinearSchrodinger::numberOfTransmittedVariables() const
{
  //9 variables (3 scalar + 2*3 vector)
  return 9;
}

//***************************************************************************

void PhaseNonLinearSchrodinger::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_density;
  buffer[++counter] = m_omega;
  buffer[++counter] = m_eta;

  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();

  buffer[++counter] = m_vectorP.getX();
  buffer[++counter] = m_vectorP.getY();
  buffer[++counter] = m_vectorP.getZ();
}

//***************************************************************************

void PhaseNonLinearSchrodinger::fillBuffer(std::vector<double>& dataToSend) const
{
  dataToSend.push_back(m_density);
  dataToSend.push_back(m_omega);
  dataToSend.push_back(m_eta);
  
  dataToSend.push_back(m_velocity.getX());
  dataToSend.push_back(m_velocity.getY());
  dataToSend.push_back(m_velocity.getZ());
    
  dataToSend.push_back(m_vectorP.getX());
  dataToSend.push_back(m_vectorP.getY());
  dataToSend.push_back(m_vectorP.getZ());
}

//***************************************************************************

void PhaseNonLinearSchrodinger::getBuffer(double* buffer, int& counter, Eos** /*eos*/)
{
  m_density = buffer[++counter];
  m_omega = buffer[++counter];
  m_eta = buffer[++counter];
  
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);

  m_vectorP.setX(buffer[++counter]);
  m_vectorP.setY(buffer[++counter]);
  m_vectorP.setZ(buffer[++counter]);
}

//***************************************************************************

void PhaseNonLinearSchrodinger::getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** /*eos*/)
{
  m_density = dataToReceive[counter++];
  m_omega = dataToReceive[counter++];
  m_eta = dataToReceive[counter++];

  m_velocity.setX(dataToReceive[counter++]);
  m_velocity.setY(dataToReceive[counter++]);
  m_velocity.setZ(dataToReceive[counter++]);

  m_vectorP.setX(dataToReceive[counter++]);
  m_vectorP.setY(dataToReceive[counter++]);
  m_vectorP.setZ(dataToReceive[counter++]);
}

//***************************************************************************