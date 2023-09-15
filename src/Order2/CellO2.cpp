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

#include "CellO2.h"
#include "../Models/Phase.h"

//***********************************************************************

CellO2::CellO2() : Cell(), m_vecPhasesO2(0), m_mixtureO2(0), m_vecTransportsO2(0), m_consSauvegarde(0), m_consTransportsSauvegarde(0) {}

//***********************************************************************

CellO2::CellO2(int lvl) : Cell(lvl), m_vecPhasesO2(0), m_mixtureO2(0), m_vecTransportsO2(0), m_consSauvegarde(0), m_consTransportsSauvegarde(0) {}

//***********************************************************************

CellO2::~CellO2()
{
  for (int k = 0; k < numberPhases; k++) {
    delete m_vecPhasesO2[k];
  }
  delete[] m_vecPhasesO2;
  if (m_vecTransportsO2 != 0) delete[] m_vecTransportsO2;
  delete m_mixtureO2;
  delete m_consSauvegarde;
	if (m_consTransportsSauvegarde != 0) delete[] m_consTransportsSauvegarde;
}

//***********************************************************************

void CellO2::allocate(const std::vector<AddPhys*>& addPhys)
{
  m_vecPhases = new Phase*[numberPhases];
  m_vecPhasesO2 = new Phase*[numberPhases];
  for (int k = 0; k < numberSolids; k++) {
    model->allocatePhaseSolid(&m_vecPhases[k]);
    model->allocatePhaseSolid(&m_vecPhasesO2[k]);
  }
  for (int k = numberSolids; k < numberPhases; k++) {
    model->allocatePhase(&m_vecPhases[k]);
    model->allocatePhase(&m_vecPhasesO2[k]);
  }
  model->allocateMixture(&m_mixture);
  model->allocateMixture(&m_mixtureO2);
  model->allocateCons(&m_cons);
  model->allocateCons(&m_consSauvegarde);
  if (numberTransports > 0) {
    m_vecTransports = new Transport[numberTransports];
    m_consTransports = new Transport[numberTransports];
    m_consTransportsSauvegarde = new Transport[numberTransports];
    m_vecTransportsO2 = new Transport[numberTransports];
  }
  for (unsigned int k = 0; k < addPhys.size(); k++) {
    addPhys[k]->addQuantityAddPhys(this);
  }
}

//***********************************************************************

void CellO2::copyPhase(const int& phaseNumber, Phase* phase)
{
  m_vecPhases[phaseNumber]->copyPhase(*phase);
  m_vecPhasesO2[phaseNumber]->copyPhase(*phase);
}

//***********************************************************************

void CellO2::saveCons()
{
  m_consSauvegarde->setCons(m_cons);
  for (int k = 0; k < numberTransports; k++) { m_consTransportsSauvegarde[k].setValue(m_consTransports[k].getValue()); }
}

//***********************************************************************

void CellO2::getBackCons()
{
  m_cons->setCons(m_consSauvegarde);
  for (int k = 0; k < numberTransports; k++) { m_consTransports[k].setValue(m_consTransportsSauvegarde[k].getValue()); }
}

//***********************************************************************

void CellO2::predictionOrdre2(const double& dt, Symmetry* symmetry)
{
  m_cons->setBufferFlux(*this);                   //fluxBuff receive conservative variables at time n: Un
  symmetry->addSymmetricTerms(this);              //m_cons (sum of fluxes) is incremented by the cylindrical or spherical symmetric terms from primitive variables at time n
  m_cons->multiply(0.5*dt);                       //m_cons is multiplied by dt/2
  m_cons->schemeCorrection(*this);                //Specific correction for non-conservative models (using Un in fluxBuff and fluxes in m_cons)
  m_cons->addFlux(1.);                            //Adding the buffer fluxBuff (Un) to obtain Un+1/2 in m_cons
  m_cons->buildPrim(m_vecPhasesO2, m_mixtureO2);  //Reconstructing m_vecPhasesO2 and m_mixtureO2 from m_cons (Un+1/2)
  
  //Same process for transport (Un construction not needed)
  for (int k = 0; k < numberTransports; k++) {
    m_consTransports[k].multiply(0.5*dt);
    m_vecTransportsO2[k].setValue(m_vecTransports[k].getValue());
    m_vecTransportsO2[k].add(m_consTransports[k].getValue());
  }

  //Relaxations and correction of energies
  //FP//Derniere news, ceci doit etre exclu de cette routine car pas generique selon modele: A mettre dans Run::solveHyperbolicO2 si besoin
  model->relaxations(this, 0.5*dt, vecPhasesO2);
  m_mixtureO2->totalEnergyToInternalEnergy(m_vecQuantitiesAddPhys); //We build the internal energy from the total energy
  m_cons->correctionEnergy(this, vecPhasesO2);
  model->fulfillState(m_vecPhasesO2, m_mixtureO2);
}

//***********************************************************************

void CellO2::fulfillState(Prim type)
{
  //Complete thermodynamical variables
  switch (type) {
  case vecPhases: case restart: //Identical to cell first order
    model->fulfillState(m_vecPhases, m_mixture);
    //This routine is used in different configurations and a short note correspond to each one:
    //- Riemann solver: No need to reconstruct the total energy there because it isn't grabbed during the Riemann problem. The total energy is directly reconstruct there.
    //The reason is to avoid calculations on the gradients of additional physics which are not necessary and furthermore wrongly computed.
    //Note that the capillary energy is not required during the Riemann problem because the models are splitted.
    //- Parallel: No need to reconstruct the total energy there because it is already communicated.
    //Note that the gradients of additional physics would also be wrongly computed if done in the ghost cells.
    //- Relaxation or correction: The total energy doesn't have to be updated there.
    break;
  case vecPhasesO2: //Only usefull for the parallel with second order
    model->fulfillState(m_vecPhasesO2, m_mixtureO2);
    break;
  default: break;
  }
}

//***********************************************************************

Phase* CellO2::getPhase(const int& phaseNumber, Prim type) const
{
  switch (type){
    case vecPhases: return m_vecPhases[phaseNumber]; break;
    case vecPhasesO2: return m_vecPhasesO2[phaseNumber]; break;
    default: return 0; break;
  }
}

//***********************************************************************

Phase** CellO2::getPhases(Prim type) const
{
  switch (type) {
    case vecPhases: return m_vecPhases; break;
    case vecPhasesO2: return m_vecPhasesO2; break;
    default: return 0; break;
  }
}

//***********************************************************************

Mixture* CellO2::getMixture(Prim type) const
{
  switch (type) {
    case vecPhases: return m_mixture; break;
    case vecPhasesO2: return m_mixtureO2; break;
  default: return 0; break;
  }
}

//***********************************************************************

Transport& CellO2::getTransport(const int& numTransport, Prim type) const
{
	switch (type) {
	case vecPhases: return m_vecTransports[numTransport]; break;
	case vecPhasesO2: return m_vecTransportsO2[numTransport]; break;
  default: return m_vecTransports[numTransport]; break; //FP//TODO// trouver un moyen plus intelligent de faire les renvoi par defaut sur objets.
	}
}

//***********************************************************************

Transport* CellO2::getTransports(Prim type) const
{
  switch (type) {
    case vecPhases: return m_vecTransports; break;
    case vecPhasesO2: return m_vecTransportsO2; break;
    default: return 0; break;
  }
}

//***********************************************************************

void CellO2::setTransport(double value, int& numTransport, Prim type)
{
  switch (type) {
  case vecPhases: m_vecTransports[numTransport].setValue(value); break;
  case vecPhasesO2: m_vecTransportsO2[numTransport].setValue(value); break;
  default: break;
  }
}

//***********************************************************************