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

#include "RelaxationPTMu.h"

using namespace tinyxml2;

//***********************************************************************

RelaxationPTMu::RelaxationPTMu(XMLElement* element, std::vector<std::string> const& nameEOS, std::string fileName)
{
  XMLElement* subElement(element->FirstChildElement("dataPTMu"));
  if (subElement == NULL) throw ErrorXMLElement("dataPTMu", fileName, __FILE__, __LINE__);
  
  //Collecting attributes
  //---------------------
  std::string liqEosName(subElement->Attribute("liquid"));
  if (liqEosName == "") throw ErrorXMLAttribut("liquid", fileName, __FILE__, __LINE__);
  std::string vapEosName(subElement->Attribute("vapor"));
  if (vapEosName == "") throw ErrorXMLAttribut("vapor", fileName, __FILE__, __LINE__);

  if (nameEOS.size() > 2) throw ErrorXMLMessage("Only two phases can be used with PTMu relaxation", fileName, __FILE__, __LINE__);

  for (unsigned int k = 0; k < nameEOS.size(); k++) {
    if (nameEOS[k] == liqEosName) { m_liq = k; }
    else if (nameEOS[k] == vapEosName) { m_vap = k; }
    else { throw ErrorXMLElement("dataPTMu", fileName, __FILE__, __LINE__); }
  }
}

//***********************************************************************

RelaxationPTMu::~RelaxationPTMu(){}

//***********************************************************************

void RelaxationPTMu::initializeCriticalPressure(Cell *cell)
{
  m_pcrit = cell->getMixture()->computeCriticalPressure(
    cell->getPhase(m_liq)->getEos(), cell->getPhase(m_vap)->getEos()
  );
}

//***********************************************************************

void RelaxationPTMu::relaxation(Cell* cell, const double& /*dt*/, Prim type)
{
  Phase* phase(0);

  if (numberPhases > 2) {
    errors.push_back(Errors("More than 2-phase calculation with evaporation not implemented in RelaxationPTMu::relaxation", __FILE__, __LINE__));
  }

  //Initial state
  double pStar(0.), Tsat;
  for (int k = 0; k < numberPhases; k++)
  {
    phase = cell->getPhase(k, type);
    TB->ak[k] = phase->getAlpha();
    TB->Yk[k] = phase->getMassFraction();
    TB->pk[k] = phase->getPressure();
    TB->rhok[k] = phase->getDensity();
    pStar += TB->ak[k] * TB->pk[k];
    //phase->verifyPhase();
  }
  if (pStar > m_pcrit) {
    warnings.push_back(Errors("Pressure higher than critical pressure in relaxPTMu", __FILE__, __LINE__));
    return;
  }
  //cell->extendedCalculus(numberPhases);
  double rho = cell->getMixture(type)->getDensity();
  double rhoe = rho * cell->getMixture(type)->getEnergy();
  double dTsat(0.);

  // Pure vapor hypothesis
  double rhov, ev, pv, Tv;
  rhov = cell->getMixture(type)->getDensity();
  ev = cell->getMixture(type)->getEnergy();
  pv = cell->getPhase(m_vap, type)->getEos()->computePressure(rhov, ev);
  Tv = cell->getPhase(m_vap, type)->getEos()->computeTemperature(rhov, pv);

  if (pv > 0) {
    Tsat = cell->getMixture(type)->computeTsat(cell->getPhase(m_liq, type)->getEos(), cell->getPhase(m_vap, type)->getEos(), pv, &dTsat);
    if (Tv >= Tsat) {
      // Hypothesis verified
      cell->getPhase(m_vap, type)->setAlpha(1.);
      cell->getPhase(m_liq, type)->setAlpha(0.);
      cell->getPhase(m_vap, type)->setDensity(rhov);

      cell->getPhase(m_vap, type)->setEnergy(ev);
      cell->getPhase(m_vap, type)->setPressure(pv);
      cell->fulfillState(type);
      return;
    }
  }

  // Pure liquid hypothesis
  double rhol, el, pl, Tl;
  rhol = cell->getMixture(type)->getDensity();
  el = cell->getMixture(type)->getEnergy();
  pl = cell->getPhase(m_liq, type)->getEos()->computePressure(rhol, el);
  Tl = cell->getPhase(m_liq, type)->getEos()->computeTemperature(rhol, pl);

  if (pl > 0.) {
    Tsat = cell->getMixture(type)->computeTsat(cell->getPhase(m_liq, type)->getEos(), cell->getPhase(m_vap, type)->getEos(), pl, &dTsat);
    if (Tl <= Tsat) {
      // Hypothesis verified
      cell->getPhase(m_liq, type)->setAlpha(1.);
      cell->getPhase(m_vap, type)->setAlpha(0.);
      cell->getPhase(m_liq, type)->setDensity(rhol);

      cell->getPhase(m_liq, type)->setEnergy(el);
      cell->getPhase(m_liq, type)->setPressure(pl);
      cell->fulfillState(type);
      return;
    }
  }

  //Iterative process for relaxed state determination
  double rhoLSat, rhoVSat, drhoLSat, drhoVSat;
  double aLSat, aVSat, daLSat, daVSat;
  double rhoeLSat, rhoeVSat, drhoeLSat, drhoeVSat;
  int iteration(0);
  double f(0.), df(1.);
  do {
    pStar -= f / df; iteration++;
    
    if (iteration > 50) {
      errors.push_back(Errors("Number of iterations too large in relaxPTMu", __FILE__, __LINE__));
      break;
    }
    
    //Physical pressure?
    for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(pStar); }
    
    //Liquid-vapor densities calculus using phases' EOS
    Tsat = cell->getMixture(type)->computeTsat(cell->getPhase(m_liq, type)->getEos(), cell->getPhase(m_vap, type)->getEos(), pStar, &dTsat);
    rhoLSat = TB->eos[m_liq]->computeDensitySaturation(pStar, Tsat, dTsat, &drhoLSat);
    rhoVSat = TB->eos[m_vap]->computeDensitySaturation(pStar, Tsat, dTsat, &drhoVSat);
    
    //Liquid-vapor volume fraction calculus using mass conservation within a cell
    aLSat = (rho - rhoVSat) / (rhoLSat - rhoVSat);
    daLSat = (-drhoVSat * (rhoLSat - rhoVSat) - (rho - rhoVSat)*(drhoLSat - drhoVSat)) / ((rhoLSat - rhoVSat)*(rhoLSat - rhoVSat));
    aVSat = 1. - aLSat;
    daVSat = -daLSat;

    f = rhoe; df = 0.;
    rhoeLSat = TB->eos[m_liq]->computeDensityEnergySaturation(pStar, rhoLSat, drhoLSat, &drhoeLSat);
    rhoeVSat = TB->eos[m_vap]->computeDensityEnergySaturation(pStar, rhoVSat, drhoVSat, &drhoeVSat);
    
    // Iterative process based on energy conservation within a cell
    f -= (aLSat*rhoeLSat + aVSat * rhoeVSat);
    df -= (daLSat*rhoeLSat + aLSat * drhoeLSat + daVSat * rhoeVSat + aVSat * drhoeVSat);
    f /= rhoe;
    df /= rhoe;
  } while (std::fabs(f) > 1e-10);

  //Cell update
  phase = cell->getPhase(m_liq, type);
  phase->setAlpha(aLSat);
  phase->setDensity(rhoLSat);
  phase->setPressure(pStar);

  phase = cell->getPhase(m_vap, type);
  phase->setAlpha(aVSat);
  phase->setDensity(rhoVSat);
  phase->setPressure(pStar);

  cell->fulfillState(type);
}

//***********************************************************************

