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

//! \file      RelaxationPTMU.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      October 16 2018

#include "RelaxationPTMu.h"

using namespace tinyxml2;

//***********************************************************************

RelaxationPTMu::RelaxationPTMu(){}

//***********************************************************************

RelaxationPTMu::RelaxationPTMu(XMLElement *element, std::string fileName)
{
	XMLElement *sousElement(element->FirstChildElement("dataPTMu"));
	if (sousElement == NULL) throw ErrorXMLElement("dataPTMu", fileName, __FILE__, __LINE__);
	//Collecting attributes
	//---------------------
	m_liq = 0;
	m_vap = 1;
	
}

//***********************************************************************

RelaxationPTMu::~RelaxationPTMu(){}

//***********************************************************************

void RelaxationPTMu::stiffRelaxation(Cell *cell, const int &numberPhases, Prim type) const
{
 	Phase *phase(0);

	if (numberPhases > 2) errors.push_back(Errors("More than 2-phase calculation with evaporation not implemented in RelaxationPTMu::stiffRelaxation", __FILE__, __LINE__));

	//Initial state
	double pStar(0.), Tsat;
	for (int k = 0; k < numberPhases; k++)
	{
		phase = cell->getPhase(k, type);
		TB->ak[k] = phase->getAlpha();
		TB->pk[k] = phase->getPressure();
		TB->rhok[k] = phase->getDensity();
		pStar += TB->ak[k] * TB->pk[k];
    phase->verifyAndCorrectPhase();
		//phase->verifyPhase();
	}
	//cell->extendedCalculus(numberPhases);
	double rho = cell->getMixture()->getDensity();
	double rhoe = rho * cell->getMixture()->getEnergy();

	//Saturation temperature determination
	double dTsat(0.);
	Tsat = cell->getMixture()->computeTsat(cell->getPhase(m_liq)->getEos(), cell->getPhase(m_vap)->getEos(), pStar, &dTsat);

	//evap or not?
	double TL = TB->eos[m_liq]->computeTemperature(TB->rhok[m_liq], pStar);
	double TV = TB->eos[m_vap]->computeTemperature(TB->rhok[m_vap], pStar);
	//if (TB->ak[m_liq] < 1.e-6) return;
	//if (TB->ak[m_vap] < 1.e-6) return;
	//if (TL < Tsat) return;

	//Iterative process for relaxed pressure determination
	double rhoLSat, rhoVSat, drhoLSat, drhoVSat;
	double aLSat, aVSat, daLSat, daVSat;
	double rhoeLSat, rhoeVSat, drhoeLSat, drhoeVSat;
	int iteration(0);
	double f(0.), df(1.);
	do {
		pStar -= f / df; iteration++;
		if (iteration > 50) {
			errors.push_back(Errors("Number of iterations too large in relaxPTMu", __FILE__, __LINE__));
			std::cout << "info cell problematic" << std::endl;
			std::cout << "Liq " << TL << " " << TB->rhok[m_liq] << " " << cell->getMixture()->getPressure() << std::endl;
			std::cout << "Vap " << TV << " " << TB->rhok[m_vap] << " " << cell->getMixture()->getPressure() << std::endl;
			std::cout << Tsat << " " << pStar << std::endl;
			break;
		}
		//Physical pressure?
		for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(pStar); }
		//Liquid-vapor densities calculus
		Tsat = cell->getMixture()->computeTsat(cell->getPhase(m_liq)->getEos(), cell->getPhase(m_vap)->getEos(), pStar, &dTsat);
		rhoLSat = TB->eos[m_liq]->computeDensitySaturation(pStar, Tsat, dTsat, &drhoLSat);
		rhoVSat = TB->eos[m_vap]->computeDensitySaturation(pStar, Tsat, dTsat, &drhoVSat);
		//limit values
		if (rhoLSat <= rho) {
			rhoLSat = rho + 1e-6;
			drhoLSat = 0.;
		}
		if (rhoVSat >= rho) {
			rhoVSat = rho - 1e-6;
			drhoVSat = 0.;
		}
		//Liquid-vapor volume fraction calculus
		aLSat = (rho - rhoVSat) / (rhoLSat - rhoVSat);
		daLSat = (-drhoVSat * (rhoLSat - rhoVSat) - (rho - rhoVSat)*(drhoLSat - drhoVSat)) / ((rhoLSat - rhoVSat)*(rhoLSat - rhoVSat));
		//limit values
		if (aLSat <= 0.) aLSat = 1e-8;
		if (aLSat >= 1.) aLSat = 1. - 1e-8;
		aVSat = 1. - aLSat;
		daVSat = -daLSat;

		f = rhoe; df = 0.;
		//Liquid-vapor couple only
		rhoeLSat = TB->eos[m_liq]->computeDensityEnergySaturation(pStar, rhoLSat, drhoLSat, &drhoeLSat);
		rhoeVSat = TB->eos[m_vap]->computeDensityEnergySaturation(pStar, rhoVSat, drhoVSat, &drhoeVSat);
		f -= (aLSat*rhoeLSat + aVSat * rhoeVSat);
		df -= (daLSat*rhoeLSat + aLSat * drhoeLSat + daVSat * rhoeVSat + aVSat * drhoeVSat);
		f /= rhoe;
		df /= rhoe;
	} while (std::fabs(f) > 1e-10);

	//Cell update
	phase = cell->getPhase(m_liq);
	phase->setAlpha(aLSat);
	phase->setDensity(rhoLSat);
	phase->setPressure(pStar);

	phase = cell->getPhase(m_vap);
	phase->setAlpha(aVSat);
	phase->setDensity(rhoVSat);
	phase->setPressure(pStar);

	cell->fulfillState();
}