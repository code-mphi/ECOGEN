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

#ifndef INPUT_H
#define INPUT_H

#include <sstream>
#include <cassert>

#include "../libTierces/tinyxml2.h"
#include "../Eos/HeaderEquationOfState.h"
#include "../Geometries/HeaderGeometricalDomain.h"
#include "../Sources/HeaderSources.h"
#include "../Symmetries/HeaderSymmetry.h"
#include "../Order2/HeaderLimiter.h"
#include "../BoundConds/HeaderBoundCond.h"
#include "../AdditionalPhysics/HeaderQuantitiesAddPhys.h"
#include "../Errors.h"
#include "../Tools.h"

class Input;

#include "../Run.h"
#include "Output.h"

class Input
{
  public:
    Input(Run *run);
    virtual ~Input();

    void lectureInputXML(std::vector<GeometricalDomain*>& domains, std::vector<BoundCond*>& boundCond);

    void entreeMain(std::string casTest);
    void entreeMesh(std::string casTest);
    void entreeModel(std::string casTest);
    Eos* entreeEOS(std::string EOS, int& numberEOS);
    void entreeConditionsInitiales(std::string casTest, std::vector<GeometricalDomain*>& domains, std::vector<BoundCond*>& boundCond);

	//Accesseur
	std::string getMain() const { return m_nameMain; };
	std::string getMesh() const { return m_nameMesh; };
	std::string getCI() const { return m_nameCI; };
	std::string getModel() const { return m_nameModel; };
  Run *getRun() const { return m_run; }

private:
  Run *m_run;    //pointeur vers run

	int m_vMain;		//!< Number de version file entree main
	int m_vMesh;    //!< Number de version file entree mesh
	int m_vCI;          //!< Number de version file entree initialConditions
	int m_vModel;      //!< Number de version file entree model

	std::string m_nameMain;		//!< Name du file main
	std::string m_nameMesh;   //!< Name du file mesh
	std::string m_nameCI;         //!< Name du file initialConditions
	std::string m_nameModel;     //!< Name du file model
};

#endif // INPUT_H