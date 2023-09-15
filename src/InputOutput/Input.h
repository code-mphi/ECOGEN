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
#include "../Gradients/HeaderGradient.h"
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

    void readInputXML(std::vector<GeometricalDomain*>& domains, std::vector<BoundCond*>& boundCond, std::vector<GeometricalDomain*>& solidDomains);

    void inputMain(std::string testCase);
    void inputMesh(std::string testCase);
    void inputModel(std::string testCase);
    Eos* inputEOS(std::string EOS, int& numberEOS);
    void inputInitialConditions(std::string testCase, std::vector<GeometricalDomain*>& domains, std::vector<BoundCond*>& boundCond, std::vector<GeometricalDomain*>& solidDomains);
    void verifyCompatibilityInput(std::string testCase);

	//Accesseur
	std::string getMain() const { return m_nameMain; };
	std::string getMesh() const { return m_nameMesh; };
	std::string getCI() const { return m_nameCI; };
	std::string getModel() const { return m_nameModel; };
  Run *getRun() const { return m_run; }

private:
  Run *m_run; // Pointer to run

	std::string m_nameMain;   //!< Name of main file
	std::string m_nameMesh;   //!< Name of mesh file
	std::string m_nameCI;     //!< Name of initialConditions file
	std::string m_nameModel;  //!< Name of model file
};

#endif // INPUT_H