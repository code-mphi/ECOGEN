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

#ifndef OutputBoundaryFluxGNU_H
#define OutputBoundaryFluxGNU_H

#include "OutputGNU.h"

class OutputBoundaryFluxGNU : public OutputGNU
{
public:
  OutputBoundaryFluxGNU(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input* entree);
  virtual ~OutputBoundaryFluxGNU();

  virtual void prepareSortieSpecifique(std::vector<CellInterface*>* cellInterfacesLvl);
  virtual void ecritSolution(std::vector<CellInterface*>* cellInterfacesLvl);
  virtual void extractFlux(std::vector<CellInterface*>* cellInterfacesLvl);

  //Accessors
  virtual double getNextTime() { return m_nextAcq; };

protected:
  std::string m_fluxType; //!< Flux type could be either massflow or power flux
  double m_flux;          //!< Flux recorded through boundary either massflow (kg.s-1) or power flux (W)
  int m_numPhys;          //!< Physical number of the boundary to record
  std::vector<int> cellInterfaceIndexes; //!< Indexes of cellInterfaces on the boundary (speed up searching process)

  double m_acqFreq; //!< Acquisition time frequency 
  double m_nextAcq; //!< Next acquisition time
};

#endif // OutputBoundaryFluxGNU_H