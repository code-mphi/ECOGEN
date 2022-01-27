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

#ifndef OUTPUTBOUNDARYFLUXGNU_H
#define OUTPUTBOUNDARYFLUXGNU_H

#include "OutputBoundaryGNU.h"

enum FluxType { MASSFLOW, POWERFLUX };

class OutputBoundaryFluxGNU : public OutputBoundaryGNU
{
public:
  OutputBoundaryFluxGNU(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input* entree);
  virtual ~OutputBoundaryFluxGNU();

  // Virtual methods
  virtual void initializeSpecificOutputBound();
  virtual void writeResults(std::vector<CellInterface*>* cellInterfacesLvl);

protected:
  //! \brief  Get flux either massflow or enthalpy through the boundary 
  double getFlux(std::vector<CellInterface*>* cellInterfacesLvl);

  //! \brief  Extract massflow throught the whole boundary surface
  double extractMassflow(std::vector<CellInterface*>* cellInterfacesLvl);

  //! \brief  Extract enthalpy flux throught the whole boundary surface
  double extractEnthalpyFlux(std::vector<CellInterface*>* cellInterfacesLvl);
  
  //! \brief  Compute the massflow contribution of a single cell interface
  double computeMassflowFace(CellInterface *bound);

  //! \brief  Compute the enthalpy flux contribution of a single cell interface
  double computeTotalEnthalpyFluxFace(CellInterface *bound);

  //! \brief  Compute the enthalpy flux contribution of a single cell interface when MRF is activated
  double computeTotalEnthalpyFluxFaceMRF(CellInterface *bound);

  FluxType m_fluxType;    //!< Flux type could be either massflow or power flux
  double m_flux;          //!< Flux recorded through boundary either massflow (kg.s-1) or power flux (W)
};

#endif // OUTPUTBOUNDARYFLUXGNU_H