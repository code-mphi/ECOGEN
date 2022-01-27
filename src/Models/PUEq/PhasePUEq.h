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

#ifndef PHASEPUEQ_H
#define PHASEPUEQ_H

#include "../UEq/PhaseUEq.h"
#include "../../Eos/Eos.h"
#include <fstream>

//! \class     PhasePUEq
//! \brief     Phase variables for the pressure-velocity-equilibrium (mechanical equilibrium) system of equations (Kapila)
class PhasePUEq : public PhaseUEq
{
  public:
    PhasePUEq();
    //! \brief     Phase constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!           ex:  <dataFluid alpha="0.5" density="1.0" pressure="1.e5"/> 
    //! \param     material           XML element to read for phase data
    //! \param     eos                EOS pointer to compute thermodynamic variables
    //! \param     pressure           Pressure from mixture read
    //! \param     fileName           string name of readed XML file
    PhasePUEq(tinyxml2::XMLElement* material, Eos* eos, const double& pressure, std::string fileName);
    virtual ~PhasePUEq();

    virtual void allocateAndCopyPhase(Phase** vecPhase);

    //Specific methods for data printing
    //----------------------------------
    virtual int getNumberScalars() const { return numberScalarsPhase; };
};

#endif // PHASEPUEQ_H
