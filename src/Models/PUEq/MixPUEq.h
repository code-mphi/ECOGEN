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

#ifndef MIXPUEQ_H
#define MIXPUEQ_H

#include <vector>
#include "../UEq/MixUEq.h"

//! \class     MixPUEq
//! \brief     Mixture variables for the pressure-velocity-equilibrium (mechanical equilibrium) system of equations (Kapila)
class MixPUEq : public MixUEq
{
    public:
      MixPUEq();
      //! \brief     Mixture constructor from a XML format reading
      //! \details   Reading data from XML file under the following format:
      //!           ex: <mixture>
      //!                 <velocity x = "0." y = "0." z = "0." />
      //!               </mixture>
      //! \param     state           XML element to read for mixture data
      //! \param     fileName       string name of readed XML file
      MixPUEq(tinyxml2::XMLElement* state, std::string fileName);
      virtual ~MixPUEq();

      virtual void allocateAndCopyMixture(Mixture** mixture);
};

#endif // MIXPUEQ_H
