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

#ifndef PHASENONLINEARSCHRODINGER_H
#define PHASENONLINEARSCHRODINGER_H

#include "../EulerKorteweg/PhaseEulerKorteweg.h"
#include "../../Eos/Eos.h"
#include <fstream>

//! \class     PhaseNonLinearSchrodinger
//! \brief     Phase variables for Augmented Euler--Korteweg equations (single phase)
class PhaseNonLinearSchrodinger : public PhaseEulerKorteweg
{
  public:
    PhaseNonLinearSchrodinger();
    //! \brief     Phase constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!            ex: <dataFluid density = "10.0">
    //!                  <velocity x = "1000." y = "1000." z = "0." / >
    //!                </dataFluid>
    //! \param     material           XML element to read for phase data
    //! \param     eos                EOS pointer to compute thermodynamic variables
    //! \param     fileName           string name of readed XML file
    PhaseNonLinearSchrodinger(tinyxml2::XMLElement* material, Eos* eos, std::string fileName);
    virtual ~PhaseNonLinearSchrodinger();

    virtual void allocateAndCopyPhase(Phase** vecPhase);

    //Specific methods for data printing
    //----------------------------------
    virtual int getNumberScalars() const { return 3; };

    //Specific methods for parallel computing
    //---------------------------------------
    virtual int numberOfTransmittedVariables() const;
    virtual void fillBuffer(double* buffer, int& counter) const;
    virtual void fillBuffer(std::vector<double>& dataToSend) const;
    virtual void getBuffer(double* buffer, int& counter, Eos** /*eos*/);
    virtual void getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** /*eos*/);

    //Verifications
    //-------------
    virtual void verifyAndCorrectDensityMax() {};
};

#endif // PHASENONLINEARSCHRODINGER_H
