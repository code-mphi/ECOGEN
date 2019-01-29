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

#ifndef TOOLS_H
#define TOOLS_H

//! \file      Tools.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      December 6 2018

#include <string>
#include <cmath>
#include "Eos/Eos.h"

//! \brief     Enumeration for the axe (X, Y and Z for the axe in the x-, y- and z-direction)
typedef enum Axe{ X, Y, Z } Axe;
//! \brief     Enumeration for the type of mesh (REC: rectilinear, UNS: unstructured, AMR: adaptative mesh refinement)
typedef enum TypeM { REC, UNS, AMR } TypeM;
//! \brief     Enumeration for the type of data (FLOAT, DOUBLE, INT, CHAR)
typedef enum TypeData { FLOAT, DOUBLE, INT, CHAR } TypeData;
//! \brief     Enumeration for the type of geometric object (VERTEX, LINE, PLAN)
typedef enum TypeGO { VERTEX, LINE, PLAN } TypeGO;

//! \class     Tools
//! \brief     Class for tools
class Tools
{
  public:
    Tools();
    //! \brief     Generic model constructor
    //! \param     numberPhases         number of phases
    Tools(const int &numberPhases);
    virtual ~Tools();

    //! \brief     Modify the string of characters to uppercase it
    //! \param     string               string of characters
    static void uppercase(std::string &string);
    //! \brief     Return the value of pi
    static double pi();

    double m_numberPhases;
    double* ak;
    double* Yk;
    double* rhok;
    double* pk;
    double* akS;
    double* rhokS;
    double* rhokStar;
    double* pkStar;
    double* ekStar;
    double* vkStar;
    double* YkStar;
    double* Deltaek;
    Eos** eos;

    double* Hk0;
    double* Yk0;

};

extern Tools *TB;

#endif // TOOLS_H