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

#ifndef TOOLS_H
#define TOOLS_H

#include <string>
#include <cmath>
#include <vector>
#include <list>
#include "Eos/Eos.h"

//#define DEBUG

//! \brief     Enumeration for the primitive-variable type (usefull for second order, slopes, etc.)
enum Prim { vecPhases, vecPhasesO2, vecSlopes, restart };
//! \brief     Enumeration for the axes (X, Y and Z for the axes in the x-, y- and z-direction)
enum Axis { X, Y, Z };
//! \brief     Enumeration for the type of mesh (REC: rectilinear, UNS: unstructured, AMR: adaptative mesh refinement)
enum TypeM { REC, UNS, AMR };
//! \brief     Enumeration for the type of data (FLOAT, DOUBLE, INT, CHAR)
enum TypeData { FLOAT, DOUBLE, INT, CHAR };
//! \brief     Enumeration for the type of geometric object (VERTEX, LINE, PLAN)
enum TypeGO { VERTEX, LINE, PLAN };
//! \brief     Enumeration for the type of boundary (INJ, NONREFLECTING, OUTFLOW, SUBINJ, SYMMETRY, TANK, WALL)
enum TypeBC { INJ = 4, NONREFLECTING = 1, OUTFLOW = 3, SUBINJ = 7, SYMMETRY = 6, TANK = 5, WALL = 2 };
//! \brief     Enumeration for the heat type of wall boundary (ADIABATIC, IMPOSEDTEMP, IMPOSEDFLUX)
enum TypeBCHeat { ADIABATIC = 0, IMPOSEDTEMP = 1, IMPOSEDFLUX = 2};
//! \brief     Enumeration for the type of relaxation (P, PT, PTMu)
enum TypeRelax { P = 1, PT = 2, PTMU = 3 };
//! \brief     Template for the type of the mesh container (std::list for now, but may change to something else if wanted)
template<class Type>
using TypeMeshContainer=std::vector<Type>;

//! \class     Tools
//! \brief     Class for tools
class Tools
{
  public:
    //! \brief     Generic model constructor
    //! \param     numberPhases         number of phases
    Tools(const int& numberPhases);
    ~Tools();

    //! \brief     Modify the string of characters to uppercase it
    //! \param     string               string of characters
    static void uppercase(std::string& string);
    //! \brief     Modify the string of characters to lowercase it
    //! \param     string               string of characters
	  static void lowercase(std::string& string);

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
    double* Deltapk;           //!< Pressure differences, one for each phase
    double* zk;                //!< Acoustic impedance of each phase
    Eos** eos;

    double* Hk0;
    double* Yk0;

};

extern Tools *TB;

#endif // TOOLS_H