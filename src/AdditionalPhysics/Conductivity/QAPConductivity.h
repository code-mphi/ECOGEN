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

#ifndef QAPCONDUCTIVITY_H
#define QAPCONDUCTIVITY_H

//! \file      QAPConductivity.h
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "../QuantitiesAddPhys.h"

//! \class     QAPConductivity
//! \brief     General class for thermal conductive quantities
class QAPConductivity : public QuantitiesAddPhys
{
    public:
    QAPConductivity();
    QAPConductivity(AddPhys* addPhys, const int &numberPhases);
    virtual ~QAPConductivity();

    virtual void computeQuantities(Cell* cell);

    //Accessors
    virtual void setGrad(const Coord &grad, int num = -1);
    virtual Coord getGrad(int num = -1) const;

    protected:
    Coord *m_gradTk;          //!< Vector of the temperature gradient of each phase of the cell

    private:
};

#endif // QAPCONDUCTIVITY_H
