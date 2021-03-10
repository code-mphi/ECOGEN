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

#ifndef QAPVISCOSITY_H
#define QAPVISCOSITY_H

#include "../QuantitiesAddPhys.h"

//! \class     QAPViscosity
//! \brief     General class for viscous quantities
class QAPViscosity : public QuantitiesAddPhys
{
    public:
    QAPViscosity(AddPhys* addPhys);
    virtual ~QAPViscosity();

    virtual void computeQuantities(Cell* cell);

    //Accessors
    virtual void setGrad(const Coord& grad, const int& num = -1);                       //1:U, 2:V, 3:W
    virtual const Coord& getGrad(const int& num = -1) const { return m_grads[num-1]; }; //1:U, 2:V, 3:W

    protected:
    std::vector<Coord> m_grads;    //!< Gradient vectors of the velocities of the cell in x-, y- and z-directions (ex: m_grads[0] = (du/dx, du/dy, du/dz))

    private:
};

#endif // QAPVISCOSITY_H
