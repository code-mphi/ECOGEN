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

#ifndef CELLGHOST_H
#define CELLGHOST_H

//! \file      CellGhost.h
//! \author    K. Schmidmayer, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include "Cell.h"

//! \class     CellGhost
//! \brief     Child class for a ghost mesh cell
class CellGhost : public Cell
{
    public:
        //! \brief     Ghost cell constructor for a non AMR cell
        CellGhost();
        //! \brief     Ghost cell constructor for an AMR cell
        //! \param     lvl    level of current AMR cell
        CellGhost(int lvl); //Pour AMR
        virtual ~CellGhost();

        virtual int getRankOfNeighborCPU() const;
        virtual void setRankOfNeighborCPU(int rank);
        virtual void createChildCell(const int &lvl);
        virtual bool isCellGhost() const { return true; };

    protected:
      int m_rankOfNeighborCPU; /*!< Rank of the neighbor CPU corresponding to this ghost cell */

    private:
};

#endif // CELLGHOST_H
