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

#ifndef STRETCHZONE_H
#define STRETCHZONE_H

#include <vector>
#include <string>

//! \class     stretchZone
//! \brief     managing a stretched zone in Cartesian meshes

class stretchZone
{
  public:
    //! \brief  Stretched zone constructor
    //
    //! Stretched zone constructor: compute number of cells from mesh stretching area size,
    //          size of reference cell and progression factor applied to reference cell.
    //          Store data needed to apply stretching later.
    //! \param  startAt       position of beginning of mesh stretching area
    //! \param  endAt         position of end of mesh stretching area
    //! \param  factor        mesh stretching factor: 1. means no stretching, value lower than 1. means mesh shrinking
    //                        and value greater than 1. mesh coarsening.
    //! \param  refCellSize   size of reference cell to which mesh stretching factor is applied
    //! \param  refCellPos    position of the reference cell ('end' or 'start')
    stretchZone(double startAt, double endtAt, double factor, double refCellSize, std::string refCellPos);
    virtual ~stretchZone();

    int stretching(std::vector<double>& dX, std::vector<double>& posX);

    static int verifyStretching(std::vector<stretchZone>& tabStretch, const double l, std::string fileName = "");

    // Accessors
    const int& getNumberCells() const { return m_numberCells; };

  private:
    double m_startAt;         //!< zone starting position along corresponding axis
    double m_endAt;           //!< zone ending position along corresponding axis
    double m_factor;          //!< stretching factor (1. is no stretching)
    double m_refCellSize;     //!< size of reference cell
    int m_numberCells;        //!< number of cells in the stretched zone
    std::string m_refCellPos; //!< position of reference cell (start or end)

    int32_t safe_round_to_int32(double x) const;
    void factorCorrection(double l);
};

#endif // STRETCHZONE_H
