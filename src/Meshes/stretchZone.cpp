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

#include <sstream>
#include <limits>

#include <cstdint>
#include "stretchZone.h"
#include "../Errors.h"
#include "../Parallel/Parallel.h"
#include "../Colors.h"

//***********************************************************************

stretchZone::stretchZone(double startAt, double endAt, double factor, double refCellSize, std::string refCellPos) :
  m_startAt(startAt), m_endAt(endAt), m_factor(factor), m_refCellSize(refCellSize), m_refCellPos(refCellPos)
{
  double numCell_d = 1.;

  // Check user input
  if (m_refCellSize <= 0.0) {
    throw ErrorXMLStretching("Reference cell size cannot be null or negative.", "", __FILE__, __LINE__);
  }
  if (m_factor <= 0.0) {
    throw ErrorXMLStretching("Stretching factor cannot be null or negative.", "", __FILE__, __LINE__);
  }
  if (m_endAt <= m_startAt) {
    throw ErrorXMLStretching("Ending position of mesh stretching area has to be greater than its starting position.", "", __FILE__, __LINE__);
  }

  double l{m_endAt - m_startAt};
  if (l < m_refCellSize) {
    throw ErrorXMLStretching("Reference cell size cannot be greater than mesh stretching area.", "", __FILE__, __LINE__);
  }

  /// Recall:
  ///   - geometric sequence is defined by:
  ///     + \f$ u_{n+1} = u_n * q \f$
  ///     + \f$ u_n = u_0 * q^n \f$
  ///
  ///   - the sum is given by:  \f$ sum_{k = 0}^{n} u_k = u_0 (1 - q^{n+1}) / (1 - q) \f$
  ///
  ///   - If we note l the length of the stretching area, we want:
  ///     + \f$ u_0 = l = S_0 \f$ when we have only 1 cell in the stretching area.
  ///     + i.e.,  \f$ Nc = n+1 \f$, with Nc the number of cells in the stretching area.
  ///     + and \f$ l = S_{Nc-1} \f$.
  ///
  /// In consequence, we have \f$ Nc = \frac{ ln(1-(1-q)*l/u_0) }{ ln(q) } \f$.
  if (m_factor == 1.) {
    // Constant cell size
    numCell_d = l / m_refCellSize;
  }
  else {
    double tmp = l * (1. - m_factor);

    // Check for cases where ln cannot be computed (the reference cell size and the progression factor will not allow the serie
    // to cover the patch length).
    if (m_refCellSize <= tmp) {
      throw ErrorXMLStretching("Mesh stretching area cannot be covered with this reference cell size and stretching factor.\n"
                               "Increase the reference cell size or the stretching factor.\n",
                               "",
                               __FILE__,
                               __LINE__);
    }

    // If refCellSize is smaller than l,  ln( 1. - tmp / refCellSize) / ln(m_factor) is greater than 1. so numCell_d is strictly greater than 0.
    numCell_d = std::log(1. - tmp / m_refCellSize) / std::log(m_factor);
  }

  try {
    // We want entire cells: compute the wanted cell count
    m_numberCells = safe_round_to_int32(numCell_d);
  }
  catch (ErrorECOGEN&) {
    throw ErrorXMLStretching("Mesh stretching initialization: Computed cell number cannot be casted into an integer.\n"
                             "Increase the reference cell size or the stretching factor.\n",
                             "",
                             __FILE__,
                             __LINE__);
  }

  // Compute suitable progression factor to end with the exact cell count previously computed.
  factorCorrection(l);

  if (rankCpu == 0) {
    std::cout << "     " << m_numberCells << " cell(s) in mesh stretching area [" << startAt << ";" << endAt
              << "] with (possibly modified) progression factor " << m_factor << ".\n";
  }
}

//***********************************************************************

stretchZone::~stretchZone() {}

//***********************************************************************

int32_t stretchZone::safe_round_to_int32(double x) const
{
  if (std::isnan(x)) {
    throw std::domain_error("NaN cannot be converted to int32_t");
  }
  double rounded = std::round(x);

  // INT32_MAX is exactly representable as double since it fit within 53 bits.
  constexpr double kMax = static_cast<double>(INT32_MAX);
  if (rounded > kMax) {
    throw ErrorECOGEN("Double out of int32_t range.");
  }
  return static_cast<int32_t>(rounded);
}

//***********************************************************************

void stretchZone::factorCorrection(double l)
{
  switch (m_numberCells) {
  case 1:
    // If we have 1 cell in stretching zone, we want u_0 == l. Check the normalized distance between l and u_0: we want |l-u_0|/l < \epsilon
    if (fabs(l - m_refCellSize) > std::numeric_limits<float>::epsilon() * l) {
      std::cout << "  Warning: only one cell in mesh stretching area => size of reference cell is redefined to: " << l << "\n";
    }
    m_refCellSize = l;
    m_factor      = 1.;
    break;
  case 2:
    // l = u_0 + u_0*q = (1+q) u_0 => q  = l/u_0 - 1.
    m_factor = l / m_refCellSize - 1.;
    break;
  case 3:
    // l = u_0 + u_0*q + u_0*q^2 => q = 0.5 * (\sqrt{4l/u_0 - 3} -1).
    m_factor = 0.5 * (std::sqrt(4 * l / m_refCellSize - 3) - 1.);
    break;
  default:
    // polynome of arbitrary level: use Newton-Raphson algo to solve it.
    constexpr int maxit{500};
    const double tmp{l / m_refCellSize};
    int it{0};
    double q_Nminus1{pow(m_factor, m_numberCells - 1)};         // q^{n-1}, computed once only for efficiency when compiling with -O0
    double f{m_factor * q_Nminus1 - tmp * m_factor + tmp - 1.}; // f(q) = q^n - q * l/u_0 + l/u0 -1;
    while (f > std::numeric_limits<double>::epsilon() && ++it < maxit) {
      double df{m_numberCells * q_Nminus1 - tmp}; // f'(q) = nq^{n-1} - l/u_0;
      assert(df != 0.0 && "Cannot use Newton algorithm with 0 derivative");
      double factor{m_factor - f / df};
      if (m_factor == factor) { // Stagnation
        break;
      }
      m_factor  = factor;
      q_Nminus1 = pow(m_factor, m_numberCells - 1);
      f         = m_factor * q_Nminus1 - tmp * m_factor + tmp - 1.;
    }
    if (it == maxit) {
      std::stringstream msg;
      msg << "Mesh stretching initialization: fail to compute the progression factor needed to have exactly " << m_numberCells << " cell(s).\n"
          << "  Newton algorithm failed to converge.\n";
      throw ErrorXMLStretching(msg.str(), "", __FILE__, __LINE__);
    }
    break;
  }
  return;
}

//***********************************************************************

int stretchZone::stretching(std::vector<double>& dX, std::vector<double>& posX)
{
  const size_t prevSize{dX.size()};
  const size_t newSize{prevSize + m_numberCells};
  dX.resize(newSize);
  posX.resize(newSize);

  size_t indicePos;
  int dir;
  double startPos, endPos;

  double dXi{m_refCellSize};
  if (m_refCellPos == "START") {
    indicePos = prevSize;
    dir       = 1;
    startPos  = m_startAt;
    endPos    = m_endAt;
  }
  else {
    // refCellPos == "END": apply mesh stretching starting from end of mesh stretching area.
    indicePos = newSize - 1;
    dir       = -1;
    startPos  = m_endAt;
    endPos    = m_startAt;
  }

  dX[indicePos]    = m_refCellSize;
  posX[indicePos]  = startPos + dir * 0.5 * m_refCellSize;
  indicePos       += dir;
  for (int i = 1; i < m_numberCells - 1; i++) {
    dXi             *= m_factor; // u_n = u_{n-1} * q
    dX[indicePos]    = dXi;
    posX[indicePos]  = posX[indicePos - dir] + dir * 0.5 * (dX[indicePos - dir] + dX[indicePos]);
    indicePos       += dir;
  }
  if (m_numberCells > 1) {
    // Compute last cell data from the stretching area boundary in order to avoid issues linked to numerical precision (if cell progression doesn't end up exactly on patch boundary, dX and posX may be false).
    double cellStart = posX[indicePos - dir] + dir * 0.5 * dX[indicePos - dir];
    dX[indicePos]    = dir * (endPos - cellStart);
    posX[indicePos]  = endPos - dir * 0.5 * dX[indicePos];
  }

  return m_numberCells;
}

//***********************************************************************

int stretchZone::verifyStretching(std::vector<stretchZone>& tabStretch, const double l, std::string fileName)
{
  if (tabStretch.size() == 0) return 0;
  try {

    if (tabStretch[0].m_startAt != 0.) throw ErrorXMLStretching("Mesh stretching area has to start at 0.", fileName, __FILE__, __LINE__);
    if (tabStretch[tabStretch.size() - 1].m_endAt != l)
      throw ErrorXMLStretching("Mesh stretching area has to recover the entire domain.", fileName, __FILE__, __LINE__);
    for (unsigned int i = 1; i < tabStretch.size(); i++) {
      if (tabStretch[i].m_startAt != tabStretch[i - 1].m_endAt)
        throw ErrorXMLStretching("Mesh stretching area should not contain holes.", fileName, __FILE__, __LINE__);
      //Verifying size factor between two neighbouring cells
      if (rankCpu == 0) {
        stretchZone szL{tabStretch[i - 1]};
        stretchZone szR{tabStretch[i]};
        double dXnL, dX1R;

        // Assign suitable cell sizes depending on the position of reference cell.
        if (szL.m_refCellPos == "START") {
          dXnL = szL.m_refCellSize * pow(szL.m_factor, szL.m_numberCells - 1);
        }
        else {
          dXnL = szL.m_refCellSize;
        }
        if (szR.m_refCellPos == "START") {
          dX1R = szR.m_refCellSize;
        }
        else {
          dX1R = szR.m_refCellSize * pow(szR.m_factor, szR.m_numberCells - 1);
        }

        if (dXnL > 1.1 * dX1R || dXnL < 0.9 * dX1R) {
          std::cout << ECOGEN_color_magenta;
          std::cout << "     WARNING: Large cell size ratio between zones " << i - 1 << " and " << i << " (ratio: " << dXnL / szR.m_refCellSize
                    << ").\n";
          std::cout << "              It can create numerical issues." << std::endl;
          std::cout << ECOGEN_color_reset;
        }
      }
    }
  }
  catch (ErrorXML&) {
    throw;
  } // Renvoi au niveau suivant
  return 0;
}

//***********************************************************************
