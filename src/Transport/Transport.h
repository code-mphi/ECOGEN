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

#ifndef TRANSPORT_H
#define TRANSPORT_H

//! \file      Transport.h
//! \author    K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "../libTierces/tinyxml2.h"
#include <fstream>

//! \class     Transport
//! \brief     Class for additional transport equations
//! \details   This class is used in particular for surface-tension effects computation
class Transport
{
  public:
    Transport();
    ~Transport();

    //! \brief     Set the value of the corresponding transport variable
    //! \param     value                  value of the transport variable
    void setValue(double value);
    //! \brief     Return the value of the corresponding transport variable
    const double& getValue() const { return m_value; };

    //! \brief     Cell to cell Riemann solver for the corresponding transport equation
    //! \param     transportLeft          left value of transport variable
    //! \param     transportRight         right value of transport variable
    //! \param     sM                     fluid velocity for intercell interfaces
   	void solveRiemann(double transportLeft, double transportRight, double sM);
    //! \brief     Wall half Riemann solver for the corresponding transport equation
    void solveRiemannWall();
    //! \brief     Inflow (injection) half Riemann solver for the corresponding transport equation
    //! \param     transportLeft          left value of transport variable
    //! \param     sM                     fluid velocity for intercell interfaces
    //! \param     valueTransport         injection value of transport variable
    void solveRiemannInflow(double transportLeft, double sM, double valueTransport);
    //! \brief     Tank half Riemann solver for the corresponding transport equation
    //! \param     transportLeft          left value of transport variable
    //! \param     sM                     fluid velocity for intercell interfaces
    //! \param     valueTransport         tank value of transport variable
    void solveRiemannTank(double transportLeft, double sM, double valueTransport);
    //! \brief     Outflow half Riemann solver for the corresponding transport equation
    //! \param     transportLeft          left value of transport variable
    //! \param     sM                     fluid velocity for intercell interfaces
    //! \param     valueTransport         outflow value of transport variable
    void solveRiemannOutflow(double transportLeft, double sM, double valueTransport);
    //! \brief     Add flux to the corresponding transport buffer flux
    //! \param     coefA                  possibility to multiply the flux before adding (set 1.d0 if not needed)
    //! \param     num                    number of the corresponding transport equation
    void addFlux(double coefA, const int num);
    //! \brief     Subtract flux to the corresponding transport buffer flux
    //! \param     coefA                  possibility to multiply the flux before adding (set 1.d0 if not needed)
    //! \param     num                    number of the corresponding transport equation
    void subtractFlux(double coefA, const int num);
    //! \brief     Add non conservative transport term to the flux
    //! \param     coefA                  possibility to multiply the non conservative transport term before adding (set 1.d0 if not needed)
    //! \param     transport              transport value used to approximate the non conservative transport term
    //! \param     sM                     fluid velocity for intercell interfaces
    void addNonCons(double coefA, double transport, const double sM);
    //! \brief     Subtract non conservative transport term to the flux
    //! \param     coefA                  possibility to multiply the non conservative transport term before adding (set 1.d0 if not needed)
    //! \param     transport              transport value used to approximate the non conservative transport term
    //! \param     sM                     fluid velocity for intercell interfaces
    void subtractNonCons(double coefA, double transport, const double sM);
    //! \brief     Multiply the corresponding transport value by a scalar
    //! \param     scalar                 scalar by what the transport value is multiplied
    void multiply(double scalar);
    //! \brief     Add a scalar to the corresponding transport value
    //! \param     scalar                 added scalar to the transport value
    void add(double scalar);
    //! \brief     Change the sign of the corresponding transport value
    void changeSign();
    
    //Specifique secondOrder
    //! \brief     Compute the slope at the edge of a cell
    //! \param     valueLeft              transport value of the left cell
    //! \param     valueRight             transport value of the right cell
    //! \param     distance               distance between the left and right cell centers
    void computeSlopeTransport(const double valueLeft, const double valueRight, const double &distance);
    //! \brief     Extrapolate the value of the corresponding transport equation from the center of the cell to its edge
    //! \param     slope                  value of the slope
    //! \param     distance               distance between the center and the corresponding edge of the cell
    void extrapolate(const double &slope, const double &distance);

  private:
    double m_value;     //! Value of the corresponding transport variable
};

extern Transport* fluxBufferTransport;

#endif // TRANSPORT_H