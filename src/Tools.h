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

#include "Eos/Eos.h"

//#define DEBUG

//! \brief     Enumeration for the primitive-variable type (usefull for second order, slopes, etc.)
enum Prim
{
  vecPhases,
  vecPhasesO2,
  vecSlopes,
  resume
};

//! \brief     Enumeration for the axes (X, Y and Z for the axes in the x-, y- and z-direction)
enum Axis
{
  X,
  Y,
  Z
};

//! \brief     Enumeration for the type of output (GNU, VTK)
enum TypeOutput
{
  GNU,
  VTK
};

//! \brief     Enumeration for the type of mesh (REC: rectilinear, UNS: unstructured, AMR: adaptative mesh refinement)
enum TypeM
{
  REC,
  UNS,
  AMR
};

//! \brief     Enumeration for the type of data (FLOAT, DOUBLE, INT, CHAR)
enum TypeData
{
  FLOAT,
  DOUBLE,
  INT,
  CHAR
};

//! \brief     Enumeration for the type of geometric object (VERTEX, LINE, PLAN)
enum TypeGO
{
  VERTEX,
  LINE,
  PLAN
};

//! \brief     Enumeration for the type of boundary
enum TypeBC
{
  NONREFLECTING     = 1,
  WALL              = 2,
  OUTLETPRESSURE    = 3,
  INLETINJSTAGSTATE = 4,
  INLETTANK         = 5,
  SYMMETRY          = 6,
  INLETINJTEMP      = 7,
  NULLFLUX          = 8,
  OUTLETMASSFLOW    = 9,
  PISTON            = 10
};

//! \brief     Enumeration for the variable to extract on a boundary
enum VarBoundary
{
  p,
  rho,
  velU,
  velV,
  velW,
  SIZE = (velW + 1)
};

//! \brief     Enumeration for the heat type of wall boundary (ADIABATIC, IMPOSEDTEMP, IMPOSEDFLUX)
enum TypeBCHeat
{
  ADIABATIC   = 0,
  IMPOSEDTEMP = 1,
  IMPOSEDFLUX = 2
};

//! \brief     Enumeration for the type of relaxation (P, PT, PTMu)
enum TypeRelax
{
  U    = 0,
  P    = 1,
  PT   = 2,
  PTMU = 3
};

//! \brief     Enumeration for the gradient method (Finite-Difference (FD), Green-Gauss (GG))
enum TypeGrad
{
  FD,
  GG
};

//! \brief     Enumeration for the phase index of a liquid/vapor couple
enum PhaseType
{
  LIQ = 0,
  VAP = 1
};

//! \brief     Enumeration for the type of second order limiter (NS only)
enum LimiterType
{
  NONE     = 0,
  MINMOD   = 1,
  SUPERBEE = 2
};

//! \brief     Enumeration for the flow variables
enum class Variable
{
  transport,
  pressure,
  density,
  height,
  alpha,
  velocityMag,
  velocityU,
  velocityV,
  velocityW,
  temperature,
  QPA,
  lambda,
  cobaseXX,
  cobaseXY,
  cobaseXZ,
  cobaseYX,
  cobaseYY,
  cobaseYZ,
  cobaseZX,
  cobaseZY,
  cobaseZZ
};

//! \brief     Template for the type of the mesh container (std::list for now, but may change to something else if wanted)
template <class Type> using TypeMeshContainer = std::vector<Type>;

//! \brief     Template for type-safe sign function
template <typename T> int sign(T val) { return (val < T(0)) ? -1 : 1; } //This sign function can take the value -1 or 1
// template <typename T> int sign(T val) { return (T(0) < val) - (val < T(0)); } //This sign function can take the value -1, 1 or 0

//! \brief     Template for deallocation with pointer zeroing
template <class T> inline void destroy(T*& p)
{
  delete p;
  p = nullptr;
}

//! \brief     Template for deallocation with pointer zeroing
template <class T> inline void destroy_array(T*& p)
{
  delete[] p;
  p = nullptr;
}

//! \class     Tools
//! \brief     Class for tools
class Tools
{
  public:
    //! \brief     Generic model constructor
    //! \param     numbPhases         number of phases
    //! \param     numbSolids         number of solid phases
    //! \param     numberTransports   number of additional transport equations
    Tools(const int& numbPhases, const int& numbSolids, const int& numbTransports);
    ~Tools();

    //! \brief     Modify the string of characters to uppercase it
    //! \param     string               string of characters
    static void uppercase(std::string& string);
    //! \brief     Modify the string of characters to lowercase it
    //! \param     string               string of characters
    static void lowercase(std::string& string);

    //! \brief     Swap two numbers
    //! \param     double   1st number to swap
    //! \param     double   2nd number to swap
    static void swap(double& a, double& b);

    //! \brief     Return a non-zero value of a float
    //! \param     double   initial float
    double returnNonZeroValue(double a);

    std::vector<double> ak;
    std::vector<double> Yk;
    std::vector<double> rhok;
    std::vector<double> pk;
    std::vector<double> ek;
    std::vector<double> Ek;
    std::vector<double> akS;
    std::vector<double> rhokS;
    std::vector<double> rhokStar;
    std::vector<double> pkStar;
    std::vector<double> ekStar;
    std::vector<double> EkStar;
    std::vector<double> vkStar;
    std::vector<double> YkStar;
    std::vector<double> Deltapk;       //!< Pressure differences, one for each phase
    std::vector<double> zk;            //!< Acoustic impedance of each phase
    std::vector<double> rho_cIksquare; //!< Density times interface sound speed square, one for each phase
    std::vector<Eos*> eos;
    std::vector<double> Hk0;
    std::vector<double> Yk0;
    std::vector<double> compactionPk_dkappadxi;    //!< Compaction pressure including the term d kappa_k / d xi_k, one for each phase
    std::vector<double> compactionPk_dkappadalpha; //!< Compaction pressure including the term d kappa_k / d alpha_k, one for each phase
    std::vector<double> dlambda;                   //!< Plastic compaction term, one for each phase
    std::vector<Tensor> dplast;                    //!< Plastic shear term, one for each phase
    std::vector<Tensor> LogCobase;                 //!< Log of the cobase
    std::vector<double> LogLambda;                 //!< Log of lambda
    std::vector<bool> alphaNull;                   //!< Parameter to know if we consider the volume fraction as null or not
    std::vector<bool> relaxSolidPlast;             //!< Parameter to know if we need the solid plastic relaxation

    //Variables for RK45 solver (visco-plastic relaxation)
    std::vector<double> ak4;
    std::vector<double> ak5;
    std::vector<double> ek4;
    std::vector<double> ek5;
    std::vector<Tensor> LogCobase4; //!< Log of the cobase
    std::vector<Tensor> LogCobase5; //!< Log of the cobase
    std::vector<double> LogLambda4; //!< Log of lambda
    std::vector<double> LogLambda5; //!< Log of lambda
    std::vector<Tensor> dLogCobaseDt1;
    std::vector<Tensor> dLogCobaseDt2;
    std::vector<Tensor> dLogCobaseDt3;
    std::vector<Tensor> dLogCobaseDt4;
    std::vector<Tensor> dLogCobaseDt5;
    std::vector<Tensor> dLogCobaseDt6;
    std::vector<double> dLogLambdaDt1;
    std::vector<double> dLogLambdaDt2;
    std::vector<double> dLogLambdaDt3;
    std::vector<double> dLogLambdaDt4;
    std::vector<double> dLogLambdaDt5;
    std::vector<double> dLogLambdaDt6;
    std::vector<double> dakDt1;
    std::vector<double> dakDt2;
    std::vector<double> dakDt3;
    std::vector<double> dakDt4;
    std::vector<double> dakDt5;
    std::vector<double> dakDt6;
    std::vector<double> dekDt1;
    std::vector<double> dekDt2;
    std::vector<double> dekDt3;
    std::vector<double> dekDt4;
    std::vector<double> dekDt5;
    std::vector<double> dekDt6;

    static double uselessDouble;
    double physicalTime; //!< Current physical time
};

extern Tools* TB;
extern int numberPhases;
extern int numberSolids;
extern int numberTransports;

#endif // TOOLS_H
