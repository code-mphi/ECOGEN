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

#ifndef RUN_H
#define RUN_H

//! \file      Run.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      June 27 2018

class Run;

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <sstream>
#include "Cell.h"
#include "Models/HeaderPhase.h"
#include "CellInterface.h"
#include "Parallel.h"
#include "Meshes/HeaderMesh.h"
#include "BoundConds/HeaderBoundCond.h"
#include "Eos/HeaderEquationOfState.h"
#include "Models/HeaderModel.h"
#include "Geometries/HeaderGeometricalDomain.h"
#include "Ordre2/HeaderLimiter.h"
#include "AdditionalPhysics/HeaderQuantitiesAddPhys.h"
#include "AdditionalPhysics/HeaderAddPhys.h"

#include "InputOutput/Input.h"
#include "InputOutput/Output.h"
#include "timeStats.h"

#include "Relaxations/HeaderRelaxations.h"

//! \class     Run
//! \brief     Class regrouping all information for a simulation
class Run
{
  public:
    Run(std::string nameCasTest, const int &number);
    virtual ~Run();

    //! \brief    Initialization of the simulation
    void initialize(int argc, char* argv[]);
    //! \brief    Hyperbolic resolution + Relaxations + Source terms integration
    //! \details  Hyperbolic part is solved using Finite Volume method : \f[ \frac{U^{n+1}_i-U^{n}_i}{\Delta t} =  -\sum_{faces} \vec{F}^*_f \cdot \vec{n}_f \f]
    void solver();
    //! \brief    Cleaning simulation
    //! \details  Memory desallocations
    void finalize();
    
    void resumeSimulation(int &iteration, double &dt, double &tempsPhysique);

    //Accessors
    int getNumberPhases() const;

  private:   

    //Specific solvers
    void integrationProcedure(double &dt, int lvl, double &dtMax, int &nbCellsTotalAMR);
    void advancingProcedure(double &dt, int &lvl, double &dtMax) const;
    void solveHyperbolic(double &dt, int &lvl, double &dtMax) const;
    void solveHyperbolicO2(double &dt, int &lvl, double &dtMax) const;
    void solveAdditionalPhysics(double &dt, int &lvl) const;
    void solveSourceTerms(double &dt, int &lvl) const;
    void solveRelaxations(int &lvl) const;
    void verifyErrors() const;

    int m_numTest;                             //!<Number of the simulation

    //Input attributes
    std::string m_simulationName;              //!<Name of the simulation
    bool m_controleIterations;                 //!Choice for time control mode (iteration or physical time)
    int m_nbIte, m_freq;                       //!<Requested number of final time iteration and frequency
    float m_finalPhysicalTime, m_timeFreq;     //!<Requested final physical time of the simulation and time frequency for output printing
    double m_cfl;                              //!<CFL criteria (between 0 and 1)
    int m_numberPhases;                        //!<Number of phases
    int m_numberEos;                           //!<Number of equations of states
    int m_numberTransports;                    //!<Number of additional transport variables
    int m_numberAddPhys;                       //!<Number of additional physical effects
    int m_numberSources;                       //!<Number of additional source terms
    int m_dimension;                           //!<dimension 1, 2 ou 3
    int m_MRF;                                 //!<source term for Moving Reference Frame computation index(in the list of source term)
    std::string m_order;                       //!<Precision scheme order (firstorder or secondOrder)

    //Specific to AMR method
    int m_lvlMax;                              //!<Maximum AMR level (if 0, then no AMR)
    int m_nbCellsTotalAMR;                     //!<Number de mailles total maximum durant la simulation
    std::vector<Cell *> *m_cellsLvl;           //!<Tableau de vecteurs contenant les cells de compute, un vecteur par niveau.
    std::vector<CellInterface *> *m_boundariesLvl;   //!<Tableau de vecteurs contenant les boundaries de compute, un vecteur par niveau.

    //Geometrical attributes
    bool m_parallelPreTreatment;               //!<Choice for mesh parallel pre-treatment  (needed for first simulation on a new parallel unstructured geometry)
    
    //Calcul attributes
    Mesh *m_mesh;                              //!<Mesh type object: contains all geometrical properties of the simulation
    Model *m_model;                            //!<Model type object: contains the flow model methods
    Cell **m_cells;                            //!<Array of computational cells objects: contains physical fluids states
    CellInterface **m_boundaries;              //!<Array of interfaces objects between cells (or between a cell and a boundary)
    Eos **m_eos;                               //!<Array of Equations of states: contains fluids EOS parameters
    std::vector<AddPhys*> m_addPhys;           //!<Vector of Additional physics
    Symmetry *m_symmetry;                      //!<Specific object for symmetry (cylindrical or spherical) if active
    Symmetry *m_symmetryAddPhys;               //!<Object containing the parent class of symmetry to trick the corresponding additional physics argument (avoid taking into account symmetry terms multipled times)
    std::vector<Source*> m_sources;            //!<Vector of source terms
    Limiter *m_globalLimiter;                  //!<Slope limiter type object for second order in space
    Limiter *m_interfaceLimiter;               //!<Slope limiter type object for second order in space specific to interface location
    Limiter *m_globalVolumeFractionLimiter;    //!<Slope limiter type object for second order in space specific to interface advected variables (alpha, transports)
    Limiter *m_interfaceVolumeFractionLimiter; //!<Slope limiter type object for second order in space specific to interface advected variables (alpha, transports) and to interface location
    std::vector<std::string> m_nameGTR;        //!<Vector of transport variable names
    std::vector<std::string> m_nameQPA;        //!<Vector of names of the quantities of additional physics
    std::vector<std::string> m_nameGPH;        //!<Vector of phasic variables name
    double m_dt;                               //!<Explicit time step
    double m_physicalTime;                     //!<Physical time
    int m_iteration;                           //!<time iteration number
    int m_resumeSimulation;                    //!<File number for restarting a simulation

    //Input/Output attributes
	  Input* m_input;						                 //!<Input object
    Output* m_outPut;                          //!<Main output object
    std::vector<Output *> m_cuts;              //!<Vector of output objects for cuts
    std::vector<Output *> m_probes;            //!<Vector of output objects for probes
    timeStats m_stat;                          //!<Object linked to computational time statistics
    double *m_pMax, *m_pMaxWall;             //!<Maximal pressure found between each written output and its corresponding coordinate (only for few test case)

    friend class Input;
    friend class Output;
    friend class OutputXML;
    friend class OutputGNU;
    friend class OutputProbeGNU;
    friend class Mesh;
};

#endif // RUN_H
