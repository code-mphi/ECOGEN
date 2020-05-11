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

#ifndef OUTPUTPROBEGNU_H
#define OUTPUTPROBEGNU_H

//! \file      OutputProbeGNU.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "OutputGNU.h"
#include "../Maths/GOLine.h"
#include "../Maths/GOPlan.h"

//! \class     OutputProbeGNU
//! \brief     Printing results for a probe in flow

class OutputProbeGNU : public OutputGNU
{
public:
  OutputProbeGNU();
  //! \brief     Probe output constructor from a XML format reading
  //! \details   Reading data from XML file under the following format:
  //!            ex: 	<probe name="capteur1">
  //!                   <vertex x = "0.3" y = "0.05" z = "0.05" / >
  //!                   <timeControl acqFreq = "1e-5." / >       <!-- if negative or nul, recording at each time step-->
  //!                 </probe>
  //! \param     casTest           Folder name of test case input files
  //! \param     run               Resutls folder name (defined in 'mainVX.xml')
  //! \param     element           XML element to read for probe data
  //! \param     fileName          string name of readed XML file
  //! \param     entree            Pointer to corresponding run entry object
  OutputProbeGNU(std::string casTest, std::string run, tinyxml2::XMLElement *element, std::string fileName, Input *entree);
  virtual ~OutputProbeGNU();

  virtual void locateProbeInMesh(const TypeMeshContainer<Cell *> &cells, const int &nbCells, bool localSeeking = false);
  virtual Cell* locateProbeInAMRSubMesh(std::vector<Cell*>* cells, const int &nbCells);

  virtual void prepareSortieSpecifique();
  virtual void ecritSolution(Mesh *mesh, std::vector<Cell *> *cellsLvl);

  virtual void prepareOutputInfos() {}; //nothing to print
  virtual void ecritInfos() {};

  //Accessors
  virtual double getNextTime() { return m_nextAcq; };
  virtual bool possesses() { return m_possessesProbe; };


private:
  double m_acqFreq;           //!< Acquisition time frequency
  double m_nextAcq;           //!< Next acquisition time
  Cell *m_cell;               //!< Pointer to the level 0 cell containing the probe
  GeometricObject *m_objet;   //!< To store position
  bool m_possessesProbe;      //!< True if the CPU possesses probe
};

#endif //OUTPUTPROBEGNU_H