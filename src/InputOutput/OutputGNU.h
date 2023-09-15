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

#ifndef OUTPUTGNU_H
#define OUTPUTGNU_H

#include "Output.h"
class OutputGNU : public Output
{
public:
  
  //! \brief   Default constructor for specific output without specific needs
  OutputGNU();
  
  //! \brief   Main constructor for datasets
  //! \param   casTest   Test case name (defined in "main.xml")  
  //! \param   nameRun   Folder to store results
  //! \param   element   XML outputMode element
  //! \param   fileName  Full path to main.xml of current test case
  //! \param   entree    Input pointer to access run pointer and its information
  OutputGNU(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input *entree);
  
  //! \brief   Constructor for specific derived GNU outputs (boundary, probe, cut)
  //! \param   element   XML GNU output element to get stream precision
  OutputGNU(tinyxml2::XMLElement* element);
  
  virtual ~OutputGNU();

  virtual void initializeSpecificOutput() {};  // Nothing to do for this output
  virtual void initializeSpecificOutput(std::vector<CellInterface*>* /*cellInterfacesLvl*/) {}; // Nothing to do for this output
  virtual void writeResults(Mesh *mesh, std::vector<Cell*>* cellsLvl);
  virtual void writeResults(std::vector<CellInterface*>* /*cellInterfacesLvl*/) { try { throw ErrorECOGEN("writeResults not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };

protected:

  void writeScriptGnuplot(const int& dim);
  void writeScriptGnuplot(const std::string& varName);

  std::string createFilenameGNU(const char* name, int lvl = -1, int proc = -1, int numFichier = -1, std::string nameVariable = "defaut") const;
  void printBlocGnuplot(std::ofstream& fileStream, int& index, const int& dim);

  std::string formatVarNameStyle(std::string const& strToFormat) const;

  std::string m_fileNameVisu;
  std::string m_folderScriptGnuplot;
};

#endif //OUTPUTGNU_H