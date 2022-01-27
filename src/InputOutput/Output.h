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

#ifndef OUTPUT_H
#define OUTPUT_H

//Macro for system interactions (creation/destruction of folders)
#ifdef WIN32
  #include <direct.h>
#else
  #include <sys/types.h>
  #include <sys/stat.h>
#endif

#include <fstream>
#include "../libTierces/tinyxml2.h"
#include "../Errors.h"
#include "../Meshes/HeaderMesh.h"
#include "../Order1/Cell.h"
#include "IO.h"

class Output;

#include "Input.h"

class Output
{
  public:
    //! \brief   Default constructor for specific output without specific needs
    Output();

    //! \brief   Main constructor for datasets used for OutputXML and OutputGNU according to outputMode
    //! \param   casTest   Test case name (defined in "main.xml")  
    //! \param   nameRun   Folder to store results
    //! \param   element   XML outputMode element
    //! \param   fileName  Full path to mainVX.xml of current test case
    //! \param   entree    Input pointer to access run pointer and its information
    Output(std::string casTest, std::string nameRun, tinyxml2::XMLElement* element, std::string fileName, Input* entree);

    //! \brief   Constructor for specific derived GNU outputs (boundary, probe, cut)
    //! \param   element   XML GNU output element to get stream precision
    Output(tinyxml2::XMLElement* element);

    virtual ~Output();

    virtual void locateProbeInMesh(const TypeMeshContainer<Cell*>& /*cells*/, const int& /*nbCells*/, bool /*localSeeking*/ = false) { try { throw ErrorECOGEN("locateProbeInMesh not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };
    virtual Cell* locateProbeInAMRSubMesh(std::vector<Cell*>* /*cells*/, const int& /*nbCells*/) { try { throw ErrorECOGEN("locateProbeInMesh not available for requested output format"); } catch (ErrorECOGEN&) { throw; } return 0; };

    void initializeOutput(const Cell& cell);
    void initializeOutput(std::vector<CellInterface*>* cellInterfacesLvl); //!< Currently only used for OutputBoundaryMassflowGNU
    virtual void initializeOutputInfos();
    virtual void writeResults(Mesh* /*mesh*/, std::vector<Cell*>* /*cellsLvl*/) { try { throw ErrorECOGEN("writeResults not available for requested output format"); } catch (ErrorECOGEN&) { throw; }};
    virtual void writeResults(std::vector<CellInterface*>* /*cellInterfacesLvl*/) { try { throw ErrorECOGEN("writeResults not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };
    void printTree(Mesh* mesh, std::vector<Cell*>* cellsLvl, int m_restartAMRsaveFreq);
    virtual void writeInfos();
    void saveInfosMailles() const;

    virtual void initializeSpecificOutput() { try { throw ErrorECOGEN("initializeSpecificOutput not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };
    virtual void initializeSpecificOutput(std::vector<CellInterface*>* /*cellInterfacesLvl*/) { try { throw ErrorECOGEN("initializeSpecificOutput not available for requested output format"); } catch (ErrorECOGEN&) { throw; } }; //!< Currently only used for OutputBoundaryMassflowGNU

    void readInfos();
    virtual void readResults(Mesh* /*mesh*/, std::vector<Cell*>* /*cellsLvl*/) { try { throw ErrorECOGEN("readResutls not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };
    void readDomainDecompostion(Mesh* mesh);
    void readTree(Mesh *mesh, TypeMeshContainer<Cell*>* cellsLvl, TypeMeshContainer<Cell*>* cellsLvlGhost, TypeMeshContainer<CellInterface*>* cellInterfacesLvl,
        const std::vector<AddPhys*>& addPhys, int& nbCellsTotalAMR);

    //Accessor
    int getNumSortie() const { return m_numFichier; };
    virtual double getNextTime() { try { throw ErrorECOGEN("getNextTime not available for requested output format"); } catch (ErrorECOGEN&) { throw; } return 0.; }
    virtual bool possesses() { try { throw ErrorECOGEN("possesses not available for requested output format"); } catch (ErrorECOGEN&) { throw; } return false; };
    const std::string& getFolderOutput(){ return m_folderOutput;}
    const TypeOutput& getType() const { return m_type; };
    bool getReducedOutput() const { return m_reducedOutput; };

  protected:

    //General data
    void afficheInfoEcriture() const;
    void saveInfos() const;
    std::string createFilename(const char* name, int lvl = -1, int proc = -1, int numFichier = -1) const;

    void writeDataset(std::vector<double> jeuDonnees, std::ofstream& fileStream, TypeData typeData);
    void getDataset(std::istringstream& data, std::vector<double>& jeuDonnees);

    Input* m_input;                                     //!<Pointer to input
    Run* m_run;                                         //!<Pointer to run
    TypeOutput m_type;                                  //!<Type of output

    //Attributes names file/folder
    std::string m_simulationName;                       //!<Test case name (defined in "main.xml")
    std::string m_infoCalcul;                           //!<Filename to save useful info of computation
    std::string m_infoMesh;                             //!<Filename of mesh info file
    std::string m_treeStructure;                        //!<Filename for tree structure backup
    std::string m_domainDecomposition;                  //!<Filename for domain decomposition backup
    std::string m_fileNameResults;                      //!<Filename of result file
    std::string m_filenameCollectionParaview;           //!<Name of the collection containing the results files (for Paraview)
    std::string m_filenameCollectionVisIt;              //!<Name of the collection containing the results files (for VisIt)
    std::string m_folderOutput;                         //!<Folder to store results
    std::string m_folderSavesInput;                     //!<Folder to store a copy of input files
    std::string m_folderDatasets;                       //!<Folder to save the datasets
    std::string m_folderInfoMesh;                       //!<Folder to store mesh info
    std::string m_folderCuts;                           //!<Cuts results folder location
    std::string m_folderProbes;                         //!<Probes results folder location
    std::string m_folderGlobalQuantities;               //!<Global quantity (e.g. mass) results folder location
    std::string m_folderBoundaries;                     //!<Boundaries flux results folder location
    std::string m_fileCollectionParaview;               //!<Chemin du file collection regroupant les fichiers resultats (for Paraview)
    std::string m_fileCollectionVisIt;                  //!<Chemin du file collection regroupant les fichiers resultats (for VisIt)
    std::string m_folderErrorsAndWarnings;              //!<File path for errors and warnings
     
    //Attributes of print parameters
    bool m_writeBinary;                                 //!<Choice to write binary/ASCII
    bool m_splitData;                                   //!<Choix print donnees dans des fichiers separes
    int m_precision;                                    //!<Output files precision (number of digits) //default: 0
    bool m_reducedOutput;                               //!<Choice of reduced number of output variables when possible (depends on the model)

    int m_numFichier; 
    std::string m_endianMode;
    
    //Useful to print cell data
    Cell m_cellRef;                                     //!<Reference cell to extract variables name
};

#endif //OUTPUT_H
