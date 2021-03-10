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

//Macro pour les interactions systeme (creation/destruction repertoires)
#ifdef WIN32
  #include <direct.h>
#else
  #include <sys/types.h>
  #include <sys/stat.h>
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
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
    Output();
    Output(std::string casTest, std::string nameRun, tinyxml2::XMLElement* element, std::string fileName, Input* entree);
    virtual ~Output();

    virtual void locateProbeInMesh(const TypeMeshContainer<Cell*>& /*cells*/, const int& /*nbCells*/, bool /*localSeeking*/ = false) { try { throw ErrorECOGEN("locateProbeInMesh not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };
    virtual Cell* locateProbeInAMRSubMesh(std::vector<Cell*>* /*cells*/, const int& /*nbCells*/) { try { throw ErrorECOGEN("locateProbeInMesh not available for requested output format"); } catch (ErrorECOGEN&) { throw; } return 0; };

    void prepareOutput(const Cell& cell);
    void prepareOutput(std::vector<CellInterface*>* cellInterfacesLvl); //!< Currently only used for OutputBoundaryMassflowGNU
    virtual void prepareOutputInfos();
    virtual void ecritSolution(Mesh* /*mesh*/, std::vector<Cell*>* /*cellsLvl*/) { try { throw ErrorECOGEN("ecritSolution not available for requested output format"); } catch (ErrorECOGEN&) { throw; }};
    virtual void ecritSolution(std::vector<CellInterface*>* /*cellInterfacesLvl*/) { try { throw ErrorECOGEN("ecritSolution not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };
    void printTree(Mesh* mesh, std::vector<Cell*>* cellsLvl, int m_restartAMRsaveFreq);
    virtual void ecritInfos();
    void saveInfosMailles() const;

    virtual void prepareSortieSpecifique() { try { throw ErrorECOGEN("prepareSortieSpecifique not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };
    virtual void prepareSortieSpecifique(std::vector<CellInterface*>* /*cellInterfacesLvl*/) { try { throw ErrorECOGEN("prepareSortieSpecifique not available for requested output format"); } catch (ErrorECOGEN&) { throw; } }; //!< Currently only used for OutputBoundaryMassflowGNU

    void readInfos();
    virtual void readResults(Mesh* /*mesh*/, std::vector<Cell*>* /*cellsLvl*/) { try { throw ErrorECOGEN("readResutls not available for requested output format"); } catch (ErrorECOGEN&) { throw; } };
    void readDomainDecompostion(Mesh* mesh);
    void readTree(Mesh *mesh, TypeMeshContainer<Cell*>* cellsLvl, TypeMeshContainer<Cell*>* cellsLvlGhost, TypeMeshContainer<CellInterface*>* cellInterfacesLvl,
        const std::vector<AddPhys*>& addPhys, Model* model, int& nbCellsTotalAMR);

    //Accesseur
    int getNumSortie() const { return m_numFichier; };
    virtual double getNextTime() { try { throw ErrorECOGEN("getNextTime not available for requested output format"); } catch (ErrorECOGEN&) { throw; } return 0.; }
    virtual bool possesses() { try { throw ErrorECOGEN("possesses not available for requested output format"); } catch (ErrorECOGEN&) { throw; } return false; };
    const std::string& getFolderOutput(){ return m_folderOutput;}

  protected:

    //Donnees generales
    void afficheInfoEcriture() const;
    void saveInfos() const;
    std::string creationNameFichier(const char* name, int lvl = -1, int proc = -1, int numFichier = -1) const;

    void ecritJeuDonnees(std::vector<double> jeuDonnees, std::ofstream& fileStream, TypeData typeData);
    void getJeuDonnees(std::istringstream& data, std::vector<double>& jeuDonnees);

    Input* m_input;                                     //!<Pointeur vers entree
    Run* m_run;                                         //!<pointeur vers run

    //attribut name fichiers/dossiers
    std::string m_simulationName;                       //!<Name du cas test (defini dans "main.xml")
    std::string m_infoCalcul;                           //!<Name file pour saver les infos utiles du compute
    std::string m_infoMesh;                             //!<Name fichiers pour stocker les infos de mesh
    std::string m_treeStructure;                        //!<File name for tree structure backup
    std::string m_domainDecomposition;                  //!<File name for domain decomposition backup
    std::string m_fileNameResults;                      //!<Name du file de sortie resultat
    std::string m_fileNameCollectionParaview;           //!<Name de la collection regroupant les fichiers resultats (for Paraview)
    std::string m_fileNameCollectionVisIt;              //!<Name de la collection regroupant les fichiers resultats (for VisIt)
    std::string m_folderOutput;                         //!<Dossier pour enregistrement des resultats
    std::string m_folderSavesInput;                     //!<Dossier pour copier les fichiers entrees
    std::string m_folderDatasets;                       //!<Folder to save the datasets
    std::string m_folderInfoMesh;                       //!<Dossier pour stocker les infos de mesh
    std::string m_folderCuts;                           //!<Cuts results folder location
    std::string m_folderProbes;                         //!<Probes results folder location
    std::string m_folderGlobalQuantities;               //!<Global quantity (e.g. mass) results folder location
    std::string m_folderBoundariesFlux;                 //!<Boundaries flux results folder location
    std::string m_fichierCollectionParaview;            //!<Chemin du file collection regroupant les fichiers resultats (for Paraview)
    std::string m_fichierCollectionVisIt;               //!<Chemin du file collection regroupant les fichiers resultats (for VisIt)
    std::string m_folderErrorsAndWarnings;               //!File path for errors and warnings
     
    //attribut parametres d print
    bool m_ecritBinaire;                                //!<Choix print binary/ASCII
    bool m_donneesSeparees;                             //!<Choix print donnees dans des fichiers separes
    int m_precision;                                    //!<Output files precision (number of digits) //default: 0

    int m_numFichier; 
    std::string m_endianMode;
    
    //Utile pour print des donnees de cells
    Cell m_cellRef;                                     //!<cell de reference pour recup�rer les names des variables
};

#endif //OUTPUT_H
