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

#include "OutputXML.h"
#include "../Run.h"

using namespace tinyxml2;

//***********************************************************************

OutputXML::OutputXML(std::string casTest, std::string run, XMLElement* element, std::string fileName, Input *entree) :
  Output(casTest, run, element, fileName, entree)
{}

//***********************************************************************

OutputXML::~OutputXML(){}

//***********************************************************************

void OutputXML::prepareSortieSpecifique()
{
  try {
    std::ofstream fileStream;
    //Creation du file de sortie collection Paraview
    m_fichierCollectionParaview = m_folderOutput + creationNameFichierXML(m_fileNameCollectionParaview.c_str());
    fileStream.open(m_fichierCollectionParaview.c_str(), std::ios::trunc);
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_fichierCollectionParaview, __FILE__, __LINE__); }
    fileStream << "<?xml version=\"1.0\"?>" << std::endl;
    fileStream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"";
    if (!m_ecritBinaire) { fileStream << "LittleEndian\" "; }
    else { fileStream << m_endianMode.c_str() << "\" "; }
    fileStream << "compressor=\"vtkZLibDataCompressor\">";
    fileStream << std::endl << "  <Collection>" << std::endl;
    fileStream << std::endl << "  </Collection>" << std::endl;
    fileStream << "</VTKFile>" << std::endl;
    fileStream.close();
    //Creation du file de sortie collection VisIt
    m_fichierCollectionVisIt = m_folderOutput + creationNameFichierXML(m_fileNameCollectionVisIt.c_str(), 0, -1, -1, "visit");
    fileStream.open(m_fichierCollectionVisIt.c_str(), std::ios::trunc);
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_fichierCollectionVisIt, __FILE__, __LINE__); }
    fileStream << "!NBLOCKS " << Ncpu << std::endl;
    fileStream.close();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void OutputXML::ecritSolution(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
  try {
    //Ecriture des fichiers de sortie au format XML
    ecritSolutionXML(mesh, cellsLvl);
    //Ajout du file Collection pour grouper les niveaux, les temps, les CPU, etc.
    if (rankCpu == 0) { ecritCollectionXML(mesh); }
  }
  catch (ErrorECOGEN &) { throw; } // Renvoi au niveau suivant
  m_numFichier++;
}

//***********************************************************************

void OutputXML::readResults(Mesh *mesh, std::vector<Cell*>* cellsLvl)
{
  try {
    //1) Parsing XML file
    //-------------------
    std::stringstream fileName(m_folderDatasets + creationNameFichierXML(m_fileNameResults.c_str(), mesh, rankCpu, m_numFichier));
    XMLDocument xmlMain;
    XMLError error(xmlMain.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);
    
    //2) Entering XML file according to mesh
    //--------------------------------------
    XMLElement* nodeVTK, *nodeGrid, *nodePiece, *nodeCellData;
    nodeVTK = xmlMain.FirstChildElement("VTKFile");
    if (nodeVTK == NULL) throw ErrorXMLRacine("VTKFile", fileName.str(), __FILE__, __LINE__);
    //Depend on mesh
    switch (mesh->getType()) {
    case REC:
      nodeGrid = nodeVTK->FirstChildElement("RectilinearGrid");
      if (nodeGrid == NULL) throw ErrorXMLRacine("RectilinearGrid", fileName.str(), __FILE__, __LINE__);
      break;
    case UNS:
      nodeGrid = nodeVTK->FirstChildElement("UnstructuredGrid");
      if (nodeGrid == NULL) throw ErrorXMLRacine("UnstructuredGrid", fileName.str(), __FILE__, __LINE__);
      break;
    case AMR:
      nodeGrid = nodeVTK->FirstChildElement("UnstructuredGrid");
      if (nodeGrid == NULL) throw ErrorXMLRacine("UnstructuredGrid", fileName.str(), __FILE__, __LINE__);
      break;
    default:
      throw ErrorECOGEN("Output::readResults: unknown mesh type", __FILE__, __LINE__); break;
    }
    nodePiece = nodeGrid->FirstChildElement("Piece");
    if (nodePiece == NULL) throw ErrorXMLRacine("Piece", fileName.str(), __FILE__, __LINE__);
    nodeCellData = nodePiece->FirstChildElement("CellData");
    if (nodeCellData == NULL) throw ErrorXMLRacine("CellData", fileName.str(), __FILE__, __LINE__);

    //3) Reading fluid data
    //---------------------
    ReadDonneesPhysiquesXML(mesh, cellsLvl, nodeCellData, fileName.str());
    m_numFichier++;
  } //Fin try
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void OutputXML::ReadDonneesPhysiquesXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, XMLElement* nodeCellData, std::string fileName)
{
  int totalCellsLvlSize(0);
  for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
    totalCellsLvlSize += cellsLvl[lvl].size();
  }
  std::vector<double> scalarDataSet(totalCellsLvlSize);
  std::vector<double> vectorDataSet(3 * totalCellsLvlSize);
  XMLElement* nodeData;

  try {

    nodeData = nodeCellData->FirstChildElement("DataArray");
    if (nodeData == NULL) throw ErrorXMLRacine("DataArray", fileName, __FILE__, __LINE__);

    //1) Reading phases data
    //----------------------
    for (int phase = 0; phase < m_run->getNumberPhases(); phase++)
    {
      //Reading scalars
      for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberScalars(); var++) {
        std::istringstream data(nodeData->GetText());
        this->getJeuDonnees(data, scalarDataSet);
        mesh->setDataSet(scalarDataSet, cellsLvl, var, phase);
        nodeData = nodeData->NextSiblingElement("DataArray");
      }
      //Reading vectors
      for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++) {
        std::istringstream data(nodeData->GetText());
        this->getJeuDonnees(data, vectorDataSet);
        mesh->setDataSet(vectorDataSet, cellsLvl, -var, phase);
        nodeData = nodeData->NextSiblingElement("DataArray");
      }
    } //Fin phase

    //2) Reading mixture data
    //-----------------------
    if (m_run->m_numberPhases > 1) {
      int mixture = -1;
      //Reading scalars
      for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++) {
        std::istringstream data(nodeData->GetText());
        this->getJeuDonnees(data, scalarDataSet);
        mesh->setDataSet(scalarDataSet, cellsLvl, var, mixture);
        nodeData = nodeData->NextSiblingElement("DataArray");
      }
      //Reading vectors
      for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++) {
        std::istringstream data(nodeData->GetText());
        this->getJeuDonnees(data, vectorDataSet);
        mesh->setDataSet(vectorDataSet, cellsLvl, -var, mixture);
        nodeData = nodeData->NextSiblingElement("DataArray");
      }
    } //Fin mixture
    
    //3) Transported data
    //-------------------
    int transport = -2;
    for (int var = 1; var <= m_run->m_numberTransports; var++) {
      std::istringstream data(nodeData->GetText());
      this->getJeuDonnees(data, scalarDataSet);
      mesh->setDataSet(scalarDataSet, cellsLvl, var, transport);
      nodeData = nodeData->NextSiblingElement("DataArray");
    }

    //4) xi indicator for AMR
    //-----------------------
    if (mesh->getType() == AMR) {
      int xi = -3;
      std::istringstream data(nodeData->GetText());
      this->getJeuDonnees(data, scalarDataSet);
      mesh->setDataSet(scalarDataSet, cellsLvl, 1, xi);
      nodeData = nodeData->NextSiblingElement("DataArray");
    }
  } //Fin try
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

std::string OutputXML::creationNameFichierXML(const char* name, Mesh *mesh, int proc, int numFichier, std::string nameVariable)
{
  std::stringstream num;
  std::string prefix;
  if (proc==-2) { prefix = "p"; }
  else { prefix = ""; }

  try {
    //FP//TODO Donnees separees a faire...
    if (m_donneesSeparees) { throw ErrorECOGEN("OutputXML::creationNameFichierXML : donnees Separees non prevu", __FILE__, __LINE__); }
    num << name;
    //Gestion nameVariable
    if ((nameVariable != "defaut") && (nameVariable != "visit")) num << "_" << nameVariable << "_";
    //Gestion binary
    if (m_ecritBinaire) num << "B64";
    //Gestion cpu
    if (proc > -1) num << "_CPU" << proc;
    //Gestion number de file resultat
    if (numFichier != -1) num << "_TIME" << numFichier;
    //Gestion extension
    if(mesh == 0) {
      if (nameVariable == "defaut") num << ".pvd";   //Extension pour la collection Paraview
      if (nameVariable == "visit")  num << ".visit"; //Extension pour la collection VisIt
    }
    else {
      num << "." << prefix;
      switch (mesh->getType()) {
      case REC:
        num << "vtr"; break;
      case UNS:
        num << "vtu"; break;
      case AMR:
        num << "vtu"; break;
        //num << "vtp"; break;
      default:
        throw ErrorECOGEN("OutputXML::creationNameFichierXML : type mesh inconnu", __FILE__, __LINE__);
      }
    }
  }
  catch (ErrorECOGEN &) { throw; }

  return num.str();
}

//***********************************************************************

void OutputXML::ecritSolutionXML(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
  std::ofstream fileStream;

  try {
      
      //1) Ouverture / creation file
      //-------------------------------
      std::string file = m_folderDatasets + creationNameFichierXML(m_fileNameResults.c_str(), mesh, rankCpu, m_numFichier);
      fileStream.open(file.c_str(), std::ios::trunc);
      if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + file, __FILE__, __LINE__); }
      fileStream << "<?xml version=\"1.0\"?>" << std::endl;
      
      //2) Ecriture du mesh
      //-----------------------
      switch (mesh->getType()) {
      case REC:
        ecritMeshRectilinearXML(mesh, fileStream); break;
      case UNS:
        ecritMeshUnstructuredXML(mesh, cellsLvl, fileStream); break;
      case AMR:
        ecritMeshUnstructuredXML(mesh, cellsLvl, fileStream); break;
      default:
        throw ErrorECOGEN("Output::ecritSolutionXML : type mesh inconnu", __FILE__, __LINE__); break;
      }
      
      //3) Ecriture des donnees phases fluids
      //-------------------------------------
      ecritDonneesPhysiquesXML(mesh, cellsLvl, fileStream);
      
      //4) Finalisation file
      //-----------------------
      switch (mesh->getType()) {
      case REC:
        ecritFinFichierRectilinearXML(fileStream); break;
      case UNS:
        ecritFinFichierUnstructuredXML(fileStream); break;
      case AMR:
        ecritFinFichierUnstructuredXML(fileStream); break;
        //ecritFinFichierPolyDataXML(fileStream); break;
      default:
        throw ErrorECOGEN("Output::ecritSolutionXML : type mesh inconnu", __FILE__, __LINE__); break;
      }
      fileStream.close();

  } //Fin try
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void OutputXML::ecritCollectionXML(Mesh *mesh)
{
  try {
    std::ofstream fileStream;
    //ifstream fileStream2((m_folderOutput + m_infoCalcul).c_str()); //For real-time file name
    //double realTime, a, b, c, d, e, f, g, h, i, j, k, l, m;          //For real-time file name
    //double realTime, a, b, c, d;                                     //For real-time file name
    //Creation du file de sortie collection Paraview
    fileStream.open(m_fichierCollectionParaview.c_str(), std::ios::trunc);
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_fichierCollectionParaview, __FILE__, __LINE__); }
    fileStream << "<?xml version=\"1.0\"?>" << std::endl;
    fileStream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"";
    if (!m_ecritBinaire) { fileStream << "LittleEndian\" "; }
    else { fileStream << m_endianMode.c_str() << "\" "; }
    fileStream << "compressor=\"vtkZLibDataCompressor\">";
    fileStream << std::endl << "    <Collection>" << std::endl;
    for (int time = 0; time <= m_numFichier; time++) {
      //fileStream2 >> a >> b >> realTime >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l >> m; //For real-time file name
      //fileStream2 >> a >> b >> realTime >> c >> d;                                              //For real-time file name
      for (int p = 0; p < Ncpu; p++) {
          std::string file = "datasets/" + creationNameFichierXML(m_fileNameResults.c_str(), mesh, p, time);
          fileStream << "        <DataSet timestep=\"" << time << "\" part=\"" << p << "\" file=\"" << file.c_str() << "\"/>" << std::endl;
          //fileStream << "        <DataSet timestep=\"" << realTime << "\" part=\"" << p << "\" file=\"" << file.c_str() << "\"/>" << std::endl; //For real-time file name
      }
    }
    fileStream << "    </Collection>" << std::endl;
    fileStream << "</VTKFile>" << std::endl;
    fileStream.close();
    //fileStream2.close(); //For real-time file name

    //Creation du file de sortie collection VisIt
    fileStream.open(m_fichierCollectionVisIt.c_str(), std::ios::trunc);
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_fichierCollectionVisIt, __FILE__, __LINE__); }
    fileStream << "!NBLOCKS " << Ncpu << std::endl;
    for (int time = 0; time <= m_numFichier; time++) {
      for (int p = 0; p < Ncpu; p++) {
          std::string file = "datasets/" + creationNameFichierXML(m_fileNameResults.c_str(), mesh, p, time);
          fileStream << file.c_str() << std::endl;
      }
    }
    fileStream.close();
  }
  catch (ErrorXML &) { throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void OutputXML::ecritDonneesPhysiquesXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, std::ofstream &fileStream, bool parallel)
{
  std::vector<double> jeuDonnees;

  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  std::string format = "ascii";
  if (m_ecritBinaire) format = "binary";

  fileStream << "      <" << prefix << "CellData>" << std::endl;

  //1) Ecriture des variables des phases
  //------------------------------------
  for (int phase = 0; phase < m_run->getNumberPhases(); phase++) //For complete output
  //for (int phase = 0; phase < 1; phase++) //For reduced output
  {
    std::string eosName(m_cellRef.getPhase(phase)->getEos()->getName());
    eosName.erase(eosName.end()-4, eosName.end());
    //Ecriture des variables scalars
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberScalars(); var++) {
      fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"F" << phase << "_" << m_cellRef.getPhase(phase)->returnNameScalar(var) << "_" << eosName << "\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->recupereDonnees(cellsLvl, jeuDonnees, var, phase);
        this->ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
    }
    //Ecriture des variables vectorielles
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++)
    {
      fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"F" << phase << "_" << m_cellRef.getPhase(phase)->returnNameVector(var) << "_" << eosName << "\" NumberOfComponents=\"3\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->recupereDonnees(cellsLvl, jeuDonnees, -var, phase);
        this->ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
    }
  } //Fin phase

  //2) Ecriture des donnees mixture
  //-------------------------------
  if (m_run->m_numberPhases > 1)
  {
    int mixture = -1;
    //Ecriture des variables scalars du mixture
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++)
    {
      fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"" << m_cellRef.getMixture()->returnNameScalar(var) << "\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->recupereDonnees(cellsLvl, jeuDonnees, var, mixture);
        this->ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
    }
    //Ecriture des variables vectorielles du mixture
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++)
    {
      fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"" << m_cellRef.getMixture()->returnNameVector(var) << "\" NumberOfComponents=\"3\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->recupereDonnees(cellsLvl, jeuDonnees, -var, mixture);
        this->ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
    }
  } //Fin mixture

  //3) Ecriture des transports et autres...
  //---------------------------------------
  int transport = -2; //For complete output
  for (int var = 1; var <= m_run->m_numberTransports; var++)
  {
    fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"T" << var << "\"";
    if (!parallel) {
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->recupereDonnees(cellsLvl, jeuDonnees, var, transport);
      this->ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </" << prefix << "DataArray>" << std::endl;
    }
    else { fileStream << "\"/>" << std::endl; }
  }

  //4) Ecriture indicateur xi
  //-------------------------
  if (mesh->getType() == AMR) { //For complete output
    int xi = -3;
    fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"Xi\"";
    if (!parallel) {
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->recupereDonnees(cellsLvl, jeuDonnees, 1, xi);
      this->ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </" << prefix << "DataArray>" << std::endl;
    }
    else { fileStream << "\"/>" << std::endl; }
  }

  //5) Ecriture gradient rho
  //------------------------
  int gradRho = -4;
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"gradRho\"";
  if (!parallel) {
    fileStream << " format=\"" << format << "\">" << std::endl;
    mesh->recupereDonnees(cellsLvl, jeuDonnees, 1, gradRho);
    this->ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "\"/>" << std::endl; }

  //6) Absolute velocity printing for Moving Reference Frame computations
  //---------------------------------------------------------------------
  if (m_run->m_MRF!=-1) {
    fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"absoluteVelocityMRF\" NumberOfComponents=\"3\"";
    if (!parallel) {
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->extractAbsVelocityMRF(cellsLvl, jeuDonnees, m_run->m_sources[m_run->m_MRF]);
      this->ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </" << prefix << "DataArray>" << std::endl;
    }
    else { fileStream << "\"/>" << std::endl; }
  }

  //7) CPU rank
  //-----------
  // if (mesh->getType() == AMR) { //For complete output
  //   int CPUrank = -5;
  //   fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"CPUrank\"";
  //   if (!parallel) {
  //     fileStream << " format=\"" << format << "\">" << std::endl;
  //     mesh->recupereDonnees(cellsLvl, jeuDonnees, 1, CPUrank);
  //     this->ecritJeuDonnees(jeuDonnees, fileStream, INT);
  //     fileStream << std::endl;
  //     fileStream << "        </" << prefix << "DataArray>" << std::endl;
  //   }
  //   else { fileStream << "\"/>" << std::endl; }
  // }

  //8) Morton index
  //---------------
  // if (mesh->getType() == AMR) { //For complete output
  //   int MortonIndex = -6;
  //   fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"MortonIndex\"";
  //   if (!parallel) {
  //     fileStream << " format=\"" << format << "\">" << std::endl;
  //     mesh->recupereDonnees(cellsLvl, jeuDonnees, 1, MortonIndex);
  //     this->ecritJeuDonnees(jeuDonnees, fileStream, INT);
  //     fileStream << std::endl;
  //     fileStream << "        </" << prefix << "DataArray>" << std::endl;
  //   }
  //   else { fileStream << "\"/>" << std::endl; }
  // }

  //Fin
  fileStream << "      </" << prefix << "CellData>" << std::endl;
}

//***********************************************************************

void OutputXML::ecritMeshRectilinearXML(Mesh* mesh, std::ofstream& fileStream, bool parallel)
{
  std::vector<double> jeuDonnees;

  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }
  
  //0) Header
  //---------
  fileStream << "<VTKFile type=\"" << prefix << "RectilinearGrid\" version=\"0.1\" byte_order=\"";
  if (!m_ecritBinaire) fileStream << "LittleEndian\">" << std::endl;
  else fileStream << m_endianMode.c_str() << "\">" << std::endl;
  if (!parallel) {
    fileStream << "  <RectilinearGrid WholeExtent=\"" << mesh->recupereChaineExtent() << "\">" << std::endl;
    fileStream << "    <Piece Extent=\"" << mesh->recupereChaineExtent() << "\">" << std::endl;
  }
  else {
    fileStream << "  <PRectilinearGrid WholeExtent = \"" << mesh->recupereChaineExtent(true) << "\" GhostLevel=\"0\">" << std::endl;
  }
  
  //1) Ecriture des Coordonnees des noeuds
  //--------------------------------------
  fileStream << "      <" << prefix << "Coordinates>" << std::endl;
  //Coordonnees en X
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    jeuDonnees.clear();
    mesh->recupereCoord(jeuDonnees, X);
    ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  //Coordonnees en Y
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    jeuDonnees.clear();
    mesh->recupereCoord(jeuDonnees, Y);
    ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  //Coordonnees en Z
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    jeuDonnees.clear();
    mesh->recupereCoord(jeuDonnees, Z);
    ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  fileStream << "      </" << prefix << "Coordinates>" << std::endl;
}

//***********************************************************************

void OutputXML::ecritMeshUnstructuredXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, std::ofstream &fileStream, bool parallel)
{
  std::vector<double> jeuDonnees;

  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  //0) Header
  //---------
  fileStream << "<VTKFile type=\"" << prefix << "UnstructuredGrid\" version=\"0.1\" byte_order=\"";
  if (!m_ecritBinaire) fileStream << "LittleEndian\">" << std::endl;
  else fileStream << m_endianMode.c_str() << "\">" << std::endl;

  if (parallel) {
    fileStream << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
  }
  else {
    fileStream << "  <UnstructuredGrid>" << std::endl;
    mesh->ecritHeaderPiece(fileStream, cellsLvl);
  }
  
  //1) Ecriture des Noeuds
  //----------------------
  fileStream << "      <" << prefix << "Points>" << std::endl;
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" NumberOfComponents=\"3\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    jeuDonnees.clear();
    mesh->recupereNoeuds(jeuDonnees, cellsLvl);
    ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  fileStream << "      </" << prefix << "Points>" << std::endl;

  //2) Ecriture des Cells
  //------------------------
  fileStream << "      <" << prefix << "Cells>" << std::endl;
  //Connectivite
  fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"connectivity\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    jeuDonnees.clear();
    mesh->recupereConnectivite(jeuDonnees, cellsLvl);
    ecritJeuDonnees(jeuDonnees, fileStream, INT);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  //Offsets
  fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"offsets\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    jeuDonnees.clear();
    mesh->recupereOffsets(jeuDonnees, cellsLvl);
    ecritJeuDonnees(jeuDonnees, fileStream, INT);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  //Type de cells
  fileStream << "        <" << prefix << "DataArray type=\"UInt8\" Name=\"types\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    jeuDonnees.clear();
    mesh->recupereTypeCell(jeuDonnees, cellsLvl);
    ecritJeuDonnees(jeuDonnees, fileStream, CHAR);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  fileStream << "      </" << prefix << "Cells>" << std::endl;
}

//***********************************************************************

void OutputXML::ecritFinFichierRectilinearXML(std::ofstream &fileStream, bool parallel)
{
  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  if (!parallel) fileStream << "    </Piece>" << std::endl;
  fileStream << "  </" << prefix << "RectilinearGrid>" << std::endl;
  fileStream << "</VTKFile>" << std::endl;
}

//***********************************************************************

void OutputXML::ecritFinFichierUnstructuredXML(std::ofstream &fileStream, bool parallel)
{
  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  if (!parallel) fileStream << "    </Piece>" << std::endl;
  fileStream << "  </" << prefix << "UnstructuredGrid>" << std::endl;
  fileStream << "</VTKFile>" << std::endl;
}

//***********************************************************************

//Old

//***********************************************************************

// void OutputXML::ecritMeshPolyDataXML(Mesh *mesh, vector<Cell*>* cellsLvl, std::ofstream &fileStream, const int& lvl, bool parallel)
// {
//   vector<double> jeuDonnees;

//   string prefix;
//   if (parallel) { prefix = "P"; }
//   else { prefix = ""; }

//   //0) Header
//   //---------
//   fileStream << "<VTKFile type=\"" << prefix << "PolyData\" version=\"0.1\" byte_order=\"";
//   if (!m_ecritBinaire) fileStream << "LittleEndian\">" << endl;
//   else fileStream << m_endianMode.c_str() << "\">" << endl;

//   if (parallel) {
//     fileStream << "  <PPolyData GhostLevel=\"0\">" << endl;
//   }
//   else {
//     fileStream << "  <PolyData>" << endl;
//     mesh->ecritHeaderPiece(fileStream, cellsLvl, lvl);
//   }

//   //1) Ecriture des Noeuds
//   //----------------------
//   fileStream << "      <" << prefix << "Points>" << endl;
//   fileStream << "        <" << prefix << "DataArray type=\"Float64\" NumberOfComponents=\"3\" ";
//   if (!parallel) {
//     if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << endl << "          "; }
//     else { fileStream << "format=\"binary\">" << endl; }
//     jeuDonnees.clear();
//     mesh->recupereNoeuds(jeuDonnees, cellsLvl);
//     ecritJeuDonnees(jeuDonnees, fileStream, DOUBLE);
//     fileStream << endl;
//     fileStream << "        </" << prefix << "DataArray>" << endl;
//   }
//   else { fileStream << "/>" << endl; }
//   fileStream << "      </" << prefix << "Points>" << endl;

//   //2) Ecriture des Polys
//   //---------------------
//   fileStream << "      <" << prefix << "Polys>" << endl;
//   //Connectivite
//   fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"connectivity\" ";
//   if (!parallel) {
//     if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << endl << "          "; }
//     else { fileStream << "format=\"binary\">" << endl; }
//     jeuDonnees.clear();
//     mesh->recupereConnectivite(jeuDonnees, lvl);
//     ecritJeuDonnees(jeuDonnees, fileStream, INT);
//     fileStream << endl;
//     fileStream << "        </" << prefix << "DataArray>" << endl;
//   }
//   else { fileStream << "/>" << endl; }
//   //Offsets
//   fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"offsets\" ";
//   if (!parallel) {
//     if (!m_ecritBinaire) { fileStream << "format=\"ascii\">" << endl << "          "; }
//     else { fileStream << "format=\"binary\">" << endl; }
//     jeuDonnees.clear();
//     mesh->recupereOffsets(jeuDonnees, lvl);
//     ecritJeuDonnees(jeuDonnees, fileStream, INT);
//     fileStream << endl;
//     fileStream << "        </" << prefix << "DataArray>" << endl;
//   }
//   else { fileStream << "/>" << endl; }
//   fileStream << "      </" << prefix << "Polys>" << endl;
// }

// //***********************************************************************

// void OutputXML::ecritFinFichierPolyDataXML(std::ofstream &fileStream, bool parallel)
// {
//   string prefix;
//   if (parallel) { prefix = "P"; }
//   else { prefix = ""; }

//   if (!parallel) fileStream << "    </Piece>" << endl;
//   fileStream << "  </" << prefix << "PolyData>" << endl;
//   fileStream << "</VTKFile>" << endl;
// }

// //***********************************************************************

// void OutputXML::ecritFichierParallelXML(Mesh *mesh, vector<Cell*>* cellsLvl)
// {
//   ofstream fileStream;
//   stringstream num;

//   bool parallel(true);

//   try {
//     //Preparation AMR
//     for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {

//       //1) Ouverture / creation file
//       //-------------------------------
//       if (!m_donneesSeparees) {
//         if (!m_ecritBinaire) { num << m_folderDatasets << "/result_" << m_numFichier; }
//         else { num << m_folderDatasets << "/resultB64_" << m_numFichier; }
//         num << "_lvl" << lvl;
//       }
//       else { throw ErrorECOGEN("Output::ecritFichierParallelXML : donnees Separees non prevu", __FILE__, __LINE__); }
//       switch (mesh->getType()) {
//       case REC:
//         num << ".pvtr"; break;
//       case UNS:
//         num << ".pvtu"; break;
//       case AMR:
//         num << ".pvtu"; break;
//         //num << ".pvtp"; break;
//       default:
//         throw ErrorECOGEN("Output::ecritFichierParallelXML : type mesh inconnu", __FILE__, __LINE__); break;
//       }

//       fileStream.open((num.str()).c_str(), ios::trunc);
//       if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + num.str(), __FILE__, __LINE__); }
//       fileStream << "<?xml version=\"1.0\"?>" << endl;

//       //2) Ecriture des infos mesh
//       //------------------------------
//       switch (mesh->getType()) {
//       case REC:
//         ecritMeshRectilinearXML(mesh, fileStream, parallel); break;
//       case UNS:
//         ecritMeshUnstructuredXML(mesh, cellsLvl, fileStream, lvl, parallel); break;
//       case AMR:
//         ecritMeshUnstructuredXML(mesh, cellsLvl, fileStream, lvl, parallel); break;
//         //ecritMeshPolyDataXML(mesh, cellsLvl, fileStream, lvl, parallel); break;
//       default:
//         throw ErrorECOGEN("Output::ecritSolutionXML : type mesh inconnu", __FILE__, __LINE__); break;
//       }

//       //3) Ecriture des donnees phases fluids
//       //-------------------------------------
//       ecritDonneesPhysiquesXML(mesh, cellsLvl, fileStream, lvl, parallel);

//       //4) Ecriture des name des fichiers
//       //--------------------------------
//       for (int p = 0; p < Ncpu; p++)
//       {
//         stringstream fichTemp;
//         if (m_ecritBinaire) { fichTemp << "resultB64_"; }
//         else { fichTemp << "result_"; }
//         //fichTemp << m_numFichier << "_lvl" << lvl << "_CPU" << rank;
//         fichTemp << "CPU" << rankCpu << "_lvl" << lvl << "_" << m_numFichier;

//         switch (mesh->getType()) {
//         case REC:
//           fichTemp << ".vtr";
//           fileStream << "    <Piece Extent=\"" << mesh->recupereChaineExtent() << "\" Source=\"" << fichTemp.str() << "\"/>" << endl;
//           break;
//         case UNS:
//           fichTemp << ".vtu";
//           fileStream << "    <Piece Source=\"" << fichTemp.str() << "\"/>" << endl;
//           break;
//         case AMR:
//           fichTemp << ".vtu";
//           //fichTemp << ".vtp";
//           fileStream << "    <Piece Source=\"" << fichTemp.str() << "\"/>" << endl;
//           break;
//         default:
//           throw ErrorECOGEN("Output::ecritSolutionXML : type mesh inconnu", __FILE__, __LINE__); break;
//         }
//       }

//       //5) Finalisation file
//       //-----------------------
//       switch (mesh->getType()) {
//       case REC:
//         ecritFinFichierRectilinearXML(fileStream, parallel); break;
//       case UNS:
//         ecritFinFichierUnstructuredXML(fileStream, parallel); break;
//       case AMR:
//         ecritFinFichierUnstructuredXML(fileStream, parallel); break;
//         //ecritFinFichierPolyDataXML(fileStream, parallel); break;
//       default:
//         throw ErrorECOGEN("Output::ecritSolutionXML : type mesh inconnu", __FILE__, __LINE__); break;
//       }
//       fileStream.close();

//     }

//   } //Fin try
//   catch (ErrorECOGEN &) { throw; }
// }

//***********************************************************************
