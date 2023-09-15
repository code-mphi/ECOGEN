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
{
  m_type = TypeOutput::XML;
}

//***********************************************************************

OutputXML::OutputXML(std::string run, int fileNumberRestartMeshMapping, Input *input) :
  Output(run, fileNumberRestartMeshMapping, input)
{
  m_type = TypeOutput::XML;
}

//***********************************************************************

OutputXML::~OutputXML(){}

//***********************************************************************

void OutputXML::initializeSpecificOutput()
{
  try {
    std::ofstream fileStream;
    //Creation du file de sortie collection Paraview
    m_fileCollectionParaview = m_folderOutput + createFilenameXML(m_filenameCollectionParaview.c_str());
    fileStream.open(m_fileCollectionParaview.c_str(), std::ios::trunc);
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_fileCollectionParaview, __FILE__, __LINE__); }
    fileStream << "<?xml version=\"1.0\"?>" << std::endl;
    fileStream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"";
    if (!m_writeBinary) { fileStream << "LittleEndian\" "; }
    else { fileStream << m_endianMode.c_str() << "\" "; }
    fileStream << "compressor=\"vtkZLibDataCompressor\">";
    fileStream << std::endl << "  <Collection>" << std::endl;
    fileStream << std::endl << "  </Collection>" << std::endl;
    fileStream << "</VTKFile>" << std::endl;
    fileStream.close();
    //Creation du file de sortie collection VisIt
    m_fileCollectionVisIt = m_folderOutput + createFilenameXML(m_filenameCollectionVisIt.c_str(), 0, -1, -1, "visit");
    fileStream.open(m_fileCollectionVisIt.c_str(), std::ios::trunc);
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_fileCollectionVisIt, __FILE__, __LINE__); }
    fileStream << "!NBLOCKS " << Ncpu << std::endl;
    fileStream.close();
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void OutputXML::writeResults(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
  try {
    //Write output files with XML format
    writeResultsXML(mesh, cellsLvl);
    //Add du file Collection pour grouper les niveaux, les temps, les CPU, etc.
    if (rankCpu == 0) { writeCollectionXML(mesh); }
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
    std::stringstream fileName(m_folderDatasets + createFilenameXML(m_fileNameResults.c_str(), mesh, rankCpu, m_numFichier));
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
    ReadPhysicalDataXML(mesh, cellsLvl, nodeCellData, fileName.str());
    m_numFichier++;
  } //End try
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void OutputXML::readResultsCpu(Mesh *mesh, std::vector<Cell*>* cellsLvl, int cpu)
{
  try {
    //1) Parsing XML file
    //-------------------
    std::stringstream fileName(m_folderDatasets + createFilenameXML(m_fileNameResults.c_str(), mesh, cpu, m_numFichier));
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
    ReadPhysicalDataXML(mesh, cellsLvl, nodeCellData, fileName.str());
  } //End try
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void OutputXML::ReadPhysicalDataXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, XMLElement* nodeCellData, std::string fileName)
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
        this->getDataset(data, scalarDataSet);
        mesh->setDataSet(scalarDataSet, cellsLvl, var, phase);
        nodeData = nodeData->NextSiblingElement("DataArray");
      }
      //Reading vectors
      for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++) {
        std::istringstream data(nodeData->GetText());
        this->getDataset(data, vectorDataSet);
        mesh->setDataSet(vectorDataSet, cellsLvl, -var, phase);
        nodeData = nodeData->NextSiblingElement("DataArray");
      }
    } //End phase

    //2) Reading mixture data
    //-----------------------
    if (m_run->m_numberPhases > 1) {
      int mixture = -1;
      //Reading scalars
      for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++) {
        std::istringstream data(nodeData->GetText());
        this->getDataset(data, scalarDataSet);
        mesh->setDataSet(scalarDataSet, cellsLvl, var, mixture);
        nodeData = nodeData->NextSiblingElement("DataArray");
      }
      //Reading vectors
      for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++) {
        std::istringstream data(nodeData->GetText());
        this->getDataset(data, vectorDataSet);
        mesh->setDataSet(vectorDataSet, cellsLvl, -var, mixture);
        nodeData = nodeData->NextSiblingElement("DataArray");
      }
    } //End mixture
    
    //3) Transported data
    //-------------------
    int transport = -2;
    for (int var = 1; var <= m_run->m_numberTransports; var++) {
      std::istringstream data(nodeData->GetText());
      this->getDataset(data, scalarDataSet);
      mesh->setDataSet(scalarDataSet, cellsLvl, var, transport);
      nodeData = nodeData->NextSiblingElement("DataArray");
    }

    //4) xi indicator for AMR
    //-----------------------
    if (mesh->getType() == AMR) {
      int xi = -3;
      std::istringstream data(nodeData->GetText());
      this->getDataset(data, scalarDataSet);
      mesh->setDataSet(scalarDataSet, cellsLvl, 1, xi);
      nodeData = nodeData->NextSiblingElement("DataArray");
    }
  } //End try
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

std::string OutputXML::createFilenameXML(const char* name, Mesh *mesh, int proc, int numFichier, std::string nameVariable)
{
  std::stringstream num;
  std::string prefix;
  if (proc==-2) { prefix = "p"; }
  else { prefix = ""; }

  try {
    //FP//TODO Donnees separees a faire...
    if (m_splitData) { throw ErrorECOGEN("OutputXML::createFilenameXML : separated data not done", __FILE__, __LINE__); }
    num << name;
    //Gestion nameVariable
    if ((nameVariable != "defaut") && (nameVariable != "visit")) num << "_" << nameVariable << "_";
    //Gestion binary
    if (m_writeBinary) num << "B64";
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
        throw ErrorECOGEN("OutputXML::createFilenameXML : type mesh unknown", __FILE__, __LINE__);
      }
    }
  }
  catch (ErrorECOGEN &) { throw; }

  return num.str();
}

//***********************************************************************

void OutputXML::writeResultsXML(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
  std::ofstream fileStream;

  try {
      
      //1) Opening and creation of file
      //-------------------------------
      std::string file = m_folderDatasets + createFilenameXML(m_fileNameResults.c_str(), mesh, rankCpu, m_numFichier);
      fileStream.open(file.c_str(), std::ios::trunc);
      if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + file, __FILE__, __LINE__); }
      fileStream << "<?xml version=\"1.0\"?>" << std::endl;
      
      //2) Write mesh
      //-------------
      switch (mesh->getType()) {
      case REC:
        writeMeshRectilinearXML(mesh, fileStream); break;
      case UNS:
        writeMeshUnstructuredXML(mesh, cellsLvl, fileStream); break;
      case AMR:
        writeMeshUnstructuredXML(mesh, cellsLvl, fileStream); break;
      default:
        throw ErrorECOGEN("Output::writeResultsXML : type mesh unknown", __FILE__, __LINE__); break;
      }
      
      //3) Write fluid phase data
      //-------------------------
      writePhysicalDataXML(mesh, cellsLvl, fileStream);
      
      //4) Finalyse file
      //----------------
      switch (mesh->getType()) {
      case REC:
        writeFinFichierRectilinearXML(fileStream); break;
      case UNS:
        writeFinFichierUnstructuredXML(fileStream); break;
      case AMR:
        writeFinFichierUnstructuredXML(fileStream); break;
        //writeFinFichierPolyDataXML(fileStream); break;
      default:
        throw ErrorECOGEN("Output::writeResultsXML : type mesh unknown", __FILE__, __LINE__); break;
      }
      fileStream.close();

  } //End try
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void OutputXML::writeCollectionXML(Mesh *mesh)
{
  try {
    std::ofstream fileStream;
    //ifstream fileStream2((m_folderOutput + m_infoCalcul).c_str()); //For real-time file name
    //double realTime, a, b, c, d, e, f, g, h, i, j, k, l, m;          //For real-time file name
    //double realTime, a, b, c, d;                                     //For real-time file name
    //Creation du file de sortie collection Paraview
    fileStream.open(m_fileCollectionParaview.c_str(), std::ios::trunc);
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_fileCollectionParaview, __FILE__, __LINE__); }
    fileStream << "<?xml version=\"1.0\"?>" << std::endl;
    fileStream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"";
    if (!m_writeBinary) { fileStream << "LittleEndian\" "; }
    else { fileStream << m_endianMode.c_str() << "\" "; }
    fileStream << "compressor=\"vtkZLibDataCompressor\">";
    fileStream << std::endl << "    <Collection>" << std::endl;
    for (int time = 0; time <= m_numFichier; time++) {
      //fileStream2 >> a >> b >> realTime >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l >> m; //For real-time file name
      //fileStream2 >> a >> b >> realTime >> c >> d;                                              //For real-time file name
      for (int p = 0; p < Ncpu; p++) {
          std::string file = "datasets/" + createFilenameXML(m_fileNameResults.c_str(), mesh, p, time);
          fileStream << "        <DataSet timestep=\"" << time << "\" part=\"" << p << "\" file=\"" << file.c_str() << "\"/>" << std::endl;
          //fileStream << "        <DataSet timestep=\"" << realTime << "\" part=\"" << p << "\" file=\"" << file.c_str() << "\"/>" << std::endl; //For real-time file name
      }
    }
    fileStream << "    </Collection>" << std::endl;
    fileStream << "</VTKFile>" << std::endl;
    fileStream.close();
    //fileStream2.close(); //For real-time file name

    //Creation du file de sortie collection VisIt
    fileStream.open(m_fileCollectionVisIt.c_str(), std::ios::trunc);
    if (!fileStream) { throw ErrorECOGEN("Impossible d ouvrir le file " + m_fileCollectionVisIt, __FILE__, __LINE__); }
    fileStream << "!NBLOCKS " << Ncpu << std::endl;
    for (int time = 0; time <= m_numFichier; time++) {
      for (int p = 0; p < Ncpu; p++) {
          std::string file = "datasets/" + createFilenameXML(m_fileNameResults.c_str(), mesh, p, time);
          fileStream << file.c_str() << std::endl;
      }
    }
    fileStream.close();
  }
  catch (ErrorXML &) { throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void OutputXML::writePhysicalDataXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, std::ofstream &fileStream, bool parallel)
{
  std::vector<double> dataset;

  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  std::string format = "ascii";
  if (m_writeBinary) format = "binary";

  fileStream << "      <" << prefix << "CellData>" << std::endl;

  //1) Write des variables des phases
  //---------------------------------
  for (int phase = 0; phase < m_run->getNumberPhases(); phase++) //For complete output
  //for (int phase = 0; phase < 1; phase++) //For reduced output
  {
    std::string eosName;
    if (m_cellRef.getPhase(phase)->getEos() != nullptr) {
      eosName = m_cellRef.getPhase(phase)->getEos()->getName();
      eosName.erase(eosName.end()-4, eosName.end());
    }
    //Write des variables scalars
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberScalars(); var++) {
      fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"F" << phase << "_" << m_cellRef.getPhase(phase)->returnNameScalar(var);
      if (m_cellRef.getPhase(phase)->getEos() != nullptr) { fileStream << "_" << eosName; }
      fileStream << "\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->getData(cellsLvl, dataset, var, phase);
        this->writeDataset(dataset, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
    }
    //Write des variables vectorielles
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++)
    {
      fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"F" << phase << "_" << m_cellRef.getPhase(phase)->returnNameVector(var);
      if (m_cellRef.getPhase(phase)->getEos() != nullptr) { fileStream << "_" << eosName; }
      fileStream << "\" NumberOfComponents=\"3\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->getData(cellsLvl, dataset, -var, phase);
        this->writeDataset(dataset, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
    }
  } //End phase

  //2) Write mixture data
  //---------------------
  if (m_run->m_numberPhases > 1)
  {
    int mixture = -1;
    //Write des variables scalars du mixture
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++)
    {
      fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"" << m_cellRef.getMixture()->returnNameScalar(var) << "\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->getData(cellsLvl, dataset, var, mixture);
        this->writeDataset(dataset, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
    }
    //Write des variables vectorielles du mixture
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++)
    {
      fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"" << m_cellRef.getMixture()->returnNameVector(var) << "\" NumberOfComponents=\"3\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->getData(cellsLvl, dataset, -var, mixture);
        this->writeDataset(dataset, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
    }
  } //End mixture

  //3) Write of transports and others...
  //------------------------------------
  int transport = -2; //For complete output
  for (int var = 1; var <= m_run->m_numberTransports; var++)
  {
    fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"T" << var << "\"";
    if (!parallel) {
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->getData(cellsLvl, dataset, var, transport);
      this->writeDataset(dataset, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </" << prefix << "DataArray>" << std::endl;
    }
    else { fileStream << "\"/>" << std::endl; }
  }

  //4) Write indicateur xi
  //----------------------
  if (mesh->getType() == AMR) { //For complete output
    int xi = -3;
    fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"Xi\"";
    if (!parallel) {
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->getData(cellsLvl, dataset, 1, xi);
      this->writeDataset(dataset, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </" << prefix << "DataArray>" << std::endl;
    }
    else { fileStream << "\"/>" << std::endl; }
  }

  //5) Write gradient rho
  //---------------------
  int gradRho = -4;
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"gradRho\"";
  if (!parallel) {
    fileStream << " format=\"" << format << "\">" << std::endl;
    mesh->getData(cellsLvl, dataset, 1, gradRho);
    this->writeDataset(dataset, fileStream, DOUBLE);
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
      mesh->extractAbsVelocityMRF(cellsLvl, dataset, m_run->m_sources[m_run->m_MRF]);
      this->writeDataset(dataset, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </" << prefix << "DataArray>" << std::endl;
    }
    else { fileStream << "\"/>" << std::endl; }
  }

  //7) Cells' reference length
  //--------------------------
  if (m_run->m_extractRefLength) {
    fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"Reference_length\"";
    if (!parallel) {
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->extractReferenceLength(cellsLvl, dataset);
      this->writeDataset(dataset, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </" << prefix << "DataArray>" << std::endl;
    }
    else { fileStream << "\"/>" << std::endl; }
  }

  //8) CPU rank
  //-----------
  // if (mesh->getType() == AMR) { //For complete output
  //   int CPUrank = -5;
  //   fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"CPUrank\"";
  //   if (!parallel) {
  //     fileStream << " format=\"" << format << "\">" << std::endl;
  //     mesh->getData(cellsLvl, dataset, 1, CPUrank);
  //     this->writeDataset(dataset, fileStream, INT);
  //     fileStream << std::endl;
  //     fileStream << "        </" << prefix << "DataArray>" << std::endl;
  //   }
  //   else { fileStream << "\"/>" << std::endl; }
  // }

  //9) Morton index
  //---------------
  // if (mesh->getType() == AMR) { //For complete output
  //   int MortonIndex = -6;
  //   fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"MortonIndex\"";
  //   if (!parallel) {
  //     fileStream << " format=\"" << format << "\">" << std::endl;
  //     mesh->getData(cellsLvl, dataset, 1, MortonIndex);
  //     this->writeDataset(dataset, fileStream, INT);
  //     fileStream << std::endl;
  //     fileStream << "        </" << prefix << "DataArray>" << std::endl;
  //   }
  //   else { fileStream << "\"/>" << std::endl; }
  // }

  //9) Writing saturation pressure (specific recording)
  //---------------------------------------------------
  if (m_run->m_recordPsat) {
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" Name=\"Psat\"";
      if (!parallel) {
        fileStream << " format=\"" << format << "\">" << std::endl;
        mesh->getData(cellsLvl, dataset, 1, -7);
        this->writeDataset(dataset, fileStream, DOUBLE);
        fileStream << std::endl;
        fileStream << "        </" << prefix << "DataArray>" << std::endl;
      }
      else { fileStream << "\"/>" << std::endl; }
  }

  //Fin
  fileStream << "      </" << prefix << "CellData>" << std::endl;
}

//***********************************************************************

void OutputXML::writeMeshRectilinearXML(Mesh* mesh, std::ofstream& fileStream, bool parallel)
{
  std::vector<double> dataset;

  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }
  
  //0) Header
  //---------
  fileStream << "<VTKFile type=\"" << prefix << "RectilinearGrid\" version=\"0.1\" byte_order=\"";
  if (!m_writeBinary) fileStream << "LittleEndian\">" << std::endl;
  else fileStream << m_endianMode.c_str() << "\">" << std::endl;
  if (!parallel) {
    fileStream << "  <RectilinearGrid WholeExtent=\"" << mesh->getStringExtent() << "\">" << std::endl;
    fileStream << "    <Piece Extent=\"" << mesh->getStringExtent() << "\">" << std::endl;
  }
  else {
    fileStream << "  <PRectilinearGrid WholeExtent = \"" << mesh->getStringExtent(true) << "\" GhostLevel=\"0\">" << std::endl;
  }
  
  //1) Write node coordinates
  //-------------------------
  fileStream << "      <" << prefix << "Coordinates>" << std::endl;
  //Coordinate in X
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" ";
  if (!parallel) {
    if (!m_writeBinary) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    dataset.clear();
    mesh->getCoord(dataset, X);
    writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  //Coordinate in Y
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" ";
  if (!parallel) {
    if (!m_writeBinary) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    dataset.clear();
    mesh->getCoord(dataset, Y);
    writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  //Coordinate in Z
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" ";
  if (!parallel) {
    if (!m_writeBinary) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    dataset.clear();
    mesh->getCoord(dataset, Z);
    writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  fileStream << "      </" << prefix << "Coordinates>" << std::endl;
}

//***********************************************************************

void OutputXML::writeMeshUnstructuredXML(Mesh *mesh, std::vector<Cell*>* cellsLvl, std::ofstream &fileStream, bool parallel)
{
  std::vector<double> dataset;

  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  //0) Header
  //---------
  fileStream << "<VTKFile type=\"" << prefix << "UnstructuredGrid\" version=\"0.1\" byte_order=\"";
  if (!m_writeBinary) fileStream << "LittleEndian\">" << std::endl;
  else fileStream << m_endianMode.c_str() << "\">" << std::endl;

  if (parallel) {
    fileStream << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
  }
  else {
    fileStream << "  <UnstructuredGrid>" << std::endl;
    mesh->writeHeaderPiece(fileStream, cellsLvl);
  }
  
  //1) Write of nodes
  //-----------------
  fileStream << "      <" << prefix << "Points>" << std::endl;
  fileStream << "        <" << prefix << "DataArray type=\"Float64\" NumberOfComponents=\"3\" ";
  if (!parallel) {
    if (!m_writeBinary) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    dataset.clear();
    mesh->getNodes(dataset, cellsLvl);
    writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  fileStream << "      </" << prefix << "Points>" << std::endl;

  //2) Write des Cells
  //------------------
  fileStream << "      <" << prefix << "Cells>" << std::endl;
  //Connectivite
  fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"connectivity\" ";
  if (!parallel) {
    if (!m_writeBinary) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    dataset.clear();
    mesh->getConnectivity(dataset, cellsLvl);
    writeDataset(dataset, fileStream, INT);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  //Offsets
  fileStream << "        <" << prefix << "DataArray type=\"Int32\" Name=\"offsets\" ";
  if (!parallel) {
    if (!m_writeBinary) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    dataset.clear();
    mesh->getOffsets(dataset, cellsLvl);
    writeDataset(dataset, fileStream, INT);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  //Type de cells
  fileStream << "        <" << prefix << "DataArray type=\"UInt8\" Name=\"types\" ";
  if (!parallel) {
    if (!m_writeBinary) { fileStream << "format=\"ascii\">" << std::endl << "          "; }
    else { fileStream << "format=\"binary\">" << std::endl; }
    dataset.clear();
    mesh->getTypeCell(dataset, cellsLvl);
    writeDataset(dataset, fileStream, CHAR);
    fileStream << std::endl;
    fileStream << "        </" << prefix << "DataArray>" << std::endl;
  }
  else { fileStream << "/>" << std::endl; }
  fileStream << "      </" << prefix << "Cells>" << std::endl;
}

//***********************************************************************

void OutputXML::writeFinFichierRectilinearXML(std::ofstream &fileStream, bool parallel)
{
  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  if (!parallel) fileStream << "    </Piece>" << std::endl;
  fileStream << "  </" << prefix << "RectilinearGrid>" << std::endl;
  fileStream << "</VTKFile>" << std::endl;
}

//***********************************************************************

void OutputXML::writeFinFichierUnstructuredXML(std::ofstream &fileStream, bool parallel)
{
  std::string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  if (!parallel) fileStream << "    </Piece>" << std::endl;
  fileStream << "  </" << prefix << "UnstructuredGrid>" << std::endl;
  fileStream << "</VTKFile>" << std::endl;
}

//***********************************************************************