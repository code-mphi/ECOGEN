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

#include "OutputVTK.h"
#include "../Run.h"

using namespace tinyxml2;

//***********************************************************************

OutputVTK::OutputVTK(std::string casTest, std::string run, XMLElement* element, std::string fileName, Input* entree) :
  Output(casTest, run, element, fileName, entree)
{
  m_type = TypeOutput::VTK;
}

//***********************************************************************

OutputVTK::OutputVTK(std::string run, int fileNumberRestartMeshMapping, Input* input) : Output(run, fileNumberRestartMeshMapping, input)
{
  m_type = TypeOutput::VTK;
}

//***********************************************************************

OutputVTK::~OutputVTK() {}

//***********************************************************************

void OutputVTK::initializeSpecificOutput()
{
  try {
    std::ofstream fileStream;
    //Set path and extension for collection file in class member and create empty collection file for Paraview
    m_fileCollectionParaview = m_folderOutput + createFilenameVTK(m_filenameCollectionParaview.c_str());
    if (rankCpu == 0) {
      fileStream.open(m_fileCollectionParaview.c_str(), std::ios::trunc);
      if (!fileStream) {
        throw ErrorECOGEN("Unable to open file " + m_fileCollectionParaview, __FILE__, __LINE__);
      }
      fileStream << "<?xml version=\"1.0\"?>" << std::endl;
      fileStream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"";
      if (!m_writeBinary) {
        fileStream << "LittleEndian\" ";
      }
      else {
        fileStream << m_endianMode.c_str() << "\" ";
      }
      fileStream << "compressor=\"vtkZLibDataCompressor\">";
      fileStream << std::endl << "  <Collection>" << std::endl;
      fileStream << std::endl << "  </Collection>" << std::endl;
      fileStream << "</VTKFile>" << std::endl;
      fileStream.close();
      //Set path and extension for collection file in class member and create empty collection file for Visit
      m_fileCollectionVisIt = m_folderOutput + createFilenameVTK(m_filenameCollectionVisIt.c_str(), 0, -1, -1, "visit");
      fileStream.open(m_fileCollectionVisIt.c_str(), std::ios::trunc);
      if (!fileStream) {
        throw ErrorECOGEN("Unable to open file " + m_fileCollectionVisIt, __FILE__, __LINE__);
      }
      fileStream << "!NBLOCKS " << Ncpu << std::endl;
      fileStream.close();
    }
  }
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************

void OutputVTK::writeResults(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
  try {
    //Write output files with VTK format
    writeResultsVTK(mesh, cellsLvl);
    //Add collection file to group levels, time steps, CPUs, etc.
    if (rankCpu == 0) {
      writeCollectionVTK(mesh);
    }
  }
  catch (ErrorECOGEN&) {
    throw;
  } // Go to next level
  m_numFile++;
}

//***********************************************************************

void OutputVTK::readResults(Mesh* mesh, std::vector<Cell*>* cellsLvl, int cpu)
{
  try {
    //1) Parsing XML VTK file
    //-----------------------

    // if 'cpu' arg is not passed or is -1, the method is called in order to resume the computation
    // so we want to open the file saved by the MPI process of current cpu rank. If 'cpu' is provided and is a null or positive integer,
    // we want to restart with mesh mapping:
    int cpuIdx = (cpu == -1) ? rankCpu : cpu;

    std::stringstream fileName(m_folderDatasets + createFilenameVTK(m_fileNameResults.c_str(), mesh, cpuIdx, m_numFile));
    XMLDocument xmlMain;
    XMLError error(xmlMain.LoadFile(fileName.str().c_str())); //File is parsed here
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Entering XML file according to mesh
    //--------------------------------------
    XMLElement *nodeVTK, *nodeGrid, *nodePiece, *nodeCellData;
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
      throw ErrorECOGEN("Output::readResults: unknown mesh type", __FILE__, __LINE__);
      break;
    }
    nodePiece = nodeGrid->FirstChildElement("Piece");
    if (nodePiece == NULL) throw ErrorXMLRacine("Piece", fileName.str(), __FILE__, __LINE__);
    nodeCellData = nodePiece->FirstChildElement("CellData");
    if (nodeCellData == NULL) throw ErrorXMLRacine("CellData", fileName.str(), __FILE__, __LINE__);

    //3) Reading fluid data
    //---------------------
    ReadPhysicalDataVTK(mesh, cellsLvl, nodeCellData, fileName.str());
    if (cpu == -1) m_numFile++;
  } //End try
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************

void OutputVTK::ReadPhysicalDataVTK(Mesh* mesh, std::vector<Cell*>* cellsLvl, XMLElement* nodeCellData, std::string fileName)
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
    for (int phase = 0; phase < m_run->getNumberPhases(); phase++) {
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
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************

std::string OutputVTK::createFilenameVTK(const char* name, Mesh* mesh, int proc, int numFichier, std::string nameVariable)
{
  std::stringstream num;

  try {
    //FP//TODO Donnees separees a faire...
    if (m_splitData) {
      throw ErrorECOGEN("OutputVTK::createFilenameVTK : separated data not done", __FILE__, __LINE__);
    }
    num << name;
    // nameVariable handling
    if ((nameVariable != "defaut") && (nameVariable != "visit")) num << "_" << nameVariable << "_";
    // binary handling
    if (m_writeBinary) num << "B64";
    // cpu handling
    if (proc > -1) num << "_CPU" << proc;
    // Handle number of result files
    if (numFichier != -1) num << "_TIME" << numFichier;
    // Handle extension
    if (mesh == 0) {
      if (nameVariable == "defaut") num << ".pvd";  //Extension for Paraview colleciton
      if (nameVariable == "visit") num << ".visit"; //Extension for VisIt collection
    }
    else {
      switch (mesh->getType()) {
      case REC:
        num << ".vtr";
        break;
      case UNS:
        num << ".vtu";
        break;
      case AMR:
        num << ".vtu";
        break;
        //num << ".vtp"; break;
      default:
        throw ErrorECOGEN("OutputVTK::createFilenameVTK : type mesh unknown", __FILE__, __LINE__);
      }
    }
  }
  catch (ErrorECOGEN&) {
    throw;
  }

  return num.str();
}

//***********************************************************************

void OutputVTK::writeResultsVTK(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
  std::ofstream fileStream;

  try {

    //1) Opening and creation of file
    //-------------------------------
    std::string file = m_folderDatasets + createFilenameVTK(m_fileNameResults.c_str(), mesh, rankCpu, m_numFile);
    fileStream.open(file.c_str(), std::ios::trunc);
    if (!fileStream) {
      throw ErrorECOGEN("Unable to open file " + file, __FILE__, __LINE__);
    }
    fileStream << "<?xml version=\"1.0\"?>" << std::endl;

    //2) Write mesh
    //-------------
    switch (mesh->getType()) {
    case REC:
      writeMeshRectilinearVTK(mesh, fileStream);
      break;
    case UNS:
      writeMeshUnstructuredVTK(mesh, cellsLvl, fileStream);
      break;
    case AMR:
      writeMeshUnstructuredVTK(mesh, cellsLvl, fileStream);
      break;
    default:
      throw ErrorECOGEN("Output::writeResultsVTK : type mesh unknown", __FILE__, __LINE__);
      break;
    }

    //3) Write fluid phase data
    //-------------------------
    writePhysicalDataVTK(mesh, cellsLvl, fileStream);

    //4) Finalize file
    //----------------
    switch (mesh->getType()) {
    case REC:
      writeRectilinearVTKFileEnd(fileStream);
      break;
    case UNS:
      writeUnstructuredVTKFileEnd(fileStream);
      break;
    case AMR:
      writeUnstructuredVTKFileEnd(fileStream);
      break;
    default:
      throw ErrorECOGEN("Output::writeResultsVTK : type mesh unknown", __FILE__, __LINE__);
      break;
    }
    fileStream.close();

  } //End try
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************

void OutputVTK::writeCollectionVTK(Mesh* mesh)
{
  try {
    if (rankCpu == 0) {
      std::ofstream fileStream;
      //ifstream fileStream2((m_folderOutput + m_infoCalcul).c_str()); //For real-time file name
      //double realTime, a, b, c, d, e, f, g, h, i, j, k, l, m;          //For real-time file name
      //double realTime, a, b, c, d;                                     //For real-time file name
      //Creation du file de sortie collection Paraview
      fileStream.open(m_fileCollectionParaview.c_str(), std::ios::trunc);
      if (!fileStream) {
        throw ErrorECOGEN("Unable to open file " + m_fileCollectionParaview, __FILE__, __LINE__);
      }
      fileStream << "<?xml version=\"1.0\"?>" << std::endl;
      fileStream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"";
      if (!m_writeBinary) {
        fileStream << "LittleEndian\" ";
      }
      else {
        fileStream << m_endianMode.c_str() << "\" ";
      }
      fileStream << "compressor=\"vtkZLibDataCompressor\">";
      fileStream << std::endl << "    <Collection>" << std::endl;
      for (int time = 0; time <= m_numFile; time++) {
        //fileStream2 >> a >> b >> realTime >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l >> m; //For real-time file name
        //fileStream2 >> a >> b >> realTime >> c >> d;                                              //For real-time file name
        for (int p = 0; p < Ncpu; p++) {
          std::string file = "datasets/" + createFilenameVTK(m_fileNameResults.c_str(), mesh, p, time);
          fileStream << "        <DataSet timestep=\"" << time << "\" part=\"" << p << "\" file=\"" << file.c_str() << "\"/>" << std::endl;
          //fileStream << "        <DataSet timestep=\"" << realTime << "\" part=\"" << p << "\" file=\"" << file.c_str() << "\"/>" << std::endl; //For real-time file name
        }
      }
      fileStream << "    </Collection>" << std::endl;
      fileStream << "</VTKFile>" << std::endl;
      fileStream.close();
      //fileStream2.close(); //For real-time file name

      // Create output collection file for VisIt
      fileStream.open(m_fileCollectionVisIt.c_str(), std::ios::trunc);
      if (!fileStream) {
        throw ErrorECOGEN("Unable to open file " + m_fileCollectionVisIt, __FILE__, __LINE__);
      }
      fileStream << "!NBLOCKS " << Ncpu << std::endl;
      for (int time = 0; time <= m_numFile; time++) {
        for (int p = 0; p < Ncpu; p++) {
          std::string file = "datasets/" + createFilenameVTK(m_fileNameResults.c_str(), mesh, p, time);
          fileStream << file.c_str() << std::endl;
        }
      }
      fileStream.close();
    }
  }
  catch (ErrorXML&) {
    throw;
  } // Go to next level
}

//***********************************************************************

void OutputVTK::writePhysicalDataVTK(Mesh* mesh, std::vector<Cell*>* cellsLvl, std::ofstream& fileStream)
{
  std::vector<double> dataset;

  std::string format = "ascii";
  if (m_writeBinary) format = "binary";

  fileStream << "      <CellData>" << std::endl;

  //1) Write phases variables
  //---------------------------------
  for (int phase = 0; phase < m_run->getNumberPhases(); phase++) //For complete output
  //for (int phase = 0; phase < 1; phase++) //For reduced output
  {
    std::string eosName;
    if (m_cellRef.getPhase(phase)->getEos() != nullptr) {
      eosName = "_" + m_cellRef.getPhase(phase)->getEos()->getName();
      eosName.erase(eosName.end() - 4, eosName.end());
    }
    //Write scalar variables
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberScalars(); var++) {
      fileStream << "        <DataArray type=\"Float64\" Name=\"" + outputPhaseName(phase, m_cellRef.getPhase(phase)->returnNameScalar(var), eosName);
      fileStream << "\"";
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->getData(cellsLvl, dataset, var, phase);
      this->writeDataset(dataset, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </DataArray>" << std::endl;
    }
    //Write vectorial variables
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++) {
      fileStream << "        <DataArray type=\"Float64\" Name=\"" + outputPhaseName(phase, m_cellRef.getPhase(phase)->returnNameVector(var), eosName);
      fileStream << "\" NumberOfComponents=\"3\"";
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->getData(cellsLvl, dataset, -var, phase);
      this->writeDataset(dataset, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </DataArray>" << std::endl;
    }
  } //End phase

  //2) Write mixture data
  //---------------------
  if (m_run->m_numberPhases > 1) {
    int mixture = -1;
    //Write scalars mixture variables
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++) {
      fileStream << "        <DataArray type=\"Float64\" Name=\"" << m_cellRef.getMixture()->returnNameScalar(var) << "\"";
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->getData(cellsLvl, dataset, var, mixture);
      this->writeDataset(dataset, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </DataArray>" << std::endl;
    }
    //Write vectorial mixture variables
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++) {
      fileStream << "        <DataArray type=\"Float64\" Name=\"" << m_cellRef.getMixture()->returnNameVector(var) << "\" NumberOfComponents=\"3\"";
      fileStream << " format=\"" << format << "\">" << std::endl;
      mesh->getData(cellsLvl, dataset, -var, mixture);
      this->writeDataset(dataset, fileStream, DOUBLE);
      fileStream << std::endl;
      fileStream << "        </DataArray>" << std::endl;
    }
  } //End mixture

  //3) Write of transports and others...
  //------------------------------------
  int transport = -2; //For complete output
  for (int var = 1; var <= m_run->m_numberTransports; var++) {
    fileStream << "        <DataArray type=\"Float64\" Name=\"T" << var << "\"";
    fileStream << " format=\"" << format << "\">" << std::endl;
    mesh->getData(cellsLvl, dataset, var, transport);
    this->writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </DataArray>" << std::endl;
  }

  //4) Write xi indicator
  //----------------------
  if (mesh->getType() == AMR) { //For complete output
    int xi = -3;
    fileStream << "        <DataArray type=\"Float64\" Name=\"Xi\"";
    fileStream << " format=\"" << format << "\">" << std::endl;
    mesh->getData(cellsLvl, dataset, 1, xi);
    this->writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </DataArray>" << std::endl;
  }

  //5) Write gradient rho
  //---------------------
  int gradRho_ = -4;
  fileStream << "        <DataArray type=\"Float64\" Name=\"gradRho\"";
  fileStream << " format=\"" << format << "\">" << std::endl;
  mesh->getData(cellsLvl, dataset, 1, gradRho_);
  this->writeDataset(dataset, fileStream, DOUBLE);
  fileStream << std::endl;
  fileStream << "        </DataArray>" << std::endl;

  //6) Absolute velocity printing for Moving Reference Frame computations
  //---------------------------------------------------------------------
  if (m_run->m_MRF != -1) {
    fileStream << "        <DataArray type=\"Float64\" Name=\"absoluteVelocityMRF\" NumberOfComponents=\"3\"";
    fileStream << " format=\"" << format << "\">" << std::endl;
    mesh->extractAbsVelocityMRF(cellsLvl, dataset, m_run->m_sources[m_run->m_MRF]);
    this->writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </DataArray>" << std::endl;
  }

  //7) Cells' reference length
  //--------------------------
  if (m_run->m_extractRefLength) {
    fileStream << "        <DataArray type=\"Float64\" Name=\"Reference_length\"";
    fileStream << " format=\"" << format << "\">" << std::endl;
    mesh->extractReferenceLength(cellsLvl, dataset);
    this->writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </DataArray>" << std::endl;
  }

  //8) CPU rank
  //-----------
  // if (mesh->getType() == AMR) { //For complete output
  //   int CPUrank = -5;
  //   fileStream << "        <DataArray type=\"Int32\" Name=\"CPUrank\"";
  //   fileStream << " format=\"" << format << "\">" << std::endl;
  //   mesh->getData(cellsLvl, dataset, 1, CPUrank);
  //   this->writeDataset(dataset, fileStream, INT);
  //   fileStream << std::endl;
  //   fileStream << "        </DataArray>" << std::endl;
  // }

  //9) Morton index
  //---------------
  // if (mesh->getType() == AMR) { //For complete output
  //   int MortonIndex = -6;
  //   fileStream << "        <DataArray type=\"Int32\" Name=\"MortonIndex\"";
  //   fileStream << " format=\"" << format << "\">" << std::endl;
  //   mesh->getData(cellsLvl, dataset, 1, MortonIndex);
  //   this->writeDataset(dataset, fileStream, INT);
  //   fileStream << std::endl;
  //   fileStream << "        </DataArray>" << std::endl;
  // }

  //9) Saturation pressure (specific recording)
  //-------------------------------------------
  if (m_run->m_recordPsat) {
    fileStream << "        <DataArray type=\"Float64\" Name=\"Psat\"";
    fileStream << " format=\"" << format << "\">" << std::endl;
    mesh->getData(cellsLvl, dataset, 1, -7);
    this->writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </DataArray>" << std::endl;
  }

  //10) Turbulent viscosity (nu_T)
  //------------------------------
  if (turbulenceActivated && !m_reducedOutput) {
    fileStream << "        <DataArray type=\"Float64\" Name=\"nu_t\"";
    fileStream << " format=\"" << format << "\">" << std::endl;
    mesh->getData(cellsLvl, dataset, 1, -8);
    this->writeDataset(dataset, fileStream, DOUBLE);
    fileStream << std::endl;
    fileStream << "        </DataArray>" << std::endl;
  }

  //Fin
  fileStream << "      </CellData>" << std::endl;
}

//***********************************************************************

void OutputVTK::writeMeshRectilinearVTK(Mesh* mesh, std::ofstream& fileStream)
{
  std::vector<double> dataset;

  //0) Header
  //---------
  fileStream << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"";
  if (!m_writeBinary) fileStream << "LittleEndian\">" << std::endl;
  else fileStream << m_endianMode.c_str() << "\">" << std::endl;
  fileStream << "  <RectilinearGrid WholeExtent=\"" << mesh->getStringExtent() << "\">" << std::endl;
  fileStream << "    <Piece Extent=\"" << mesh->getStringExtent() << "\">" << std::endl;

  //1) Write node coordinates
  //-------------------------
  fileStream << "      <Coordinates>" << std::endl;
  //Coordinate in X
  fileStream << "        <DataArray type=\"Float64\" ";
  if (!m_writeBinary) {
    fileStream << "format=\"ascii\">" << std::endl << "          ";
  }
  else {
    fileStream << "format=\"binary\">" << std::endl;
  }

  dataset.clear();
  mesh->getCoord(dataset, X);
  writeDataset(dataset, fileStream, DOUBLE);
  fileStream << std::endl;
  fileStream << "        </DataArray>" << std::endl;
  //Coordinate in Y
  fileStream << "        <DataArray type=\"Float64\" ";
  if (!m_writeBinary) {
    fileStream << "format=\"ascii\">" << std::endl << "          ";
  }
  else {
    fileStream << "format=\"binary\">" << std::endl;
  }
  dataset.clear();
  mesh->getCoord(dataset, Y);
  writeDataset(dataset, fileStream, DOUBLE);
  fileStream << std::endl;
  fileStream << "        </DataArray>" << std::endl;
  //Coordinate in Z
  fileStream << "        <DataArray type=\"Float64\" ";
  if (!m_writeBinary) {
    fileStream << "format=\"ascii\">" << std::endl << "          ";
  }
  else {
    fileStream << "format=\"binary\">" << std::endl;
  }
  dataset.clear();
  mesh->getCoord(dataset, Z);
  writeDataset(dataset, fileStream, DOUBLE);
  fileStream << std::endl;
  fileStream << "        </DataArray>" << std::endl;
  fileStream << "      </Coordinates>" << std::endl;
}

//***********************************************************************

void OutputVTK::writeMeshUnstructuredVTK(Mesh* mesh, std::vector<Cell*>* cellsLvl, std::ofstream& fileStream)
{
  std::vector<double> dataset;

  //0) Header
  //---------
  fileStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"";
  if (!m_writeBinary) fileStream << "LittleEndian\">" << std::endl;
  else fileStream << m_endianMode.c_str() << "\">" << std::endl;

  fileStream << "  <UnstructuredGrid>" << std::endl;
  mesh->writeHeaderPiece(fileStream, cellsLvl);

  //1) Nodes output
  //-----------------
  fileStream << "      <Points>" << std::endl;
  fileStream << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" ";
  if (!m_writeBinary) {
    fileStream << "format=\"ascii\">" << std::endl << "          ";
  }
  else {
    fileStream << "format=\"binary\">" << std::endl;
  }
  dataset.clear();
  size_t numberNodes{mesh->getNodes(dataset, cellsLvl)};
  writeDataset(dataset, fileStream, DOUBLE);
  fileStream << std::endl;
  fileStream << "        </DataArray>" << std::endl;
  fileStream << "      </Points>" << std::endl;

  //2) Cells output
  //------------------
  fileStream << "      <Cells>" << std::endl;
  //Connectivity
  fileStream << "        <DataArray type=\"Int32\" Name=\"connectivity\" ";
  if (!m_writeBinary) {
    fileStream << "format=\"ascii\">" << std::endl << "          ";
  }
  else {
    fileStream << "format=\"binary\">" << std::endl;
  }

  dataset.clear();
  if (numberNodes > MAX_SAFE_INTEGER) { // Check that the connectivity array can be exactly casted into a double array
    throw ErrorECOGEN("Too many nodes to call current implementation of '" + (std::string) __func__ +
                      "' method.\n"
                      "Consider changing the 'ECOGEN_VTKHDF_intSize' value at CMake configuration step.");
  }
  mesh->getConnectivity(dataset, cellsLvl);
  writeDataset(dataset, fileStream, INT);
  fileStream << std::endl;
  fileStream << "        </DataArray>" << std::endl;
  //Offsets
  fileStream << "        <DataArray type=\"Int32\" Name=\"offsets\" ";
  if (!m_writeBinary) {
    fileStream << "format=\"ascii\">" << std::endl << "          ";
  }
  else {
    fileStream << "format=\"binary\">" << std::endl;
  }
  dataset.clear();
  mesh->getOffsets(dataset, cellsLvl);
  writeDataset(dataset, fileStream, INT);
  fileStream << std::endl;
  fileStream << "        </DataArray>" << std::endl;
  //Cells type
  fileStream << "        <DataArray type=\"UInt8\" Name=\"types\" ";
  if (!m_writeBinary) {
    fileStream << "format=\"ascii\">" << std::endl << "          ";
  }
  else {
    fileStream << "format=\"binary\">" << std::endl;
  }
  dataset.clear();
  mesh->getTypeCell(dataset, cellsLvl);
  writeDataset(dataset, fileStream, CHAR);
  fileStream << std::endl;
  fileStream << "        </DataArray>" << std::endl;
  fileStream << "      </Cells>" << std::endl;
}

//***********************************************************************

void OutputVTK::writeRectilinearVTKFileEnd(std::ofstream& fileStream)
{
  fileStream << "    </Piece>" << std::endl;
  fileStream << "  </RectilinearGrid>" << std::endl;
  fileStream << "</VTKFile>" << std::endl;
}

//***********************************************************************

void OutputVTK::writeUnstructuredVTKFileEnd(std::ofstream& fileStream)
{
  fileStream << "    </Piece>" << std::endl;
  fileStream << "  </UnstructuredGrid>" << std::endl;
  fileStream << "</VTKFile>" << std::endl;
}

//***********************************************************************
