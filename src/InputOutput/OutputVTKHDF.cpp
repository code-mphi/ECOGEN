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

#include <H5Tpublic.h>
#include <cstdlib>
#include <unistd.h>
#include <mpi.h>
#include <limits>
#include "OutputVTKHDF.h"
#include "../Run.h"
#include "../Colors.h"

#include "IntTypes.h" // Define preprocessor variables and aliases to appropriate integer types depending on architecture and compilation choice

constexpr int OutputVTKHDF::VTKHDF_VERSION[2];

using namespace tinyxml2;

//***********************************************************************

OutputVTKHDF::OutputVTKHDF(std::string casTest, std::string run, XMLElement* element, std::string fileName, Input* entree) :
  Output(casTest, run, element, fileName, entree)
{
  m_h5_file_floatType = H5T_IEEE_F32LE; // default value is simple precision, i.e a decimal presision of 8
  m_h5_file_intType   = H5T_STD_I32LE;  // default value is 32-bit integers

  m_h5_mem_floatType = H5T_NATIVE_FLOAT;
  m_h5_mem_intType   = H5T_NATIVE_INT32;

  int intSize{32}; // Default integer size if not specified

  element->QueryIntAttribute("intSize", &intSize);

  if (intSize == 32) {
    m_h5_file_intType = H5T_STD_I32LE; // 32-bit integers
    m_h5_mem_intType  = H5T_NATIVE_INT32;
  }
  else if (intSize == 64) {
    m_h5_file_intType = H5T_STD_I64LE; // 64-bit integers
    m_h5_mem_intType  = H5T_NATIVE_INT64;
  }
  else { // Unexpected value
    throw ErrorXMLMessage("VTKHDF file format only compatible with 32 or 64-bits integers.\n"
                          "          Please, modify the 'intSize' attribute of your input file.");
  }

  // Check precision value provided by user: only 8 (for Float32) and 16 (for Float64) values are autorized
  if (m_precision != 0) {
    if (m_precision == 16) {
      m_h5_file_floatType = H5T_IEEE_F64LE;
      m_h5_mem_floatType  = H5T_NATIVE_DOUBLE;
    }
    else if (m_precision == 8) {
      m_h5_file_floatType = H5T_IEEE_F32LE;
      m_h5_mem_floatType  = H5T_NATIVE_FLOAT;
    }
    else {
      throw ErrorXMLMessage("VTKHDF file format only compatible with decimal precision value of 8 or 16.\n"
                            "          Please, modify the 'precision' attribute of your input file.");
    }
  }

  if (!m_writeBinary) {
    throw ErrorXMLMessage("VTKHDF file format only compatible with binary outputs.\n"
                          "          Please, modify the 'binary' attribute of your input file.");
  }

  // Set parallel mode for raw data transfer
  const char* transferMode_c;
  transferMode_c = element->Attribute("transferMode"); // Get string defining transfer mode

  H5FD_mpio_xfer_t transfer_mode;
  if (transferMode_c == nullptr) {
    transfer_mode = H5FD_MPIO_COLLECTIVE; // default mode when user doesn't specify what he wants
  }
  else {
    std::string transferMode_str{transferMode_c};
    Tools::uppercase(transferMode_str);

    if (transferMode_str == "H5FD_MPIO_COLLECTIVE") {
      transfer_mode = H5FD_MPIO_COLLECTIVE;
    }
    else if (transferMode_str == "H5FD_MPIO_INDEPENDENT") {
      transfer_mode = H5FD_MPIO_INDEPENDENT;
    }
    else {
      throw ErrorXMLMessage("VTKHDF file format only compatible with 'H5FD_MPIO_INDEPENDENT' or"
                            " 'H5FD_MPIO_COLLECTIVE' transfer mode.");
    }
  }

  // Initialization of output type and hdf5 properties.
  this->initializeProperties(transfer_mode);
}

//***********************************************************************

OutputVTKHDF::OutputVTKHDF(std::string run, int fileNumberRestartMeshMapping, Input* input, const H5FD_mpio_xfer_t& transfer_mode) :
  Output(run, fileNumberRestartMeshMapping, input)
{
  // This constructor is called when restarting a simulation with mesh mapping: an incomplete output object is created
  // in order to read the output results of the previous simulation. Most of output parameters readed in 'main.xml' file are not
  // needed (such as for example, the output precision, its ASCII/binary type, etc.). In case of hdf5 outputs, the collective or independent
  // transfer mode (that is readed in the 'main.xml' file) is needed. In order to allow OutputVTKHDF object instanciation
  // without reading again this file, the transfer mode is recovered from previously constructed 'outputVTKHDF' object, passed by argument to this constructor
  // and initialized by 'initializeProperties' method.
  this->initializeProperties(transfer_mode);
}

//***********************************************************************

OutputVTKHDF::~OutputVTKHDF()
{
  if (m_group_creation_properties >= 0) H5Pclose(m_group_creation_properties);
  if (m_dataset_creation_properties >= 0) H5Pclose(m_dataset_creation_properties);
  if (m_parallel_x_properties >= 0) H5Pclose(m_parallel_x_properties);
}

//***********************************************************************

void OutputVTKHDF::initializeProperties(const H5FD_mpio_xfer_t& transfer_mode)
{
  // -- Properties needed when writing and reading output files.

  // Output type
  m_type = TypeOutput::VTKHDF;

  // Hdf5 properties
  m_group_creation_properties   = H5Pcreate(H5P_GROUP_CREATE);
  m_dataset_creation_properties = H5Pcreate(H5P_DATASET_CREATE);
  m_parallel_x_properties       = H5Pcreate(H5P_DATASET_XFER);

  if (m_group_creation_properties < 0 || m_parallel_x_properties < 0 || m_dataset_creation_properties < 0) {
    throw ErrorECOGEN("Unable to create hdf5 variables", __FILE__, __LINE__);
  }

  // Enable attribute tracking and indexing of creation order and fasten attribute traversal and metadata analysis (not sure that it is useful)
  H5Pset_attr_creation_order(m_group_creation_properties, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
  H5Pset_fill_time(m_dataset_creation_properties, H5D_FILL_TIME_NEVER); // do not initialize dataset values

  // Assign m_transfer_mode and set I/Os for raw data transfer to, either collective, or independent
  m_transfer_mode = transfer_mode;
  assert((m_transfer_mode == H5FD_MPIO_COLLECTIVE || m_transfer_mode == H5FD_MPIO_INDEPENDENT) && "Wrong transfer mode for HDF5");
  H5Pset_dxpl_mpio(m_parallel_x_properties, m_transfer_mode);
}

//***********************************************************************

void OutputVTKHDF::initializeSpecificOutput()
{
  try {
    std::ofstream fileStream;
    //Add path for timeseries output file to m_fileVTKHDFTimeseries class member
    m_filenameVTKHDFTimeseries = m_folderOutput + m_filenameVTKHDFTimeseries.c_str();
    //Create an empty file to ensure to erase possible pre-existing time serie file.
    //This empty file will be erased if we succeed to save the timestep 0.
    if (rankCpu == 0) {
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
      if (m_run->m_mesh->getType() == REC) {
        // Temporary solution to save output file at cartesian format while as long as vtkhdf format is not yet available for rectilinear meshes
        if (m_run->m_resumeSimulation == 0) { // Create a new empty file and fill header
          writeXDMF(nullptr, -1);
        }
      }
      else {
#endif
        // Create a timeserie file (possible improvement: save all timesteps in one file only and static mesh only once)
        writeTimeseriesVTKHDF(-1);
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
      }
#endif
    }
  }
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************

void OutputVTKHDF::writeResults(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
  try {
    //Write output files with VTKHDF format
    writeResultsVTKHDF(mesh, cellsLvl);
    //Add timeseries file to group time steps.
    if (rankCpu == 0) {
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
      if (mesh->getType() == REC) {
        // Temporary solution to save output file at cartesian format while as long as vtkhdf format is not yet available for rectilinear meshes
        writeXDMF(mesh, m_numFile);
      }
      else {
#endif
        // Create a timeserie file (possible improvement: save all timesteps in one file only and static mesh only once)
        writeTimeseriesVTKHDF(m_numFile);
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
      }
#endif
    }
  }
  catch (ErrorECOGEN&) {
    throw;
  } // Go to next level
  m_numFile++;
}

//***********************************************************************

int OutputVTKHDF::readScalarVTKHDF(hid_t group_id, const std::string& dataset_name, std::vector<size_t>& scalar_array) const
{
  hid_t dataspace{-1};
  hid_t hdf_dataset{-1};

  int status{SUCCESS_VALUE};

  hsize_t dimsf[1];

  // Read the wanted dataset
  hdf_dataset = H5Dopen2(group_id, dataset_name.c_str(), H5P_DEFAULT);

  if (hdf_dataset < 0) {
    Errors::errorMessage("Missing '" + dataset_name + "' dataset in VTKHDF file.", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  if (status != ERROR_VALUE) {
    dataspace = H5Dget_space(hdf_dataset);
    if (dataspace < 0) {
      Errors::errorMessage("Unable to get dataspace associated to '" + dataset_name + "' dataset of VTKHDF file.", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  if (status != ERROR_VALUE) {
    // Get dataset dimensions: in case of mesh mapping, the number of MPI processes may have changed so we cannot assume than
    // dimsf == Ncpu
    H5Sget_simple_extent_dims(dataspace, dimsf, NULL);
  }

  // Read the complete array on each rank (as it provide number of entities saved by each mpi process when writting the file,
  // we will need the complete array to compute the offset of this process and read the wanted part of the data arrays.
  if (status != ERROR_VALUE) {
    scalar_array.resize(dimsf[0]);
    herr_t status_read = H5Dread(hdf_dataset, ECOGEN_H5_SIZE_T, H5S_ALL, H5S_ALL, m_parallel_x_properties, scalar_array.data());
    if (status_read < 0) {
      Errors::errorMessage("Error when reading '" + dataset_name + "' dataset in VTKHDF file.", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  // Close hdf5 objects
  if (dataspace >= 0) H5Sclose(dataspace);
  if (hdf_dataset >= 0) H5Dclose(hdf_dataset);

  return status;
}

//***********************************************************************
template <size_t N>
int OutputVTKHDF::readArrayVTKHDF(
  hid_t group_id, std::string dataset_name, std::vector<double>& dataset, const hsize_t (&local_dims)[N], const hsize_t (&offsets)[N])
{

  int status{SUCCESS_VALUE};

  hid_t hdf_dataset{H5Dopen2(group_id, dataset_name.c_str(), H5P_DEFAULT)};
  if (hdf_dataset < 0) {
    Errors::errorMessage("Missing '" + dataset_name + "' dataset in VTKHDF file.", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  hid_t dataspace{-1};
  if (status != ERROR_VALUE) {
    dataspace = H5Dget_space(hdf_dataset);
    if (dataspace < 0) {
      Errors::errorMessage("Unable to get dataspace associated to '" + dataset_name + "' dataset of VTKHDF file.", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  if (status != ERROR_VALUE && H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, local_dims, NULL) < 0) {
    Errors::errorMessage("Unable to create hdf5 hyperslab to read '" + dataset_name + "' dataset.", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Create a memspace to store array of size line_count x column_count
  hid_t memspace{H5Screate_simple(N, local_dims, NULL)};
  if (dataspace < 0) {
    Errors::errorMessage("Unable to create hdf5 memspace", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Read value
  if (status != ERROR_VALUE && H5Dread(hdf_dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, m_parallel_x_properties, dataset.data()) < 0) {
    Errors::errorMessage("Error when reading '" + dataset_name + "' dataset in VTKHDF file.", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Close hdf5 objects
  if (dataspace >= 0) H5Sclose(dataspace);
  if (memspace >= 0) H5Sclose(memspace);
  if (hdf_dataset >= 0) H5Dclose(hdf_dataset);

  return status;
}

//***********************************************************************

void OutputVTKHDF::readResults(Mesh* mesh, std::vector<Cell*>* cellsLvl, int cpu)
{
  int status = SUCCESS_VALUE;

  // if 'cpu' arg is not passed or is -1, the method is called in order to resume the computation
  // so we want to read the part of file saved by the MPI process of current cpu rank. If 'cpu' is provided and is a null or positive integer,
  // we want to restart with mesh mapping and to read the part of file saved by another cpu rank.
  const int cpuIdx = (cpu == -1) ? rankCpu : cpu;

  //1) Parsing VTKHDF file hierarchy
  //--------------------------------
  std::string filename = m_folderDatasets + createFilenameVTKHDF(m_fileNameResults.c_str(), m_numFile);

  hid_t plist_id{H5Pcreate(H5P_FILE_ACCESS)}; // define a variable to store file access properties
  if (plist_id < 0) {
    throw ErrorECOGEN("Unable to create hdf5 file access property variable", __FILE__, __LINE__);
  }
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL); // store MPI IO communicator infos to plist_id
  H5Pset_all_coll_metadata_ops(plist_id, 1);                 // enable MPI-IO collective read operation for metadata
                                                             // (to be used when size of metadata created by the app is large)

  hid_t file_id{H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id)};

  if (file_id < 0) {
    Errors::errorMessage("Unable to open '" + filename + "' file.", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  hid_t vtkhdf_group_id{-1};
  if (status != ERROR_VALUE) {
    vtkhdf_group_id = H5Gopen2(file_id, "VTKHDF", H5P_DEFAULT);
    if (vtkhdf_group_id < 0) {
      Errors::errorMessage("Unable to find 'vtkhdf' group in hierarchy of '" + filename + "' file.", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  // For unstructured and AMR meshes, check that the file contains the expected cell count.
  std::vector<size_t> numberCells_r;
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
  if (mesh->getType() != REC) {
#endif
    // Security check: count the local number of cells and check that it maches the cell count saved in VTKHDF file
    if (status != ERROR_VALUE) {
      status = readScalarVTKHDF(vtkhdf_group_id, "NumberOfCells", numberCells_r);
    }

    if (status != ERROR_VALUE) {
      // count the local number of cells and compare it with the value that has been readed
      size_t numberCellsLocal{0};
      if (mesh->getLvlMax() == 0) { // Non AMR case: get the number of cells
        numberCellsLocal = mesh->getNumberCells();
      }
      else {
        for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
          for (size_t i = 0; i < cellsLvl[lvl].size(); ++i) {
            if (!cellsLvl[lvl][i]->getSplit()) { // Do not count Cells with children
              ++numberCellsLocal;
            }
          }
        }
      }
      if (numberCellsLocal != numberCells_r[cpuIdx]) {
        Errors::errorMessage("Unexpected: local cell count mismatches NumberOfCells for MPI rank " + std::to_string(cpuIdx) + " in '" + filename +
                               "' file.",
                             __FILE__,
                             __LINE__);
        status = ERROR_VALUE;
      }
    }
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
  }
#endif

  //2) Reading fluid data
  //---------------------
  hid_t celldata_group_id{-1};
  if (status != ERROR_VALUE) {
    celldata_group_id = H5Gopen2(vtkhdf_group_id, "CellData", H5P_DEFAULT);
    if (celldata_group_id < 0) {
      Errors::errorMessage("Unable to find 'CellData' group in hierarchy of '" + filename + "' file.", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  if (status != ERROR_VALUE) status = ReadPhysicalDataVTKHDF(mesh, cellsLvl, celldata_group_id, cpuIdx, numberCells_r);

  //4) Finalize file
  //----------------
  // Close all opened HDF5 objects
  if (celldata_group_id >= 0) H5Gclose(celldata_group_id);
  if (vtkhdf_group_id >= 0) H5Gclose(vtkhdf_group_id);
  if (plist_id >= 0) H5Pclose(plist_id);

  herr_t status_close = H5Fclose(file_id);
  if (status_close < 0) {
    // We try to close the file but it fails... nothing more to do automatically
    throw ErrorECOGEN("Error when closing '" + filename + "' file.", __FILE__, __LINE__);
  }
  if (status == ERROR_VALUE) {
    // From here we can throw errors because the file has been closed cleanly
    throw ErrorECOGEN("Error when reading '" + filename + "' file.", __FILE__, __LINE__);
  }

  if (cpu == -1) m_numFile++;
}

//***********************************************************************

template <size_t N>
int OutputVTKHDF::ReadPhysicalDataVTKHDF(
  Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t group_id, const hsize_t (&count_scalar)[N], const hsize_t (&start_scalar)[N])
{
  int status = SUCCESS_VALUE;

  // Create arrays storing dimensions and offset for vectorial data parsing: dimensions will be the
  // same than for scalar data parsing plus one additional range of dimension DIM_OUTPUT to store
  // the vector fields
  hsize_t count_vector[N + 1];
  hsize_t start_vector[N + 1];
  for (size_t i = 0; i < N; ++i) {
    count_vector[i] = count_scalar[i];
    start_vector[i] = start_scalar[i];
  }
  count_vector[N] = DIM_OUTPUT;
  start_vector[N] = 0;

  // Allocate vector to store data to read
  size_t count{count_scalar[0]};
  for (size_t i = 1; i < N; ++i) {
    count *= count_scalar[i];
  }
  std::vector<double> scalarDataSet(count);
  std::vector<double> vectorDataSet(DIM_OUTPUT * count);

  //1) Reading phases data
  //----------------------
  for (int phase = 0; phase < m_run->getNumberPhases(); phase++) {

    std::string eosName; //empty string
    if (m_cellRef.getPhase(phase)->getEos() != nullptr) {
      eosName = "_" + m_cellRef.getPhase(phase)->getEos()->getName();
      eosName.erase(eosName.end() - 4, eosName.end());
    }

    //Reading scalars
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberScalars(); var++) {
      if (status != ERROR_VALUE) {
        std::string dataName{outputPhaseName(phase, m_cellRef.getPhase(phase)->returnNameScalar(var), eosName)};
        status = readArrayVTKHDF(group_id, dataName, scalarDataSet, count_scalar, start_scalar);
        mesh->setDataSet(scalarDataSet, cellsLvl, var, phase);
      }
    }
    //Reading vectors
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++) {
      if (status != ERROR_VALUE) {
        std::string dataName{outputPhaseName(phase, m_cellRef.getPhase(phase)->returnNameVector(var), eosName)};
        status = readArrayVTKHDF(group_id, dataName, vectorDataSet, count_vector, start_vector);
        mesh->setDataSet(vectorDataSet, cellsLvl, -var, phase);
      }
    }
  } //End phase

  //2) Reading mixture data
  //-----------------------
  if (m_run->m_numberPhases > 1) {
    int mixture = -1;
    //Reading scalars
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++) {
      if (status != ERROR_VALUE) {
        std::string dataName{m_cellRef.getMixture()->returnNameScalar(var)};
        status = readArrayVTKHDF(group_id, dataName, scalarDataSet, count_scalar, start_scalar);
        mesh->setDataSet(scalarDataSet, cellsLvl, var, mixture);
      }
    }
    //Reading vectors
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++) {
      if (status != ERROR_VALUE) {
        std::string dataName{m_cellRef.getMixture()->returnNameVector(var)};
        status = readArrayVTKHDF(group_id, dataName, vectorDataSet, count_vector, start_vector);
        mesh->setDataSet(vectorDataSet, cellsLvl, -var, mixture);
      }
    }
  } //End mixture

  //3) Transported data
  //-------------------
  int transport = -2;
  for (int var = 1; var <= m_run->m_numberTransports; var++) {
    if (status != ERROR_VALUE) {
      std::string dataName{"T" + std::to_string(var)};
      status = readArrayVTKHDF(group_id, dataName, scalarDataSet, count_scalar, start_scalar);
      mesh->setDataSet(scalarDataSet, cellsLvl, var, transport);
    }
  }

  //4) xi indicator for AMR
  //-----------------------
  if (mesh->getType() == AMR) {
    int xi = -3;
    if (status != ERROR_VALUE) {
      status = readArrayVTKHDF(group_id, "Xi", scalarDataSet, count_scalar, start_scalar);
      mesh->setDataSet(scalarDataSet, cellsLvl, 1, xi);
    }
  }

  return status;
}

//***********************************************************************

int OutputVTKHDF::ReadPhysicalDataVTKHDF(Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t group_id, int cpu, const std::vector<size_t>& numberCells)
{
  int status = SUCCESS_VALUE;

  switch (mesh->getType()) {
  case REC:
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
  {
    // Each rank reads in a location that matches the geometrical representation of mesh partitionning
    auto numberCellsX{mesh->getNumberCellsX()};
    auto numberCellsY{mesh->getNumberCellsY()};
    auto numberCellsZ{mesh->getNumberCellsZ()};

    auto offsetX{mesh->getOffsetX()};
    auto offsetY{mesh->getOffsetY()};
    auto offsetZ{mesh->getOffsetZ()};

    hsize_t count_scalar[DIM_OUTPUT] = {(hsize_t)numberCellsZ, (hsize_t)numberCellsY, (hsize_t)numberCellsX};
    hsize_t start_scalar[DIM_OUTPUT] = {(hsize_t)offsetZ, (hsize_t)offsetY, (hsize_t)offsetX};

    status = ReadPhysicalDataVTKHDF(mesh, cellsLvl, group_id, count_scalar, start_scalar);
    break;
  }
#endif
  case UNS:
  case AMR:
    // Always deal with 2-dimensional arrays. In case of scalars, second dimension size will be 1.
    hsize_t count_scalar[1];
    hsize_t start_scalar[1];

    start_scalar[0] = 0;
    for (int i = 0; i < cpu; ++i) {
      start_scalar[0] += numberCells[i];
    }
    count_scalar[0] = numberCells[cpu];

    status = ReadPhysicalDataVTKHDF(mesh, cellsLvl, group_id, count_scalar, start_scalar);

    break;
  default:

    Errors::errorMessage("Output::readPhysicalDataVTKHDF: type mesh unknown", __FILE__, __LINE__);
    break;
  }

  return status;
}

//***********************************************************************

std::string OutputVTKHDF::createFilenameVTKHDF(const char* name, int numFile, std::string nameVariable) const
{
  std::stringstream num;

  try {
    //FP//TODO Donnees separees a faire...
    if (m_splitData) {
      throw ErrorECOGEN("OutputVTKHDF::createFilenameVTKHDF : separated data not done", __FILE__, __LINE__);
    }
    num << name;

    // nameVariable handling
    if (nameVariable != "defaut") num << "_" << nameVariable;
    // Handle number of result files
    if (numFile != -1) num << "_TIME" << numFile << ".vtkhdf";
  }
  catch (ErrorECOGEN&) {
    throw;
  }

  return num.str();
}

//***********************************************************************

int OutputVTKHDF::openFileVTKHDF(hid_t& file_id, const bool parallel_access, const bool create) const
{
  int status = SUCCESS_VALUE;

  std::string filename = m_folderDatasets + createFilenameVTKHDF(m_fileNameResults.c_str(), m_numFile);

  hid_t plist_id{H5Pcreate(H5P_FILE_ACCESS)}; // define a variable to store file access properties
  if (plist_id < 0) {
    throw ErrorECOGEN("Unable to create hdf5 file access property variable", __FILE__, __LINE__);
  }

  if (parallel_access) {
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL); // store MPI IO communicator infos to plist_id
    H5Pset_coll_metadata_write(plist_id, 1);                   // enable MPI-IO collective write operation for metadata
                                                               // (to be used when size of metadata created by the app is large)
    H5Pset_all_coll_metadata_ops(plist_id, 1);                 // enable MPI-IO collective read operation for metadata
                                                               // (to be used when size of metadata created by the app is large)
  }

  if (create) { // File has to be created
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  }
  else { // File already exists and has to be re-opened for completion (rectilinear case)
    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
  }

  if (file_id < 0) {
    Errors::errorMessage("Unable to create and open hdf5 file '" + filename + "'.", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Close file access properties object
  if (plist_id >= 0) H5Pclose(plist_id);

  return status;
}

//***********************************************************************

void OutputVTKHDF::writeResultsVTKHDF(Mesh* mesh, std::vector<Cell*>* cellsLvl)
{
  int status = SUCCESS_VALUE;

  hid_t file_id{-1};

  //1) Open file and write mesh: file handler is stored in 'file_id' so file can be closed later
  //--------------------------------------------------------------------------------------------
  size_t numberCells{0};
  hid_t vtkhdf_group_id{-1};

  if (status != ERROR_VALUE) {
    switch (mesh->getType()) {
    case REC:
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
      // Temporary output while VTKHDF file format is not yet specified. For now, uses format specified on VTK discourse forum:
      // https://discourse.vtk.org/t/support-for-rectilineargrids-and-structuredgrids-in-vtkhdf-file-format/14283/6
      // An additional xmdf file is required for Paraview to interpret data.
      status = writeMeshRectilinearVTKHDF(mesh, file_id, vtkhdf_group_id, numberCells);
      break;
#endif
    case UNS:
    case AMR:
      status = writeMeshUnstructuredVTKHDF(mesh, cellsLvl, file_id, vtkhdf_group_id, numberCells);
      break;
    default:
      Errors::errorMessage("Output::writeResultsVTKHDF: type mesh unknown", __FILE__, __LINE__);
      break;
    }
  }

  //2) Write fluid phase data
  //-------------------------
  if (status != ERROR_VALUE) status = writePhysicalDataVTKHDF(file_id, mesh, cellsLvl, vtkhdf_group_id, numberCells);

  //3) Finalize file
  //----------------
  // Close all opened HDF5 objects and hdf5 file
  if (vtkhdf_group_id >= 0) H5Gclose(vtkhdf_group_id);

  herr_t status_close = H5Fclose(file_id);
  if (status_close < 0) {
    // We try to close the file but it fails... nothing more to do automatically
    throw ErrorECOGEN("Error when closing VTKHDF file.", __FILE__, __LINE__);
  }
  if (status == ERROR_VALUE) {
    // From here we can throw errors because the file has been closed cleanly
    throw ErrorECOGEN("Error when saving output results at VTKHDF file format.", __FILE__, __LINE__);
  }
}

//***********************************************************************

void OutputVTKHDF::writeTimeseriesVTKHDF(int numFile)
{
  try {
    std::ofstream fileStream;

    fileStream.open(m_filenameVTKHDFTimeseries + ".vtkhdf.series", std::ios::trunc);
    if (!fileStream) {
      throw ErrorECOGEN("Unable to open file " + m_filenameVTKHDFTimeseries, __FILE__, __LINE__);
    }
    fileStream << "{\n \"file-series-version\" : \"1.0\",\n  \"files\" : [" << std::endl;
    for (int time = 0; time <= numFile; time++) {
      std::string file = "datasets/" + createFilenameVTKHDF(m_fileNameResults.c_str(), time);
      fileStream << "    { \"name\" : \"" << file.c_str() << "\", \"time\" : " << time << " }," << std::endl;
    }
    fileStream << "  ]\n}" << std::endl;
    fileStream.close();
  }
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************

void OutputVTKHDF::writeDataXDMF(std::ofstream& fileStream,
                                 const std::string& dataFile,
                                 const std::string& dataName,
                                 const std::string& dataType,
                                 int numberCellsX,
                                 int numberCellsY,
                                 int numberCellsZ,
                                 int dimData)
{
  fileStream << "\t<Attribute Center=\"Cell\" Name=\"" << dataName << "\" Type=\"" << dataType << "\">\n";
  fileStream << "\t\t<DataItem DataType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << numberCellsZ << " " << numberCellsY << " "
             << numberCellsX << " " << dimData << "\">\n";
  fileStream << "\t\t\t" << dataFile << ":/VTKHDF/CellData/" << dataName << "\n";
  fileStream << "\t\t</DataItem>\n\t</Attribute>\n";
}

//***********************************************************************

void OutputVTKHDF::writeXDMF(Mesh* mesh, int numFile)
{
  const std::string closureString = "</Grid>\n</Domain>\n</Xdmf>\n";

  assert(DIM_OUTPUT == 3 && "Output not implemented in 1D and 2D");

  try {
    int numberCellsGlobalX{0};
    int numberCellsGlobalY{0};
    int numberCellsGlobalZ{0};
    if (mesh != nullptr) {
      numberCellsGlobalX = mesh->getNumberCellsXGlobal();
      numberCellsGlobalY = mesh->getNumberCellsYGlobal();
      numberCellsGlobalZ = mesh->getNumberCellsZGlobal();
    }

    auto numberPointsGlobalX{numberCellsGlobalX + 1};
    auto numberPointsGlobalY{numberCellsGlobalY + 1};
    auto numberPointsGlobalZ{numberCellsGlobalZ + 1};

    std::ofstream fileStream;
    std::string filename = m_filenameVTKHDFTimeseries + ".xdmf";

    if (numFile == -1) {
      // Open file in 'trunc' mode to delete eventually pre-existing file
      fileStream.open(filename.c_str(), std::ios::trunc);
      if (!fileStream) {
        throw ErrorECOGEN("Unable to open file " + filename, __FILE__, __LINE__);
      }

      // Write file header
      fileStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<Xdmf Version=\"3.0\">\n<Domain>\n";
      fileStream << "<Grid CollectionType=\"Temporal\" GridType=\"Collection\" Name=\"" << m_simulationName << "\">\n";
    }
    else {
      // Open file and seek position where to append time step data
      // 'in' mode is required to preserve data content, 'out' mode is required to seek and write in file
      fileStream.open(filename.c_str(), std::ios::in | std::ios::out);
      if (!fileStream) {
        throw ErrorECOGEN("Unable to open file " + filename, __FILE__, __LINE__);
      }
      // Set cursor to end of file
      fileStream.seekp(0, std::ios::end);
      // Find cursor position value
      long position = fileStream.tellp();
      // Shift cursor position from length of closure string to find position at which we want to write data
      fileStream.seekp(position - closureString.length());

      // 1. Save time step
      fileStream << "\t<Time Value=\"" << numFile << "\"/>\n";

      // 2. Pointers to mesh data
      fileStream << "<Grid Name=\"mesh_TIME" << numFile << "\" GridType=\"Uniform\">\n";
      // name of VTKHDF file of current iteration (contains physical data)
      std::string datafile = "datasets/" + createFilenameVTKHDF(m_fileNameResults.c_str(), numFile);
      // name of VTKHDF file of iteration 0 (the only one to contains the mesh definition)
      std::string datafile0 = "datasets/" + createFilenameVTKHDF(m_fileNameResults.c_str(), 0);

      fileStream << "\t<Topology Type=\"3DRectMesh\" Dimensions=\"" << numberPointsGlobalX << " " << numberPointsGlobalY << " " << numberPointsGlobalZ
                 << "\"/>\n";
      fileStream << "\t<Geometry GeometryType=\"VXVYVZ\">\n";
      // Xdmf format allow 4 or 8 precision. Changing this value and visualizing data in Paraview shows that precision '4' gives ~5 decimal digits and precision 8, 9 decimal digits.
      // As ECOGEN allows for 8 or 16 decimal digits only, we choose to use only precision 8 of xdmf.
      fileStream << "\t\t<DataItem DataType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << numberPointsGlobalX << "\">\n";
      fileStream << "\t\t\t" << datafile0 << ":/VTKHDF/XCoordinates\n";
      fileStream << "\t\t</DataItem>\n";
      fileStream << "\t\t<DataItem DataType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << numberPointsGlobalY << "\">\n";
      fileStream << "\t\t\t" << datafile0 << ":/VTKHDF/YCoordinates\n";
      fileStream << "\t\t</DataItem>\n";
      fileStream << "\t\t<DataItem DataType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << numberPointsGlobalZ << "\">\n";
      fileStream << "\t\t\t" << datafile0 << ":/VTKHDF/ZCoordinates\n";
      fileStream << "\t\t</DataItem>\n\t</Geometry>\n";

      // 3. Pointers to data fields
      // 3.1 Phase variables
      for (int phase = 0; phase < m_run->getNumberPhases(); phase++) //For complete output
      //for (int phase = 0; phase < 1; phase++) //For reduced output
      {
        std::string eosName; //empty string
        if (m_cellRef.getPhase(phase)->getEos() != nullptr) {
          eosName = "_" + m_cellRef.getPhase(phase)->getEos()->getName();
          eosName.erase(eosName.end() - 4, eosName.end());
        }

        // Pointers to scalar variables
        for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberScalars(); var++) {
          std::string dataName{outputPhaseName(phase, m_cellRef.getPhase(phase)->returnNameScalar(var), eosName)};
          writeDataXDMF(fileStream, datafile, dataName, "Scalar", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, 1);
        }
        // Pointers to vectorial variables
        for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++) {
          std::string dataName{outputPhaseName(phase, m_cellRef.getPhase(phase)->returnNameVector(var), eosName)};
          writeDataXDMF(fileStream, datafile, dataName, "Vector", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, DIM_OUTPUT);
        }
      }
      // 3.2 Mixture data
      if (m_run->m_numberPhases > 1) {
        for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++) {
          std::string dataName{m_cellRef.getMixture()->returnNameScalar(var)};
          writeDataXDMF(fileStream, datafile, dataName, "Scalar", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, 1);
        }
        for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++) {
          std::string dataName{m_cellRef.getMixture()->returnNameVector(var)};
          writeDataXDMF(fileStream, datafile, dataName, "Vector", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, DIM_OUTPUT);
        }
      }

      //3.3 Transports and others...
      for (int var = 1; var <= m_run->m_numberTransports; var++) {
        std::string dataName{"T" + std::to_string(var)};
        writeDataXDMF(fileStream, datafile, dataName, "Scalar", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, 1);
      }

      //3.4 xi indicator
      if (mesh->getType() == AMR) { //For complete output
        writeDataXDMF(fileStream, datafile, "Xi", "Scalar", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, 1);
      }

      //3.5 Gradient rho
      writeDataXDMF(fileStream, datafile, "gradRho", "Scalar", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, 1);

      //3.6 Absolute velocity  for Moving Reference Frame computations
      if (m_run->m_MRF != -1) {
        writeDataXDMF(fileStream, datafile, "absoluteVelocityMRF", "Vector", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, DIM_OUTPUT);
      }

      //3.7 Cells' reference length
      if (m_run->m_extractRefLength) {
        writeDataXDMF(fileStream, datafile, "Reference_length", "Scalar", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, 1);
      }

      //3.8 CPU rank
      writeDataXDMF(fileStream, datafile, "CPUrank", "Scalar", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, 1);

      //3.9 Morton index
      //if (mesh->getType() == AMR) { //For complete output
      //  writeDataXDMF(fileStream,datafile,"MortonIndex","Scalar", numberCellsGlobal,1);
      //}

      //3.10 Saturation pressure (specific recording)
      if (m_run->m_recordPsat) {
        writeDataXDMF(fileStream, datafile, "Psat", "Scalar", numberCellsGlobalX, numberCellsGlobalY, numberCellsGlobalZ, 1);
      }

      fileStream << "</Grid>\n";
    } // End if (numFile==-1)

    // End of file
    fileStream << closureString;
    fileStream.close();
  }
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************

template <size_t N, typename T>
int OutputVTKHDF::writePhysicalDataVTKHDF(const Mesh* mesh,
                                          std::vector<Cell*>* cellsLvl,
                                          const hid_t vtkhdf_group_id,
                                          const size_t numberCells,
                                          const hsize_t (&dimsf_scalar)[N],
                                          const hsize_t (&count_scalar)[N],
                                          const hsize_t (&start_scalar)[N],
                                          const hsize_t (&chunk_scalar)[N]) const
{
  int status{SUCCESS_VALUE};
  std::vector<T> dataset_scalar, dataset_vector;

  // CellData Group
  hid_t celldata_group_id{H5Gcreate2(vtkhdf_group_id, "CellData", H5P_DEFAULT, m_group_creation_properties, H5P_DEFAULT)};
  if (celldata_group_id < 0) {
    Errors::errorMessage("Unable to create 'CellData' group", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Create arrays storing dimensions and offset for vectorial data: dimensions will be the
  // same than for scalar parsing plus one additional range of dimension DIM_OUTPUT
  hsize_t dimsf_vector[N + 1];
  hsize_t count_vector[N + 1];
  hsize_t start_vector[N + 1];
  hsize_t chunk_vector[N + 1];

  for (size_t i = 0; i < N; ++i) {
    dimsf_vector[i] = dimsf_scalar[i];
    count_vector[i] = count_scalar[i];
    start_vector[i] = start_scalar[i];
    chunk_vector[i] = chunk_scalar[i];
  }
  dimsf_vector[N] = DIM_OUTPUT;
  count_vector[N] = DIM_OUTPUT;
  start_vector[N] = 0;
  chunk_vector[N] = DIM_OUTPUT;

  // Write integers values first in case it improve efficiency (it can minimize time gap between processors)

  //1) CPU rank
  //-----------
  if (status != ERROR_VALUE) {
    std::vector<int> dataset_i(numberCells, rankCpu); // Allocate dataset_i vector at size 'numberCells' and initialize it at value 'rankCpu'

    status = writeArrayVTKHDF(celldata_group_id,
                              "CPUrank",
                              H5T_NATIVE_INT,
                              H5T_NATIVE_INT,
                              m_parallel_x_properties,
                              dataset_i,
                              dimsf_scalar,
                              count_scalar,
                              start_scalar,
                              chunk_scalar);
  }

  //2) Write phases variables
  //---------------------------------
  for (int phase = 0; phase < m_run->getNumberPhases(); phase++) //For complete output
  //for (int phase = 0; phase < 1; phase++) //For reduced output
  {
    std::string eosName; //empty string
    if (m_cellRef.getPhase(phase)->getEos() != nullptr) {
      eosName = "_" + m_cellRef.getPhase(phase)->getEos()->getName();
      eosName.erase(eosName.end() - 4, eosName.end());
    }

    //Write scalar variables
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberScalars(); var++) {

      if (status != ERROR_VALUE) {
        mesh->getData(cellsLvl, dataset_scalar, var, phase);
        assert(numberCells == dataset_scalar.size() && "Mismatch number of Cells and number of CellData");

        std::string dataName{outputPhaseName(phase, m_cellRef.getPhase(phase)->returnNameScalar(var), eosName)};
        status = writeArrayVTKHDF(celldata_group_id,
                                  dataName,
                                  m_h5_mem_floatType,
                                  m_h5_file_floatType,
                                  m_parallel_x_properties,
                                  dataset_scalar,
                                  dimsf_scalar,
                                  count_scalar,
                                  start_scalar,
                                  chunk_scalar);
      }
    }

    //Write vectorial variables
    for (int var = 1; var <= m_cellRef.getPhase(phase)->getNumberVectors(); var++) {
      if (status != ERROR_VALUE) {
        mesh->getData(cellsLvl, dataset_vector, -var, phase);
        assert(numberCells == dataset_vector.size() / DIM_OUTPUT && "Mismatch number of Cells and number of CellData");

        std::string dataName{outputPhaseName(phase, m_cellRef.getPhase(phase)->returnNameVector(var), eosName)};
        status = writeArrayVTKHDF(celldata_group_id,
                                  dataName,
                                  m_h5_mem_floatType,
                                  m_h5_file_floatType,
                                  m_parallel_x_properties,
                                  dataset_vector,
                                  dimsf_vector,
                                  count_vector,
                                  start_vector,
                                  chunk_vector);
      }
    }
  } //End phase

  //3) Write mixture data
  //---------------------
  if (m_run->m_numberPhases > 1) {
    int mixture{-1};
    //Write scalars mixture variables
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberScalars(); var++) {
      if (status != ERROR_VALUE) {
        mesh->getData(cellsLvl, dataset_scalar, var, mixture);
        assert(numberCells == dataset_scalar.size() && "Mismatch number of Cells and number of CellData");
        std::string dataName{m_cellRef.getMixture()->returnNameScalar(var)};
        status = writeArrayVTKHDF(celldata_group_id,
                                  dataName,
                                  m_h5_mem_floatType,
                                  m_h5_file_floatType,
                                  m_parallel_x_properties,
                                  dataset_scalar,
                                  dimsf_scalar,
                                  count_scalar,
                                  start_scalar,
                                  chunk_scalar);
      }
    }
    //Write vectorial mixture variables
    for (int var = 1; var <= m_cellRef.getMixture()->getNumberVectors(); var++) {
      if (status != ERROR_VALUE) {
        mesh->getData(cellsLvl, dataset_vector, -var, mixture);
        assert(numberCells == dataset_vector.size() / DIM_OUTPUT && "Mismatch number of Cells and number of CellData");
        std::string dataName{m_cellRef.getMixture()->returnNameVector(var)};
        status = writeArrayVTKHDF(celldata_group_id,
                                  dataName,
                                  m_h5_mem_floatType,
                                  m_h5_file_floatType,
                                  m_parallel_x_properties,
                                  dataset_vector,
                                  dimsf_vector,
                                  count_vector,
                                  start_vector,
                                  chunk_vector);
      }
    }
  } //End mixture

  //4) Write of transports and others...
  //------------------------------------
  int transport = -2; //For complete output
  for (int var = 1; var <= m_run->m_numberTransports; var++) {
    if (status != ERROR_VALUE) {
      mesh->getData(cellsLvl, dataset_scalar, var, transport);
      assert(numberCells == dataset_scalar.size() && "Mismatch number of Cells and number of CellData");
      std::string dataName{"T" + std::to_string(var)};
      status = writeArrayVTKHDF(celldata_group_id,
                                dataName,
                                m_h5_mem_floatType,
                                m_h5_file_floatType,
                                m_parallel_x_properties,
                                dataset_scalar,
                                dimsf_scalar,
                                count_scalar,
                                start_scalar,
                                chunk_scalar);
    }
  }

  //5) Write xi indicator
  //----------------------
  if (mesh->getType() == AMR) { //For complete output
    if (status != ERROR_VALUE) {
      int xi{-3};
      mesh->getData(cellsLvl, dataset_scalar, 1, xi);
      assert(numberCells == dataset_scalar.size() && "Mismatch number of Cells and number of CellData");
      status = writeArrayVTKHDF(celldata_group_id,
                                "Xi",
                                m_h5_mem_floatType,
                                m_h5_file_floatType,
                                m_parallel_x_properties,
                                dataset_scalar,
                                dimsf_scalar,
                                count_scalar,
                                start_scalar,
                                chunk_scalar);
    }
  }

  //6) Write gradient rho
  //---------------------
  if (status != ERROR_VALUE) {
    int gradRho_{-4};
    mesh->getData(cellsLvl, dataset_scalar, 1, gradRho_);
    assert(numberCells == dataset_scalar.size() && "Mismatch number of Cells and number of CellData");
    status = writeArrayVTKHDF(celldata_group_id,
                              "gradRho",
                              m_h5_mem_floatType,
                              m_h5_file_floatType,
                              m_parallel_x_properties,
                              dataset_scalar,
                              dimsf_scalar,
                              count_scalar,
                              start_scalar,
                              chunk_scalar);
  }

  //7) Absolute velocity printing for Moving Reference Frame computations
  //---------------------------------------------------------------------
  if (m_run->m_MRF != -1) {
    if (status != ERROR_VALUE) {
      mesh->extractAbsVelocityMRF(cellsLvl, dataset_vector, m_run->m_sources[m_run->m_MRF]);
      assert(numberCells == dataset_vector.size() / DIM_OUTPUT && "Mismatch number of Cells and number of CellData");
      status = writeArrayVTKHDF(celldata_group_id,
                                "absoluteVelocityMRF",
                                m_h5_mem_floatType,
                                m_h5_file_floatType,
                                m_parallel_x_properties,
                                dataset_vector,
                                dimsf_vector,
                                count_vector,
                                start_vector,
                                chunk_vector);
    }
  }

  //8) Cells' reference length
  //--------------------------
  if (m_run->m_extractRefLength) {
    if (status != ERROR_VALUE) {
      mesh->extractReferenceLength(cellsLvl, dataset_scalar);
      assert(numberCells == dataset_scalar.size() && "Mismatch number of Cells and number of CellData");
      status = writeArrayVTKHDF(celldata_group_id,
                                "Reference_length",
                                m_h5_mem_floatType,
                                m_h5_file_floatType,
                                m_parallel_x_properties,
                                dataset_scalar,
                                dimsf_scalar,
                                count_scalar,
                                start_scalar,
                                chunk_scalar);
    }
  }

  //9) Morton index
  //---------------
  //if (mesh->getType() == AMR) { //For complete output
  //  if ( status != ERROR_VALUE ) {
  //    int MortonIndex = -6;
  //    mesh->getData(cellsLvl, dataset_scalar, 1, MortonIndex);
  //    assert ( numberCells == dataset_scalar.size() && "Mismatch number of Cells and number of CellData") ;
  //
  //    std::vector<int64_t> dataset_i;
  //    dataset_i.reserve(numberCells);
  //    for (double d : dataset_scalar) {
  //      // Warning: as this part of code is commented, I have not implemented the check of the validity of the cast
  //      dataset_i.push_back(static_cast<int64_t>(std::round(d)));
  //    }
  //    status = writeArrayVTKHDF(celldata_group_id,"MortonIndex",H5T_NATIVE_INT64,H5T_STD_I64LE,m_parallel_x_properties,dataset_i,dimsf_scalar,count_scalar,start_scalar);
  //  }
  //}

  //10) Writing saturation pressure (specific recording)
  //---------------------------------------------------
  if (m_run->m_recordPsat) {
    if (status != ERROR_VALUE) {
      mesh->getData(cellsLvl, dataset_scalar, 1, -7);
      assert(numberCells == dataset_scalar.size() && "Mismatch number of Cells and number of CellData");
      status = writeArrayVTKHDF(celldata_group_id,
                                "Psat",
                                m_h5_mem_floatType,
                                m_h5_file_floatType,
                                m_parallel_x_properties,
                                dataset_scalar,
                                dimsf_scalar,
                                count_scalar,
                                start_scalar,
                                chunk_scalar);
    }
  }

  // Close all opened HDF5 objects
  if (celldata_group_id >= 0) H5Gclose(celldata_group_id);

  return status;
}

//***********************************************************************

int OutputVTKHDF::writePhysicalDataVTKHDF(
  hid_t& file, const Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t& vtkhdf_group_id, const size_t numberCells) const
{

  int status{SUCCESS_VALUE};

  switch (mesh->getType()) {
  case REC:
#ifndef ECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS
  {
    // Each rank writes in a location that allow to have a serial output.
    // It is necessary to save array in 'natural' ranks and dimension because hdf stores data in row-major mode while ECOGEN flattened data using column-major order.
    auto numberCellsXGlobal{mesh->getNumberCellsXGlobal()};
    auto numberCellsYGlobal{mesh->getNumberCellsYGlobal()};
    auto numberCellsZGlobal{mesh->getNumberCellsZGlobal()};

    auto numberCellsX{mesh->getNumberCellsX()};
    auto numberCellsY{mesh->getNumberCellsY()};
    auto numberCellsZ{mesh->getNumberCellsZ()};

    auto offsetX{mesh->getOffsetX()};
    auto offsetY{mesh->getOffsetY()};
    auto offsetZ{mesh->getOffsetZ()};

    hsize_t dimsf_scalar[DIM_OUTPUT] = {(hsize_t)numberCellsZGlobal, (hsize_t)numberCellsYGlobal, (hsize_t)numberCellsXGlobal};
    hsize_t count_scalar[DIM_OUTPUT] = {(hsize_t)numberCellsZ, (hsize_t)numberCellsY, (hsize_t)numberCellsX};
    hsize_t start_scalar[DIM_OUTPUT] = {(hsize_t)offsetZ, (hsize_t)offsetY, (hsize_t)offsetX};

    // Compute chunking size: we have to ensure that chunk size is the same on all ranks (so we cannot use local_dims directly).
    // From our tests, when all partitions doesn't have the same sizes, it is more efficient to use the largest partition size as chunk size than the smallest one.
    int maxNumberCells[3] = {numberCellsX, numberCellsY, numberCellsZ};
    MPI_Allreduce(MPI_IN_PLACE, &maxNumberCells, 3, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    hsize_t chunk_scalar[DIM_OUTPUT] = {(hsize_t)maxNumberCells[2], (hsize_t)maxNumberCells[1], (hsize_t)maxNumberCells[0]};

    if (m_numFile == 0) {
      // At time step 0, file has been initially opened (and closed) in serial mode.
      // Re-open file in parallel mode.
      constexpr bool parallel_access = true;
      constexpr bool create          = false;
      status                         = openFileVTKHDF(file, parallel_access, create);

      if (status != ERROR_VALUE) {
        // Re-open 'VTKHDF' group
        vtkhdf_group_id = H5Gopen2(file, "VTKHDF", H5P_DEFAULT);
      }

      if (vtkhdf_group_id < 0) {
        Errors::errorMessage("Unable to create 'VTKHDF' group", __FILE__, __LINE__);
        status = ERROR_VALUE;
      }
    }

    if (status != ERROR_VALUE) {
      if (m_h5_file_floatType == H5T_IEEE_F32LE) {
        status = writePhysicalDataVTKHDF<DIM_OUTPUT, float>(
          mesh, cellsLvl, vtkhdf_group_id, numberCells, dimsf_scalar, count_scalar, start_scalar, chunk_scalar);
      }
      else {
        status = writePhysicalDataVTKHDF<DIM_OUTPUT, double>(
          mesh, cellsLvl, vtkhdf_group_id, numberCells, dimsf_scalar, count_scalar, start_scalar, chunk_scalar);
      }
    }
    break;
  }
#else
    (void)file; // Trick to avoid warning about unused parameter
#endif
  case UNS:
  case AMR:
    hsize_t dimsf_scalar[1];
    hsize_t count_scalar[1];
    hsize_t start_scalar[1];
    hsize_t chunk_scalar[1];

    // Each rank writes in an area of size numberCells  ( x DIM_OUTPUT for vectorial data)
    unstructuredFirstRangeLayout(numberCells, dimsf_scalar, count_scalar, start_scalar);

    // disable chunking for unstructured data
    chunk_scalar[0] = 0;

    if (m_h5_file_floatType == H5T_IEEE_F32LE) {
      status = writePhysicalDataVTKHDF<1, float>(mesh, cellsLvl, vtkhdf_group_id, numberCells, dimsf_scalar, count_scalar, start_scalar, chunk_scalar);
    }
    else {
      status =
        writePhysicalDataVTKHDF<1, double>(mesh, cellsLvl, vtkhdf_group_id, numberCells, dimsf_scalar, count_scalar, start_scalar, chunk_scalar);
    }

    break;
  default:
    Errors::errorMessage("Output::writePhysicalDataVTKHDF: type mesh unknown", __FILE__, __LINE__);
    break;
  }
  return status;
}

//***********************************************************************
//
int OutputVTKHDF::writeArrayAttribute(const hid_t& group_id, const std::string& attribute_name, const hsize_t& nitem, const int* const array) const
{
  int status{SUCCESS_VALUE};

  const hsize_t dim[1] = {nitem}; // Dimension of array to write

  hid_t dataspace{-1}, attr{-1};
  dataspace = H5Screate_simple(1, dim, NULL); // dataspace for a 1-dimensional array of dimension dim
  if (dataspace < 0) {
    Errors::errorMessage("Unable to create datapace for '" + attribute_name + "' attribute", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // create the attribute variable
  if (status != ERROR_VALUE) {
    attr = H5Acreate2(group_id, attribute_name.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
      Errors::errorMessage("Unable to create hdf5 '" + attribute_name + "' attribute variable", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  // Write attribute value
  if (status != ERROR_VALUE && H5Awrite(attr, H5T_NATIVE_INT, array) < 0) {
    Errors::errorMessage("Error when writing '" + attribute_name + "' attribute in VTKHDF file", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Close opened HDF5 objects
  if (dataspace >= 0) H5Sclose(dataspace);
  if (attr >= 0) H5Aclose(attr);

  return status;
}
//***********************************************************************

int OutputVTKHDF::writeHeaderVTKHDF(const hid_t& group_id, const std::string& mesh_type) const
{
  hid_t str_type{-1};
  hid_t scalar_dataspace{-1};
  hid_t char_encoding{-1};
  hid_t type_attr{-1};

  int status{SUCCESS_VALUE};

  // 1. Write "Type" attribute
  size_t len{mesh_type.length()};

  str_type         = H5Tcreate(H5T_STRING, len);      // create a new string datatype
  scalar_dataspace = H5Screate(H5S_SCALAR);           // define a scalar dataspace
  char_encoding    = H5Pcreate(H5P_ATTRIBUTE_CREATE); // variable to store char encoding property list
  if (str_type < 0 || scalar_dataspace < 0 || char_encoding < 0) {
    Errors::errorMessage("Unable to create hdf5 variables", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Set UTF-8 encoding
  if (status != ERROR_VALUE && H5Pset_char_encoding(char_encoding, H5T_CSET_UTF8) < 0) {
    Errors::errorMessage("Unable to set UTF-8 encoding for hdf5 'Type' attribute", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // create the "Type" attribute variable
  if (status != ERROR_VALUE) {
    type_attr = H5Acreate2(group_id, "Type", str_type, scalar_dataspace, char_encoding, H5P_DEFAULT);
    if (type_attr < 0) {
      Errors::errorMessage("Unable to create hdf5 'Type' attribute variable", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  // Write "Type" attribute value
  if (status != ERROR_VALUE && H5Awrite(type_attr, str_type, mesh_type.c_str()) < 0) {
    Errors::errorMessage("Error when writing 'Type' attribute in VTKHDF file", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // 2. Write "Version" attribute
  if (status != ERROR_VALUE) {
    status = writeArrayAttribute(group_id, "Version", 2, VTKHDF_VERSION);
    if (status == ERROR_VALUE) {
      Errors::errorMessage("Unable to create datapace for 'Version' attribute", __FILE__, __LINE__);
    }
  }

  // Close opened HDF5 objects
  if (str_type >= 0) H5Tclose(str_type);
  if (scalar_dataspace >= 0) H5Sclose(scalar_dataspace);
  if (char_encoding >= 0) H5Pclose(char_encoding);
  if (type_attr >= 0) H5Aclose(type_attr);

  return status;
}

//***********************************************************************

int OutputVTKHDF::writeScalarVTKHDF(const hid_t& group_id, const std::string& dataset_name, const size_t& scalar_val) const
{
  hid_t dataspace{-1};
  hid_t hdf_dataset{-1};
  hid_t memspace{-1};

  int status{SUCCESS_VALUE};

  hsize_t dimsf[1], count[1], start[1];

  // Create a dataspace to store a 1-dimensional array of size Ncpu
  dimsf[0]  = Ncpu; // Global dimension of the array to write across MPI ranks
  dataspace = H5Screate_simple(1, dimsf, NULL);
  if (dataspace < 0) {
    Errors::errorMessage("Unable to create hdf5 dataspace", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  if (status != ERROR_VALUE) {
    // Create the hdf dataset
    // H5T_STD_*BE and  H5T_STD_*LE dataset types are supported by VTKHDF file format so we can use hdf5 native datatypes.
    hdf_dataset = H5Dcreate(group_id, dataset_name.c_str(), ECOGEN_H5_SIZE_T, dataspace, H5P_DEFAULT, m_dataset_creation_properties, H5P_DEFAULT);
    if (hdf_dataset < 0) {
      Errors::errorMessage("Unable to create hdf5 dataset to save " + dataset_name, __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  // Select the wanted hyperslab region in which current rank will write
  count[0] = 1;       // number of data to write for the current rank
  start[0] = rankCpu; // as we write 1 data per rank, each MPI rank will write at the position given by its index
  if (status != ERROR_VALUE && H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
    Errors::errorMessage("Unable to create hdf5 hyperslab to save " + dataset_name, __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Create a memory space for the current rank
  memspace = H5Screate_simple(1, count, NULL);
  if (memspace < 0) {
    Errors::errorMessage("Unable to create hdf5 memory space", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Write value
  if (status != ERROR_VALUE && H5Dwrite(hdf_dataset, ECOGEN_H5_SIZE_T, memspace, dataspace, m_parallel_x_properties, &scalar_val) < 0) {
    Errors::errorMessage("Unable to create hdf5 memory space", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Close hdf5 objects
  if (dataspace >= 0) H5Sclose(dataspace);
  if (memspace >= 0) H5Sclose(memspace);
  if (hdf_dataset >= 0) H5Dclose(hdf_dataset);

  return status;
}

//***********************************************************************

template <size_t N, typename T>
int OutputVTKHDF::writeArrayVTKHDF(const hid_t& group_id,
                                   const std::string& dataset_name,
                                   hid_t hdf_mem_dtype,
                                   hid_t hdf_file_dtype,
                                   hid_t parallel_x_properties,
                                   const std::vector<T>& dataset,
                                   const hsize_t (&global_dims)[N],
                                   const hsize_t (&local_dims)[N],
                                   const hsize_t (&offsets)[N],
                                   const hsize_t (&chunk_dims)[N]) const
{
  hid_t dataspace{-1};
  hid_t memspace{-1};
  hid_t hdf_dataset{-1};

  int status{SUCCESS_VALUE};

  dataspace = H5Screate_simple(N, global_dims, NULL);
  if (dataspace < 0) {
    Errors::errorMessage("Unable to create hdf5 dataspace", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  if (chunk_dims[0] != 0) { // Use chunks only if chunk_dims is different from 0
    if (status != ERROR_VALUE && H5Pset_chunk(m_dataset_creation_properties, N, chunk_dims) < 0) {
      Errors::errorMessage("Unable to set chunking", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  // Create the hdf dataset
  if (status != ERROR_VALUE) {
    // Both big endian and little endian dataset types seems to be supported by VTKHDF format so we can use hdf5 native
    // type alias required by the 'H5Dwrite' function and use only on hdf5 type (stored in the 'hdf_type' variable)
    hdf_dataset = H5Dcreate(group_id, dataset_name.c_str(), hdf_file_dtype, dataspace, H5P_DEFAULT, m_dataset_creation_properties, H5P_DEFAULT);
    if (hdf_dataset < 0) {
      Errors::errorMessage("Unable to create hdf5 dataset to save " + dataset_name, __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  if (status != ERROR_VALUE && H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, local_dims, NULL) < 0) {
    Errors::errorMessage("Unable to create hdf5 hyperslab to save " + dataset_name, __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Create a memory space for the current rank
  memspace = H5Screate_simple(N, local_dims, NULL);
  if (memspace < 0) {
    Errors::errorMessage("Unable to create hdf5 memory space", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Write value
  if (status != ERROR_VALUE && H5Dwrite(hdf_dataset, hdf_mem_dtype, memspace, dataspace, parallel_x_properties, &dataset[0]) < 0) {
    Errors::errorMessage("Unable to create hdf5 memory space", __FILE__, __LINE__);
    status = ERROR_VALUE;
  }

  // Close hdf5 objects
  if (dataspace >= 0) H5Sclose(dataspace);
  if (memspace >= 0) H5Sclose(memspace);
  if (hdf_dataset >= 0) H5Dclose(hdf_dataset);

  return status;
}

//***********************************************************************

void OutputVTKHDF::unstructuredFirstRangeLayout(const size_t& line_count, hsize_t* dimsf, hsize_t* count, hsize_t* start)
{
  size_t globalCount{0}, initialPosition{0};
  std::vector<size_t> countPerRank(Ncpu);
  // Determine array distribution accros MPI ranks (each rank contains 'line_count' lines of the array):
  // compute global line count and initial position for current process
  MPI_Allgather(&line_count, 1, ECOGEN_MPI_SIZE_T, countPerRank.data(), 1, ECOGEN_MPI_SIZE_T, MPI_COMM_WORLD);

  for (int i = 0; i < rankCpu; ++i) {
    initialPosition += countPerRank[i];
  }
  globalCount = initialPosition;
  for (int i = rankCpu; i < Ncpu; i++) {
    globalCount += countPerRank[i];
  }

  // Create a dataspace to store array of size globalCount x column_count
  dimsf[0] = globalCount; // Global size across MPI ranks of the first dimension of the array to write

  // Select the wanted hyperslab region in which current rank will write
  count[0] = line_count; // number of lines of data to write for the current rank

  start[0] = initialPosition; // Line index from which current rank will start to write
}

//***********************************************************************

template <typename T>
int OutputVTKHDF::writeUnstructuredArrayVTKHDF(const hid_t& group_id,
                                               const std::string& dataset_name,
                                               const hid_t hdf_mem_dtype,
                                               const hid_t hdf_file_dtype,
                                               const std::vector<T>& dataset,
                                               const size_t& line_count,
                                               const size_t column_count) const
{
  hsize_t dimsf[2], count[2], start[2], chunk[2];

  // Compute dimensions and offset for first range (line count) for an array for which each MPI
  // process writes its 'line_count' lines in order of cpu rank (rank 0 writes its 'line_count'
  // lines first, then rank 1, etc.)
  unstructuredFirstRangeLayout(line_count, dimsf, count, start);

  // Fill dimensions and offset for second range (column count): all columns are written by each MPI process
  dimsf[1] = column_count; // Size of the second dimension of the array to write
  count[1] = column_count; // number of columns to write (each rank writes always all the columns in our case)
  start[1] = 0;            // column index from which current rank will start to write (always 0 in our case)

  // Disable cunking for unstructured arrays
  chunk[0] = chunk[1] = 0;
  return writeArrayVTKHDF(group_id, dataset_name, hdf_mem_dtype, hdf_file_dtype, m_parallel_x_properties, dataset, dimsf, count, start, chunk);
}

//***********************************************************************

int OutputVTKHDF::writeMeshRectilinearVTKHDF(Mesh* mesh, hid_t& file, hid_t& vtkhdf_group_id, size_t& numberCells)

{
  int status{SUCCESS_VALUE};

  // We will save the grid only for time step 0. In this case, a serial access is needed.
  // For other time steps, the file can directly be opened in parallel
  const bool parallel_access = (m_numFile == 0) ? false : true;
  constexpr bool create      = true;

  //1) Create and open file
  //-----------------------
  if ((m_numFile == 0 && rankCpu == 0) || m_numFile != 0) { // Non-0 CPU should not open the file at time 0 because file access is serial.
    status = openFileVTKHDF(file, parallel_access, create);

    //2) Root group
    //-------------
    // VTKHDF Group
    if (status != ERROR_VALUE) {
      vtkhdf_group_id = H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, m_group_creation_properties, H5P_DEFAULT);
      if (vtkhdf_group_id < 0) {
        Errors::errorMessage("Unable to create 'VTKHDF' group", __FILE__, __LINE__);
        status = ERROR_VALUE;
      }
    }

    //3) Header
    //---------
    if (status != ERROR_VALUE) {
      status = writeHeaderVTKHDF(vtkhdf_group_id, "RectilinearGrid"); // This field is useful at all time steps to identify the type of VTKHDF file
                                                                      // (Unstructured versus Rectilinear).
                                                                      // It is used for example by 'compare_vtk.py' script.
    }
  }

  //4) Grid is obtained by tensor product of vectors of axis: get coordinates of these vectors
  //------------------------------------------------------------------------------------------
  // Grid is saved at time 0 and by rank 0 only, with serial access.
  //
  //  Write "WholeExtent" specific attribute
  if (m_numFile == 0) {
    if (rankCpu == 0) {

      std::array<int, 6> extent = mesh->getArrayGlobalExtent();
      if (status != ERROR_VALUE) {
        status = writeArrayAttribute(vtkhdf_group_id, "WholeExtent", 6, extent.data());
        if (status == ERROR_VALUE) {
          Errors::errorMessage("Unable to create datapace for 'WholeExtent' attribute", __FILE__, __LINE__);
        }
      }

      auto numberCellsX{mesh->getNumberCellsXGlobal()};
      auto numberCellsY{mesh->getNumberCellsYGlobal()};
      auto numberCellsZ{mesh->getNumberCellsZGlobal()};
      bool global = true;

      std::vector<double> dataset;
      hsize_t count[1], start[1], chunk[1];

      hsize_t dimsf[1] = {static_cast<hsize_t>(numberCellsX + 1)};
      mesh->getCoord(dataset, X, global);
      count[0] = static_cast<hsize_t>(numberCellsX + 1);
      start[0] = 0;
      chunk[0] = 0; // chunking is disabled

      hid_t parallel_x_properties = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(parallel_x_properties, H5FD_MPIO_INDEPENDENT); // As file is opened in serial mode, COLLECTIVE mode would create an HDF5 error
      if (status != ERROR_VALUE) {
        status = writeArrayVTKHDF(
          vtkhdf_group_id, "XCoordinates", H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, parallel_x_properties, dataset, dimsf, count, start, chunk);
      }

      dimsf[0] = static_cast<hsize_t>(numberCellsY + 1);
      dataset.clear();
      mesh->getCoord(dataset, Y, global);
      count[0] = static_cast<hsize_t>(numberCellsY + 1);
      start[0] = static_cast<hsize_t>(0);

      if (status != ERROR_VALUE) {
        status = writeArrayVTKHDF(
          vtkhdf_group_id, "YCoordinates", H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, parallel_x_properties, dataset, dimsf, count, start, chunk);
      }

      dimsf[0] = static_cast<hsize_t>(numberCellsZ + 1);
      dataset.clear();
      mesh->getCoord(dataset, Z, global);
      count[0] = static_cast<hsize_t>(numberCellsZ + 1);
      start[0] = static_cast<hsize_t>(0);

      if (status != ERROR_VALUE) {
        status = writeArrayVTKHDF(
          vtkhdf_group_id, "ZCoordinates", H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, parallel_x_properties, dataset, dimsf, count, start, chunk);
      }

      // Close hdf objects and file
      if (vtkhdf_group_id >= 0) H5Gclose(vtkhdf_group_id);
      herr_t status_close = H5Fclose(file);
      if (status_close < 0) {
        // We try to close the file but it fails... nothing more to do automatically
        Errors::errorMessage("Error when closing VTKHDF file.", __FILE__, __LINE__);
        status = ERROR_VALUE;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // Store local cell count in numberCells
  numberCells = mesh->getNumberCells();

  return status;
}

//***********************************************************************

template <typename T>
int OutputVTKHDF::writeConnectivity(
  Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t& vtkhdf_group_id, size_t& numberCells, const size_t& numberNodes) const
{
  int status{SUCCESS_VALUE};
  std::vector<T> dataset_i;

  if (numberNodes > std::numeric_limits<T>::max()) {
    Errors::errorMessage("Too many nodes to call current implementation of '" + (std::string) __func__ +
                         "' method.\n"
                         "Consider changing the 'intSize' input variable to 64 or reducing your mesh size.");
    status = ERROR_VALUE;
  }
  if (status != ERROR_VALUE) {
    numberCells = mesh->getConnectivity(dataset_i, cellsLvl);
  }

  // Step 1: write number of cells
  if (status != ERROR_VALUE) {
    status = writeScalarVTKHDF(vtkhdf_group_id, "NumberOfCells", numberCells);
  }

  // Step 2: write number of connectivity indices
  size_t numberConnectivityIds{dataset_i.size()};
  if (status != ERROR_VALUE) {
    status = writeScalarVTKHDF(vtkhdf_group_id, "NumberOfConnectivityIds", numberConnectivityIds);
  }

  // Step 3: write connectivity indices

  // Format 2.4 of VTKHDF files doesn't support unsigned integer data set types (H5T_STD_UI64LE or H5T_STD_UI32LE types)
  // so we have to check that our 'size_t' data ('uint32_t' or 'uint64_t' depending on the architecture) can be casted into
  // the chosen size of signed integer.

  // If the number of nodes holds in a T variable, then the values of the connectivity array holds too
  if (status != ERROR_VALUE) {
    status = writeUnstructuredArrayVTKHDF(vtkhdf_group_id, "Connectivity", m_h5_mem_intType, m_h5_file_intType, dataset_i, numberConnectivityIds);
  }

  // Step 4: write offset
  dataset_i.clear();
  dataset_i.reserve(numberCells + 1);
  dataset_i.push_back(0); // In VTKHDF format, Offsets array gives the index of the start of each cell in the 'Connectivity' array
  mesh->getOffsets(dataset_i, cellsLvl);
  assert(dataset_i.size() == numberCells + 1 && "Size of Offset array not equal to NumberOfCells+1 as expected");
  if (status != ERROR_VALUE) {
    status = writeUnstructuredArrayVTKHDF(vtkhdf_group_id, "Offsets", m_h5_mem_intType, m_h5_file_intType, dataset_i, numberCells + 1);
  }

  // Step 5: write Cells type
  std::vector<unsigned char> dataset_c;
  dataset_c.reserve(numberCells);
  mesh->getTypeCell(dataset_c, cellsLvl);
  assert(dataset_c.size() == numberCells && "Size of Types array not equal to NumberOfCells as expected");
  if (status != ERROR_VALUE) {
    status = writeUnstructuredArrayVTKHDF(vtkhdf_group_id, "Types", H5T_NATIVE_UCHAR, H5T_NATIVE_UCHAR, dataset_c, numberCells);
  }
  return status;
}

//***********************************************************************

int OutputVTKHDF::writeMeshUnstructuredVTKHDF(Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t& file, hid_t& vtkhdf_group_id, size_t& numberCells) const
{
  int status{SUCCESS_VALUE};

  //1) Create and open file
  //-----------------------
  constexpr bool parallel_access = true;
  constexpr bool create          = true;
  status                         = openFileVTKHDF(file, parallel_access, create);

  //2) Root group
  //-------------
  // VTKHDF Group
  if (status != ERROR_VALUE) {
    vtkhdf_group_id = H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, m_group_creation_properties, H5P_DEFAULT);
    if (vtkhdf_group_id < 0) {
      Errors::errorMessage("Unable to create 'VTKHDF' group", __FILE__, __LINE__);
      status = ERROR_VALUE;
    }
  }

  //3) Header
  //---------
  if (status != ERROR_VALUE) {
    status = writeHeaderVTKHDF(vtkhdf_group_id, "UnstructuredGrid");
  }

  //4) Nodes output
  //-----------------
  std::vector<double> dataset;
  size_t numberNodes{mesh->getNodes(dataset, cellsLvl)};

  // Step 1: write number of nodes
  if (status != ERROR_VALUE) {
    status = writeScalarVTKHDF(vtkhdf_group_id, "NumberOfPoints", numberNodes);
  }

  // Step 2: write points coordinates
  if (status != ERROR_VALUE) {
    status = writeUnstructuredArrayVTKHDF(vtkhdf_group_id, "Points", H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, dataset, numberNodes, DIM_OUTPUT);
  }

  //5) Cells output
  //------------------
  // If the number of nodes holds in variable of wanted output type, then the values of the connectivity array holds too
  if (m_h5_file_intType == H5T_STD_I32LE) {
    writeConnectivity<int32_t>(mesh, cellsLvl, vtkhdf_group_id, numberCells, numberNodes);
  }
  else {
    writeConnectivity<int64_t>(mesh, cellsLvl, vtkhdf_group_id, numberCells, numberNodes);
  }

  return status;
}

//***********************************************************************

void OutputVTKHDF::prepareErrorFiles() const
{
  std::string error{m_folderErrorsAndWarnings + "error.out"};
  std::string warning{m_folderErrorsAndWarnings + "warning.out"};

  if (rankCpu == 0) {
    MPI_File_delete(error.c_str(), MPI_INFO_NULL);
    MPI_File_delete(warning.c_str(), MPI_INFO_NULL);
  }
}

//***********************************************************************

void OutputVTKHDF::writeErrorInFile(std::vector<Errors>& v, const int& ErrorType, const int& iteration, bool& warn) const
{

  size_t err_count{v.size()};
  MPI_Allreduce(MPI_IN_PLACE, &err_count, 1, ECOGEN_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

  if (err_count == 0) {
    return;
  } // No error or warning in any MPI process. Exit.

  std::string fileName{m_folderErrorsAndWarnings};

  if (ErrorType == ERROR) {
    fileName += "error.out";
  }
  else if (ErrorType == WARNING) {
    fileName += "warning.out";
    warn      = true;
  }

  std::string msg;
  for (size_t e = 0; e < v.size(); e++) {
    msg += v[e].errorToString(e, ErrorType, iteration);
  }

  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD, fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  MPI_File_write_ordered(fh, msg.c_str(), msg.length(), MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);

  // Clearing the vector is required even in case of errors because if we run several tests and don't clean
  // the vector, it will contains old data.
  v.clear();

  if (ErrorType == ERROR) {
    // Some MPI rank have failed due to errors.
    std::cout << ECOGEN_color_red << msg << ECOGEN_color_reset;
    throw ErrorECOGEN("Stop code after error... not managed");
  }
}

//***********************************************************************
void OutputVTKHDF::printTree(Mesh* mesh, const std::vector<Cell*>* cellsLvl, int resumeAMRsaveFreq) const
{
  if (resumeAMRsaveFreq != 0 && (m_numFile % resumeAMRsaveFreq) == 0) {
    try {
      std::string file;
      /* Step 1: Print domain decomposition */
      if (rankCpu == 0) {
        std::ofstream fileStream;
        file = m_folderInfoMesh + createFilename(m_domainDecomposition.c_str(), -1, -1, m_numFile);
        fileStream.open(file.c_str());
        mesh->printDomainDecomposition(fileStream);
        fileStream.close();
      }
      /* Step 2: Print cell tree */
      file = m_folderInfoMesh + createFilename(m_treeStructure.c_str(), -1, -1, m_numFile);
      MPI_File fh;

      // O_TRUNC mode doesn't exist in mpiio: delete the file (and ignore error if file doesn't exist) and create a new one.
      if (rankCpu == 0) MPI_File_delete(file.c_str(), MPI_INFO_NULL);
      MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

      //Store data to write in a vector
      int cell_count{0}; // local number of data elements
      for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
        assert(static_cast<size_t>(cell_count) <= static_cast<size_t>(std::numeric_limits<int>::max()) - cellsLvl[lvl].size() &&
               "AMR cell count overflows maximal integer");
        cell_count += static_cast<int>(cellsLvl[lvl].size());
      }
      std::vector<int> dataset(cell_count);
      size_t idx{0};
      for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
        for (size_t c = 0; c < cellsLvl[lvl].size(); c++) {
          dataset[idx++] = static_cast<int>(cellsLvl[lvl][c]->getSplit());
        }
      }

      MPI_File_write_ordered(fh, &cell_count, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write_ordered(fh, dataset.data(), cell_count, MPI_INT, MPI_STATUS_IGNORE);

      MPI_File_close(&fh);
    }
    catch (ErrorECOGEN&) {
      throw;
    }
  }
}

//***********************************************************************

void OutputVTKHDF::readTree(Mesh* mesh,
                            TypeMeshContainer<Cell*>* cellsLvl,
                            TypeMeshContainer<Cell*>* cellsLvlGhost,
                            TypeMeshContainer<CellInterface*>* cellInterfacesLvl,
                            const std::vector<AddPhys*>& addPhys,
                            int& nbCellsTotalAMR)
{
  try {
    std::string file = m_folderInfoMesh + createFilename(m_treeStructure.c_str(), -1, -1, m_numFile);

    MPI_File fh;
    int err = MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    if (err != 0) {
      // Avoid segfault if file doesn't exist
      throw ErrorInput("failed to open file: " + file);
    }

    int cell_count;
    MPI_File_read_ordered(fh, &cell_count, 1, MPI_INT, MPI_STATUS_IGNORE);
    std::vector<int> dataset(cell_count);
    MPI_File_read_ordered(fh, dataset.data(), cell_count, MPI_INT, MPI_STATUS_IGNORE);

    size_t idx{0};
    for (int lvl = 0; lvl <= mesh->getLvlMax(); lvl++) {
      //Refine cells and cell interfaces
      for (unsigned int c = 0; c < cellsLvl[lvl].size(); c++) {
        if (dataset[idx++]) mesh->refineCellAndCellInterfaces(cellsLvl[lvl][c], addPhys, nbCellsTotalAMR);
      }

      // Refine ghost cells, update of persistent communications of cells lvl + 1, reconstruct arrays of cells and cell interfaces of lvl + 1
      readTree_common(mesh, cellsLvl, cellsLvlGhost, cellInterfacesLvl, addPhys, lvl);
    }
    MPI_File_close(&fh);
  }
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************
