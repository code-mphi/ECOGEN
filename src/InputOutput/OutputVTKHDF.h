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

#ifndef OUTPUTVTKHDF_H_
#define OUTPUTVTKHDF_H_

#include <hdf5.h>

#include "Output.h"

class OutputVTKHDF : public Output
{
  public:
    OutputVTKHDF(std::string casTest, std::string run, tinyxml2::XMLElement* element, std::string fileName, Input* entree);
    OutputVTKHDF(std::string nameRun, int fileNumberRestartMeshMapping, Input* input, const H5FD_mpio_xfer_t& transfer_mode);
    ~OutputVTKHDF() override;

    virtual void initializeSpecificOutput() override;
    virtual void writeResults(Mesh* mesh, std::vector<Cell*>* cellsLvl) override;
    virtual void readResults(Mesh* mesh, std::vector<Cell*>* cellsLvl, int cpu = -1) override;

    //! \brief   Get transfer_mode field (collective or independent hdf5 and MPI I/O mode for raw data transfer).
    virtual const void* getTransferMode() const override { return (void*)&m_transfer_mode; }

  protected:
    //! \brief   Initialize fields needed by the two constructors that are provided (m_type and m_data_transfer) and assign hdf5 properties.
    //! \details Two constructors are provided with 'OutputVTKHDF' class: a constructor to write output files and a constructor to read output files of a previous run when
    // restarting a simulation with mesh mapping. This method gather the code needed by both constructors: the output type assignement and the
    // setting of the main hdf5 properties.
    //! \param   transfer_mode hdf5 parallel mode for raw data transfer (H5FD_MPIO_COLLECTIVE or H5FD_MPIO_INDEPENDENT)
    void initializeProperties(const H5FD_mpio_xfer_t& transfer_mode);

    //! \brief   Read a scalar array using parallel HDF5
    //! \details Implements a collective read operation where each MPI rank reads the entire array
    //! \param   group_id     Identifier of hdf5 object where data will be written
    //! \param   dataset_name String containing the name of the VTKHDF dataset to write
    //! \param   scalar_array Vector to store values that are readed
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    int readScalarVTKHDF(hid_t group_id, const std::string& dataset_name, std::vector<size_t>& scalar_array) const;

    //! \brief   Read an array of dimensions \a local_dims at offset \a offsets using parallel HDF5
    //! \tparam  N            Range of array to read
    //! \param   group_id     Identifier of hdf5 object where data will be readed
    //! \param   dataset_name String containing the name of the VTKHDF dataset to read
    //! \param   local_dims   Dimensions of area to read for current MPI process
    //! \param   offsets      Offsets along each direction where starts the area to read for current MPI process
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    template <size_t N>
    int readArrayVTKHDF(
      hid_t group_id, std::string dataset_name, std::vector<double>& dataset, const hsize_t (&local_dims)[N], const hsize_t (&offsets)[N]);

    //! \brief   Read the list of ECOGEN physical data in a VTKHDF file with dimensions and offsets of area to read provided as input
    //! \tparam  N            Range of array to read
    //! \param   mesh         Mesh object of the current simulation
    //! \param   cellsLvl     Computationnal cells to be filled with results of previous simulation
    //! \param   group_id     Identifier o hdf5 object in which data will be readed
    //! \param   count_scalar Dimensions of area to read for current MPI process and for a scalar data
    //!                       (dims for vectorial data are computed from dims of scalar datas)
    //! \param   start_scalar Offsets along each direction where starts the area to read for current MPI process and for scalar data
    //!                       (offset for vectorial data are computed from scalar offsets)
    template <size_t N>
    int ReadPhysicalDataVTKHDF(
      Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t group_id, const hsize_t (&count_scalar)[N], const hsize_t (&start_scalar)[N]);

    //! \brief   Read the list of ECOGEN physical data in a VTKHDF file. Dimensions of area to parse are calculated automatically based on the mesh type (cartesian, AMR...)
    //! \param   mesh         Mesh object of the current simulation
    //! \param   cellsLvl     Computationnal cells to be filled with results of previous simulation
    //! \param   group_id     Identifier o hdf5 object in which data will be readed
    //! \param   cpu          Rank of MPI process from which we want to read data from the file (not used for cartesian meshes)
    //! \param   numberCells  Vector of cell count per MPI process when file was saved. Readed from the 'NumberOfCells' file dataset (not used for cartesian meshes)
    int ReadPhysicalDataVTKHDF(Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t group_id, int cpu, const std::vector<size_t>& numberCells);

    std::string createFilenameVTKHDF(const char* name, int numFile = -1, std::string nameVariable = "defaut") const;

    //! \brief   Create or open a hdf5 file
    //! \param   file_id         file identifier
    //! \param   parallel_access true to open file with parallel access, false for serial access
    //! \param   create          true to create a new file, false to reopen an existing file in order to add data
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    int openFileVTKHDF(hid_t& file_id, bool parallel_access, bool create) const;

    //! \brief   Write a Paraview .series file for time series visualization of VTKHDF outputs.
    //! \param   numFile Index of maximal timestep value to list (an empty list is saved if \p numFile is equal to -1).
    void writeTimeseriesVTKHDF(int numFile);

    //! \brief   Write a XDMF file to visualize cartesian VTKHDF outputs while VTKHDF file format is not yet specified.
    //! \param   numFile Index of maximal timestep value to list (an empty list is saved if \p numFile is equal to -1).
    void writeXDMF(Mesh* mesh, int numFile);

    //! \brief   Write data array at XDMF format
    //! \param   dataFile     Filename
    //! \param   dataNam      Name of data to write
    //! \param   dataType     Type ('Scalar' or 'Vector') of data to write
    //! \param   numberCellsX Number of space steps in X direction
    //! \param   numberCellsY Number of space steps in Y direction
    //! \param   numberCellsZ Number of space steps in Z direction
    //! \param   dimData      Data dimension (1 for scalars, 3 for vectors)
    static void writeDataXDMF(std::ofstream& fileStream,
                              const std::string& dataFile,
                              const std::string& dataName,
                              const std::string& dataType,
                              int numberCellsX,
                              int numberCellsY,
                              int numberCellsZ,
                              int dimData);

    //! \brief   Write results at VTKHDF file format for current time step
    //! \param   mesh     Mesh object of the current simulation
    //! \param   cellsLvl Computationnal cells to be filled with results of previous simulation
    void writeResultsVTKHDF(Mesh* mesh, std::vector<Cell*>* cellsLvl);

    //! \brief   Write physical data fields at VTKHDF file format in \a file hdf5 object
    //! \details For rectilinear meshes, data are saved using natural serial organization (data in files matches geometrical data representation).
    //! For unstructured meshes, data are saved by partitions in increasing order of MPI rank.
    //! \param   mesh            Mesh object of the current simulation
    //! \param   cellsLvl        Computationnal cells to be filled with results of previous simulation
    //! \param   vtkhdf_group_id Identifier of hdf5 object where data will be written
    //! \param   numberCells     Number of Cells
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    //! \note    Exceptions are not used to avoid the creation of a corrupted file
    //! \warning Rectilinear outputs are disabled for now because they are too inefficient
    int writePhysicalDataVTKHDF(hid_t& file, const Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t& vtkhdf_group_id, size_t numberCells) const;

    //! \brief   Helper to write physical data fields of rank N at VTKHDF file format in \a file hdf5 object
    //! \tparam  N               Range in hdf5 memory and file of scalar data to write
    //! \tparam  T               Precision at which to save data in file
    //! \param   mesh            Mesh object of the current simulation
    //! \param   cellsLvl        Computationnal cells to be filled with results of previous simulation
    //! \param   vtkhdf_group_id Identifier of hdf5 object where data will be written
    //! \param   numberCells     Number of Cells
    //! \param   dimsf_scalar    global dimensions of scalar data (global dimensions of vectorial data are computed from scalar ones)
    //! \param   count_scalar    local dimensions of scalar data (local dimensions of vectorial data are computed from scalar ones)
    //! \param   start_scalar    offset where to start to write scalar data for current MPI process (offsets of vectorial data are computed from scalar ones)
    //! \param   chunk_scalar    size of chunks to optimize data access efficiency (chunks of vectorial data are computed from scalar ones)
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    //! \note    Exceptions are not used to avoid the creation of a corrupted file
    //! \warning Rectilinear outputs are disabled for now because they are too inefficient
    template <size_t N, typename T>
    int writePhysicalDataVTKHDF(const Mesh* mesh,
                                std::vector<Cell*>* cellsLvl,
                                hid_t vtkhdf_group_id,
                                size_t numberCells,
                                const hsize_t (&dimsf_scalar)[N],
                                const hsize_t (&count_scalar)[N],
                                const hsize_t (&start_scalar)[N],
                                const hsize_t (&chunk_scalar)[N]) const;

    //! \brief   Write a distributed scalar array using parallel HDF5
    //! \details Implements a collective write operation where:
    //!   - Each MPI rank writes exactly one scalar value
    //!   - The resulting dataset spans all MPI ranks
    //!   - Total size of the array equals the number of MPI processes
    //!
    //! \param   group_id     Identifier of hdf5 object where data will be written
    //! \param   dataset_name String containing the name of the VTKHDF dataset to write
    //! \param   scalar_val   Value to write for current MPI rank
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    //! \note    Exceptions are not used to avoid the creation of a corrupted file
    int writeScalarVTKHDF(const hid_t& group_id, const std::string& dataset_name, const size_t& scalar_val) const;

    //! \brief   Compute size, position and offset of writing area from 'line_count' and 'column_count' data.
    //!
    //! \details Compute size, position and offset of writing area from 'line_count' and 'column_count' argument using a collective call ('Allgather')
    //! \param line_count number of line of array to write
    //! \param dimsf global dimension along range 0 of array to write (filled by this function)
    //! \param count local dimension along range 0 of array to write (filled by this function)
    //! \param start position from which local rank will start to write in global array along range 0 (filled by this function)
    static void unstructuredFirstRangeLayout(const size_t& line_count, hsize_t* dimsf, hsize_t* count, hsize_t* start);

    //! \brief   Get buffer to write at suitable format for hdf5 (float or double)
    //!
    //! \param double_buffer buffer containing data at double precision.
    //! \param float_buffer buffer to store double data converted at float format if it is the output format asked by user
    //! \return pointer to the buffer that has to be passed to hdf5.
    const void* getOutputBuffer(const std::vector<double>& double_buffer, std::vector<float>& float_buffer) const;

    template <typename T>
    int writeConnectivity(Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t& vtkhdf_group_id, size_t& numberCells, const size_t& numberNodes) const;

    //! \brief   Write a 2-dimensional array of an unstructured mesh using parallel HDF5
    //!
    //! \details Implements a collective write operation to dump arrays of unstructured data (so for unstructured and AMR meshes) where:
    //!   - array has C++ type \a T and hdf5 type \a hdf_type
    //!   - Each MPI rank writes its \a line_count lines of the array and all its columns
    //!   - The resulting dataset spans all MPI ranks
    //!   - Size, position and offset of writing area is locally computed using a collective call ('Allgather')
    //!
    //! \tparam  T               Type of \a dataset vector
    //! \param   group_id        Identifier of hdf5 object where data will be written
    //! \param   dataset_name    String containing the name of the VTKHDF dataset to write
    //! \param   hdf5_mem_dtype  Hdf5 data type to use for storage and transfer: must match architecture (LE, BE) and the size of the data in memory (int32_t, int64_t..)
    //! \param   hdf5_file_dtype Hdf5 data type saved in file: hdf5 will manage the conversion from the \p hdf5_mem_dtype to this type if needed
    //! \param   dataset         Reference to the local vector of type \a T to write
    //! \param   line_count      Size of the first dimension of the local array
    //! \param   column_count    Size of the second dimension of the local array (1 for a 1D array)
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    //! \note    - Exceptions are not used to avoid the creation of a corrupted file
    template <typename T>
    int writeUnstructuredArrayVTKHDF(const hid_t& group_id,
                                     const std::string& dataset_name,
                                     hid_t hdf_file_dtype,
                                     hid_t hdf_mem_dtype,
                                     const std::vector<T>& dataset,
                                     const size_t& line_count,
                                     size_t column_count = 1) const;

    //! \brief   Write a N-dimensional array of C++ type \a T and hdf5 type \a hdf_type using parallel HDF5
    //!
    //! \details Implements a collective write operation where:
    //!   - Each MPI rank writes the 'dataset' vector resized under the form of an array of rank N and dimensions 'local_dims'.
    //!   - The resulting dataset spans all MPI ranks and is written in a region of rank N and size 'global_dims'
    //!   - Size, position and offset of writing area are provided as argument
    //! Allow to write efficiently structured data and can be called by previous 'writeArrayVTKHDF' function.
    //!
    //! \tparam  N               Range in hdf5 memory and file of array to write
    //! \param   group_id        Identifier of hdf5 object where data will be written
    //! \param   dataset_name    String containing the name of the VTKHDF dataset to write
    //! \param   hdf5_mem_dtype  Hdf5 data type to use for storage and transfer: must match architecture (LE, BE) and the size of the data in memory (int32_t, int64_t..)
    //! \param   hdf5_file_dtype Hdf5 data type saved in file: hdf5 will manage the conversion from the \p hdf5_mem_dtype to this type if needed
    //! \param   parallel_x_properties Property for raw data transfer (collective or independent - collective writes are unavailable if file is not open in parallel)
    //! \param   dataset         Pointer to the array to write
    //! \param   global_dims     global dimensions of array to write in file
    //! \param   local_dims      local dimensions on the current MPI process of array to write in file
    //! \param   offsets         offsets at which write the array in file
    //! \param   chunk_dims      Wanted dimensions of file chunks (to optimize data access). If chunk_dims[0] is 0, chunking is disabled.
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    //! \note    - Exceptions are not used to avoid the creation of a corrupted file
    template <size_t N, typename T>
    int writeArrayVTKHDF(const hid_t& group_id,
                         const std::string& dataset_name,
                         hid_t hdf_mem_dtype,
                         hid_t hdf_file_dtype,
                         hid_t parallel_x_properties,
                         const std::vector<T>& dataset,
                         const hsize_t (&global_dims)[N],
                         const hsize_t (&local_dims)[N],
                         const hsize_t (&offsets)[N],
                         const hsize_t (&chunk_dims)[N]) const;

    //! \brief   Write array attribute in VTKHDF file
    //! \param   group_id  Identifier of hdf5 object where data will be written
    //! \param   attribute_name name of attribute to save
    //! \param   nitem number of item to write (array size)
    //! \param   array pointer to the array to write
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    int writeArrayAttribute(const hid_t& group_id, const std::string& attribute_name, const hsize_t& nitem, const int* array) const;

    //! \brief   Write header of VTKHDF file
    //! \param   group_id  Identifier of hdf5 object where data will be written
    //! \param   mesh_type String containing the name of VTK mesh Type
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    //! \note    - file is saved at format 2.4 for now
    //! \note    - Exceptions are not used to avoid the creation of a corrupted file
    int writeHeaderVTKHDF(const hid_t& group_id, const std::string& mesh_type) const;

    // Depends on mesh type

    //! \brief   Write cartesian mesh at VTKHDF file format
    //! \param   mesh            Mesh object of the current simulation
    //! \param   file            Identifier of hdf5 object where data will be written
    //! \param   vtkhdf_group_id Identifier of 'VTKHDF' group created by current function
    //! \param   numberCells     Number of cells in the mesh (filled by the current function)
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    //! \note    Temporary output while VTKHDF file format is not yet specified. For now, uses format specified on VTK discourse forum:
    // https://discourse.vtk.org/t/support-for-rectilineargrids-and-structuredgrids-in-vtkhdf-file-format/14283/6
    // An additional xmdf file is required for Paraview to interpret data.
    //! \warning Rectilinear outputs are disabled for now because they are too inefficient
    int writeMeshRectilinearVTKHDF(Mesh* mesh, hid_t& file, hid_t& vtkhdf_group_id, size_t& numberCells);

    //! \brief   Write unstructured mesh at VTKHDF file format
    //! \param   mesh            Mesh object of the current simulation
    //! \param   cellsLvl        Computationnal cells to be filled with results of previous simulation
    //! \param   file            Identifier of hdf5 object where data will be written
    //! \param   vtkhdf_group_id Identifier of 'VTKHDF' group created by current function
    //! \param   numberCells     Number of cells in the mesh (filled by the current function)
    //! \return  ERROR_VALUE if fail, SUCCESS_VALUE if success
    //! \note    Exceptions are not used to avoid the creation of a corrupted file
    int writeMeshUnstructuredVTKHDF(Mesh* mesh, std::vector<Cell*>* cellsLvl, hid_t& file, hid_t& vtkhdf_group_id, size_t& numberCells) const;

    void writeDataset(std::vector<double> /*dataset*/, std::ofstream& /*fileStream*/, TypeData /*typeData*/) const
    {
      throw ErrorECOGEN("'writeDataset' not available for VTKHDF output format: use 'writeArrayVTKHDF' instead.");
    };

    static void getDataset(std::istringstream& /*data*/, std::vector<double>& /*dataset*/)
    {
      throw ErrorECOGEN("'getDataset' not available for VTKHDF output format: use 'readArrayVTKHDF' instead.");
    };

    virtual void prepareErrorFiles() const override;
    virtual void writeErrorInFile(std::vector<Errors>& v, const int& ErrorType, const int& iteration, bool& warn) const override;
    //! \brief   Write 'meshInfo' file that stores for each cell its subdivision status (1 if subdivided, 0 if not)
    //! \details In case of VTKHDF outputs, 1 file per time step is saved at binary format.
    //! subdivision status is saved as a 'int8_t' so file can be inspected with 'od -v -t d4
    //! <meshInfo.out>' command.
    //! \param   mesh                Mesh object of the current simulation.
    //! \param   cellsLvl            Computationnal cells
    //! \param   m_resumeAMRsaveFreq File output frequency
    //! \remark  Useful only in AMR case.
    virtual void printTree(Mesh* mesh, const std::vector<Cell*>* cellsLvl, int m_resumeAMRsaveFreq) const override;

    virtual void readTree(Mesh* mesh,
                          TypeMeshContainer<Cell*>* cellsLvl,
                          TypeMeshContainer<Cell*>* cellsLvlGhost,
                          TypeMeshContainer<CellInterface*>* const cellInterfacesLvl,
                          const std::vector<AddPhys*>& addPhys,
                          int& nbCellsTotalAMR) override;

  private:
    hid_t m_h5_file_floatType; //!< Hdf5 floating type to write in file depending on the precision asked by user
    hid_t m_h5_file_intType;   //!< Hdf5 integer type to write in file depending on the precision asked by user

    hid_t m_h5_mem_floatType; //!< Hdf5 floating type in memory space: must match architecture (LE, BE) and the size of the data in memory (float, double...)
    hid_t m_h5_mem_intType; //!< Hdf5 integer type in memory space: must match architecture (LE, BE) and the size of the data in memory (int32_t, int64_t..)

    hid_t m_group_creation_properties;   //!< Properties for hdf5 group creation
    hid_t m_dataset_creation_properties; //!< Properties for hdf5 dataset creation
    hid_t m_parallel_x_properties;       //!< Properties for hdf5 raw data transfer

    H5FD_mpio_xfer_t m_transfer_mode; //!< Perform collective (default or H5FD_MPIO_COLLECTIVE) or independent (H5FD_MPIO_INDEPENDENT) parallel I/Os

    static constexpr size_t DIM_OUTPUT{3};           //!< Always write 3D-coordinates in output files
    static constexpr int VTKHDF_VERSION[2] = {2, 4}; //!< VTKHDF format version: set to version 2.4 in .cpp file
    static constexpr int ERROR_VALUE       = -1;     //!< Error value
    static constexpr int SUCCESS_VALUE     = 0;      //!< Success value
};

#endif // OUTPUTVTKHDF_H_
