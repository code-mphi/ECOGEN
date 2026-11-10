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

#include "Run.h"
#include "Colors.h"

using namespace tinyxml2;

//***********************************************************************

Run::Run(std::string nameCasTest, const int& number) :
  m_numTest(number),
  m_simulationName(nameCasTest),
  m_numberSolids(0),
  m_numberEos(0),
  m_numberTransports(0),
  m_MRF(-1),
  m_viscous(false),
  m_extractRefLength(false),
  m_smoothCrossSection1d(false),
  m_dt(1.e-15),
  m_physicalTime(0.),
  m_iteration(0),
  m_resumeSimulation(0),
  m_resumeAMRsaveFreq(0),
  m_recordPsat(false),
  m_warn(false)
{
  m_mesh                           = nullptr;
  m_model                          = nullptr;
  m_gradient                       = nullptr;
  m_cellsLvl                       = nullptr;
  m_cellsLvlGhost                  = nullptr;
  m_cellInterfacesLvl              = nullptr;
  m_eos                            = nullptr;
  m_symmetry                       = nullptr;
  m_globalLimiter                  = nullptr;
  m_interfaceLimiter               = nullptr;
  m_globalVolumeFractionLimiter    = nullptr;
  m_interfaceVolumeFractionLimiter = nullptr;
  m_input                          = nullptr;
  m_outPut                         = nullptr;

  // Specific for simulation restart with mesh mapping
  m_outputMeshMapping  = nullptr;
  m_restartMeshMapping = false;

  //LES statistics (used by Output::saveInfos)
  m_globalVolume     = 0.;
  m_globalUiui       = 0.;
  m_globalEnstrophy  = 0.;
  m_globalTT         = 0.;
  m_globalT          = 0.;
  m_globalThetaTheta = 0.;
  m_meanTemp         = 0.;

  m_stat.initialize();
}

//***********************************************************************

Run::~Run() {}

//***********************************************************************

void Run::initialize()
{

  //1) Reading input file (XML format)
  //----------------------------------
  std::vector<GeometricalDomain*> domains;
  std::vector<BoundCond*> boundCond;
  try {
    m_input = new Input(this);
    m_input->readInputXML(domains, boundCond, m_solidDomains);
    // Copy input files to simulation result folder
    m_outPut->copyInputFiles();
  }
  catch (ErrorInput&) {
    for (unsigned int k = 0; k < domains.size(); k++) {
      destroy(domains[k]);
    }
    throw;
  }
  TB = new Tools(m_numberPhases, m_numberSolids, m_numberTransports);

  //2) Initialization of parallel computing (also needed for 1 CPU)
  //---------------------------------------------------------------
  parallel.initialization();
  if (Ncpu > 1) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0) std::cout << "T" << m_numTest << " | Number of CPU: " << Ncpu << std::endl;
  }

  //3) Mesh data initialization
  //---------------------------
  m_mesh->assignLimits(boundCond);
  m_cellsLvl          = new TypeMeshContainer<Cell*>[m_lvlMax + 1];
  m_cellsLvlGhost     = new TypeMeshContainer<Cell*>[m_lvlMax + 1];
  m_cellInterfacesLvl = new TypeMeshContainer<CellInterface*>[m_lvlMax + 1];
  try {
    if (m_resumeSimulation > 0) {
      if (m_outPut->getType() == TypeOutput::VTK || m_outPut->getType() == TypeOutput::VTKHDF) {
        if (rankCpu == 0)
          std::cout << "T" << m_numTest << " | Resuming simulation from result file number: " << m_resumeSimulation << "..." << std::flush;
        m_outPut->readInfos();
        if (m_mesh->getType() == AMR) {
          m_stat.startAMRTime();
          try {
            if (m_resumeSimulation % m_resumeAMRsaveFreq == 0) {
              m_outPut->readDomainDecompostion(m_mesh);
            }
            else {
              throw ErrorInput("Run::resumeSimulation: AMR resume files not available");
            }
          }
          catch (ErrorInput&) {
            throw;
          }
          m_stat.endAMRTime();
        }
      }
      else {
        Errors::errorMessage("Run::resumeSimulation: Resume option only available for VTK or VTKHDF outputs");
      }
    }
    if (m_mesh->getType() == AMR) m_stat.startAMRTime();
    m_dimension = m_mesh->initializeGeometrie(
      m_cellsLvl[0], m_cellsLvlGhost[0], m_cellInterfacesLvl[0], m_resumeSimulation, m_parallelPreTreatment, simulationOrder);
    if (m_dimension == 0) {
      // Problem when initializing the geometry
      throw ErrorECOGEN("Error when initializing the problem geometry");
    }
    if (m_mesh->getType() == AMR) m_stat.endAMRTime();
  }
  catch (ErrorECOGEN&) {
    for (unsigned int k = 0; k < domains.size(); k++) {
      destroy(domains[k]);
    }
    throw;
  }

  //4) Main array initialization using model and phase number
  //---------------------------------------------------------
  m_cellsLvl[0][0]->associateExtVar(m_model, m_gradient); //Associate external variables (model, gradient method)
  for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) {
    m_cellsLvl[0][i]->allocate(m_addPhys);
  }
  for (unsigned int i = 0; i < m_cellsLvlGhost[0].size(); i++) {
    m_cellsLvlGhost[0][i]->allocate(m_addPhys);
  }

  //5) Physical data initialization: filling fluid states
  //-----------------------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) {
    m_cellsLvl[0][i]->fill(domains, m_lvlMax);
  }
  for (unsigned int i = 0; i < m_cellsLvlGhost[0].size(); i++) {
    m_cellsLvlGhost[0][i]->fill(domains, m_lvlMax);
  }
  //EOS filling
  m_cellsLvl[0][0]->allocateEos();
  //Initialize relaxation
  m_model->initializeRelaxation(m_cellsLvl[0][0]);
  //Complete fluid state with additional calculations (sound speed, energies, mixture variables, etc.)
  for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) {
    m_cellsLvl[0][i]->completeFulfillState();
  }
  //Print turbulence decay diagnostics
  if (algorithmSaadOn) {
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
      m_addPhys[pa]->computeTurbulenceParameters();
    }
  }

  //6) Allocate Sloped and buffer Cells for Riemann problems
  //--------------------------------------------------------
  int allocateSlopeLocal = 0;
  for (unsigned int i = 0; i < m_cellInterfacesLvl[0].size(); i++) {
    m_cellInterfacesLvl[0][i]->allocateSlopes(allocateSlopeLocal);
  }
  bufferCellLeft  = new Cell;
  bufferCellRight = new Cell;
  bufferCellLeft->allocate(m_addPhys);
  bufferCellRight->allocate(m_addPhys);
  domains[0]->fillIn(bufferCellLeft);
  domains[0]->fillIn(bufferCellRight);
  m_cellsLvl[0][0]->allocateSecondOrderBuffersAndGradientVectors(bufferCellLeft->getPhases(), bufferCellLeft->getMixture()); //Only for O2NS

  //7) Set immersed boundaries
  //--------------------------
  m_mesh->setImmersedBoundaries(m_cellInterfacesLvl, simulationOrder);

  //8) Intialization of persistant communications for parallel computing
  //--------------------------------------------------------------------
  m_mesh->initializePersistentCommunications(m_cellsLvl[0], simulationOrder);
  if (Ncpu > 1) {
    parallel.communicationsPrimitives(m_eos, 0);
  }

  //9) AMR initialization
  //---------------------
  if (m_mesh->getType() == AMR) m_stat.startAMRTime();
  m_mesh->procedureRaffinementInitialization(m_cellsLvl,
                                             m_cellsLvlGhost,
                                             m_cellInterfacesLvl,
                                             m_addPhys,
                                             m_nbCellsTotalAMR,
                                             domains,
                                             m_eos,
                                             m_resumeSimulation,
                                             simulationOrder,
                                             m_solidDomains);
  if (m_mesh->getType() == AMR) m_stat.endAMRTime();

  //10) Initialize MRF Riemann coupling
  //----------------------------------
  if (m_MRF != -1) {
    if (m_sources[m_MRF]->getRiemannCoupling() != false) {
      for (unsigned int i = 0; i < m_cellInterfacesLvl[0].size(); i++) {
        m_cellInterfacesLvl[0][i]->checkMrfInterface(m_sources[m_MRF]);
      }
    }
  }

  //11) Output file preparation
  //--------------------------
  m_outPut->initializeOutput(*bufferCellLeft);
  for (unsigned int c = 0; c < m_cuts.size(); c++) m_cuts[c]->initializeOutput(*bufferCellLeft);
  for (unsigned int p = 0; p < m_probes.size(); p++) m_probes[p]->initializeOutput(*bufferCellLeft);
  for (unsigned int g = 0; g < m_globalQuantities.size(); g++) m_globalQuantities[g]->initializeOutput(*bufferCellLeft);
  for (unsigned int b = 0; b < m_recordBoundaries.size(); b++) m_recordBoundaries[b]->initializeOutput(m_cellInterfacesLvl);

  //12) Resume simulation
  //----------------------
  if (m_resumeSimulation > 0) {
    try {
      this->resumeSimulation();
    }
    catch (ErrorECOGEN&) {
      for (unsigned int k = 0; k < domains.size(); k++) {
        destroy(domains[k]);
      }
      throw;
    }
  }

  //13) Mesh mapping restart
  //------------------------
  if (m_restartMeshMapping) {
    try {
      this->restartSimulationMeshMapping(domains, m_mesh);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    catch (ErrorECOGEN&) {
      for (unsigned int k = 0; k < domains.size(); k++) {
        destroy(domains[k]);
      }
      throw;
    }
  }

  for (unsigned int d = 0; d < domains.size(); d++) {
    destroy(domains[d]);
  }

  //14) Printing t0 solution
  //------------------------
  if (m_resumeSimulation == 0) {
    try {
      //Only for few test cases
      //Additional output with purpose to track the radius of a bubble over time, the maximum pressures and else.
      //To comment if not needed. Be carefull when using it, integration for bubble radius and maximum pressure at the wall are not generalized.
      //-----
      // int phase(1); //Phase wanted for total volume computation
      // m_volumePhaseK = 0.;
      // for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->computeVolumePhaseK(m_volumePhaseK, phase); }
      // if (Ncpu > 1) { parallel.computeSum(m_volumePhaseK); }
      // // m_volumePhaseK = std::pow(m_volumePhaseK, 1./3.);
      //m_pMax = new double[4];
      //m_pMaxWall = new double[4];
      //m_pMax[0] = 0.;
      //m_pMaxWall[0] = 0.;
      //for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->lookForPmax(m_pMax, m_pMaxWall); }
      //if (Ncpu > 1) { parallel.computePMax(m_pMax[0], m_pMaxWall[0]); }
      // m_massWanted = 0.; //Initial mass
      // m_alphaWanted = -1.;
      // for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->computeMass(m_massWanted, m_alphaWanted); }
      // if (Ncpu > 1) { parallel.computeSum(m_massWanted); }
      // m_massWanted = 0.9*m_massWanted; //Percentage of the initial mass we want
      // m_alphaWanted = 0.;
      //-----
      m_stat.startIOTime();
      m_outPut->initializeOutputInfos();
      if (rankCpu == 0) std::cout << "T" << m_numTest << " | Printing file number: " << m_outPut->getNumFile() << "... " << std::flush;
      if (m_mesh->getType() == AMR) {
        m_outPut->printTree(m_mesh, m_cellsLvl, m_resumeAMRsaveFreq);
      }
      for (unsigned int c = 0; c < m_cuts.size(); c++) m_cuts[c]->writeResults(m_mesh, m_cellsLvl);
      for (unsigned int p = 0; p < m_probes.size(); p++) {
        if (m_probes[p]->possesses()) m_probes[p]->writeResults(m_mesh, m_cellsLvl);
      }
      for (unsigned int g = 0; g < m_globalQuantities.size(); g++) {
        m_globalQuantities[g]->writeResults(m_mesh, m_cellsLvl);
      }
      m_outPut->writeResults(m_mesh, m_cellsLvl);
      m_outPut->prepareErrorFiles();
      m_stat.endIOTime();
      m_stat.updateComputationTime(); // Allow to print valid value if we call printWritingInfo at iteration 0
    }
    catch (ErrorXML&) {
      throw;
    }
    if (rankCpu == 0) m_outPut->writeInfos();
  }
}

//***********************************************************************

void Run::restartSimulationMeshMapping(std::vector<GeometricalDomain*>& domains, Mesh* mesh)
{
  // Get the number of CPUs
  int nbCpusPrevious(m_outputMeshMapping->readNbCpu());
  // Output init. required to set ref. cell used for results reading
  m_outputMeshMapping->initializeOutputMeshMapping(*bufferCellLeft);

  // Store dist. between cells -> to avoid min. dist. of each partition
  // Set big value to ensure that a lower distance is found
  std::vector<double> distMin(m_cellsLvl[0].size(), 1.e10);

  // Define prefix mesh filename
  std::string meshExtension(mesh->getMeshExtension());
  bool mappedMeshIsSeq(true);
  if (nbCpusPrevious != 1) {
    // Remove mesh file extension to add 'CPU' prefix
    m_meshFileMapped.resize(m_meshFileMapped.size() - meshExtension.size() - 1);
    mappedMeshIsSeq = false;
  }

  // Each cpu of new run builds the previous run mesh
  for (int cpu = 0; cpu < nbCpusPrevious; cpu++) {
    // Build array of cells for rough mesh
    TypeMeshContainer<Cell*>* cellsLvlMeshMapped;
    cellsLvlMeshMapped = new TypeMeshContainer<Cell*>[m_lvlMax + 1];

    // Build mesh of each cpu on the fly
    std::string filenameCpuMesh(m_meshFileMapped);
    if (!mappedMeshIsSeq) {
      filenameCpuMesh += "_CPU" + IO::toString(cpu) + ".msh";
    }
    Mesh* meshMapped = new MUSGmshV2(filenameCpuMesh, meshExtension);

    // Read the mesh to build the list of of nodes and elements
    // Attribute cells to elements but no need to define faces
    if (!mappedMeshIsSeq) {
      meshMapped->initCpuMeshParallel(cellsLvlMeshMapped[0], simulationOrder, cpu);
    }
    else {
      meshMapped->initCpuMeshSequential(cellsLvlMeshMapped[0], simulationOrder);
    }

    // Set phases, mixture etc. to cell according to the model
    for (auto cell : cellsLvlMeshMapped[0]) {
      cell->allocate(m_addPhys);
    }
    // Filling cell is required to set proper EOS to mesh mapped cells and avoid nullptr during cell copy
    for (auto cell : cellsLvlMeshMapped[0]) {
      cell->fill(domains, m_lvlMax);
    }

    // Read results of rough mesh simulation
    m_outputMeshMapping->readResults(meshMapped, cellsLvlMeshMapped, cpu);

    // Mapping value rough mesh to fine one
    double distCells(0.);
    int indexNearestCell(0);
    // Find nearest cell in rough mesh
    for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) {

      // Be careful to not copy state of nearest cell in partition in case it is not
      // the nearest in the all old mesh
      indexNearestCell = -1;

      for (unsigned int j = 0; j < cellsLvlMeshMapped[0].size(); j++) {
        distCells = m_cellsLvl[0][i]->distance(cellsLvlMeshMapped[0][j]);
        if (distMin[i] > distCells) {
          distMin[i]       = distCells;
          indexNearestCell = j;
        }
      }

      // Copy only if cell is the nearest in the current partition and the previous ones
      if (indexNearestCell >= 0) {
        // Copy primitive state of rough cell to fine one
        m_cellsLvl[0][i]->copyVec(cellsLvlMeshMapped[0][indexNearestCell]->getPhases(),
                                  cellsLvlMeshMapped[0][indexNearestCell]->getMixture(),
                                  cellsLvlMeshMapped[0][indexNearestCell]->getTransports());
        // Fill cell state
        m_cellsLvl[0][i]->fulfillStateResume();
        m_cellsLvl[0][i]->completeFulfillState();
      }
    }

    //Deallocation of rough mesh
    if (cellsLvlMeshMapped != nullptr) {
      for (auto cell : cellsLvlMeshMapped[0]) {
        delete cell;
      }
    }
    delete[] cellsLvlMeshMapped;
    delete meshMapped;
  }
  delete m_outputMeshMapping;
}

//***********************************************************************

void Run::resumeSimulation()
{
  //Reconstruct the AMR mesh if any and get physical data from resume point
  try {
    if (m_mesh->getType() == AMR) {
      if (m_resumeSimulation % m_resumeAMRsaveFreq == 0) {
        m_outPut->readTree(m_mesh, m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, m_addPhys, m_nbCellsTotalAMR);
      }
    }
    m_outPut->readResults(m_mesh, m_cellsLvl);
  }
  catch (ErrorECOGEN&) {
    throw;
  }

  //Communicate physical data between processors and complete fluid state with additional calculations
  //(sound speed, energies, mixture variables,etc.)
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      m_cellsLvl[lvl][i]->fulfillStateResume();
    }
  }
  if (Ncpu > 1) {
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      parallel.communicationsPrimitives(m_eos, lvl);
      parallel.communicationsTransports(lvl);
    }
  }
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      m_cellsLvl[lvl][i]->completeFulfillState();
    }
  }
  if (m_mesh->getType() == AMR) {
    for (int lvl = 0; lvl < m_lvlMax; lvl++) {
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
        m_cellsLvl[lvl][i]->averageChildrenInParent();
      }
    }
  }
  if (Ncpu > 1) {
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      parallel.communicationsPrimitives(m_eos, lvl);
    }
  }

  if (rankCpu == 0) {
    std::cout << " OK" << std::endl;
    std::cout << "T" << m_numTest << " | -------------------------------------------" << std::endl;
  }
}

//***********************************************************************

void Run::solver()
{
  int nbCellsTotalAMRMax = m_nbCellsTotalAMR;
  double dtMax;

  //-------------------
  //Time iterative loop
  //-------------------
  bool computeFini(false);
  bool print(false);
  double printSuivante(m_physicalTime + m_timeFreq);
  while (!computeFini) {
    //Errors checking
    try {
      this->verifyErrors();
    }
    catch (ErrorECOGEN&) {
      throw;
    }

    //------------------- INTEGRATION PROCEDURE -------------------

    //Setting cons variable to zero for spatial scheme on dU/dt: no need for time step at this point
    for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) {
      m_cellsLvl[0][i]->setToZeroConsGlobal();
    }
    dtMax      = 1.e10;
    int lvlDep = 0;

    try {
      this->integrationProcedure(m_dt, lvlDep, dtMax, m_nbCellsTotalAMR);
    }
    catch (ErrorECOGEN& e) {
      std::cout << e.infoError() << std::endl;
      throw;
    }
    //-------------------- CONTROL ITERATIONS/TIME ---------------------

    //Still alive...
    if (m_iteration != 0 && m_iteration % 1000 == 0 && rankCpu == 0) {
      m_outPut->writeProgress();
    }

    m_physicalTime   += m_dt;
    TB->physicalTime  = m_physicalTime;
    m_iteration++;
    //Managing output files printing / End of time iterative loop
    if (m_timeControlIterations) {
      if (m_iteration % m_freq == 0) {
        print = true;
      }
      if (m_iteration >= m_nbIte) {
        computeFini = true;
      }
    }
    else {
      if (m_physicalTime >= printSuivante) {
        print          = true;
        printSuivante += m_timeFreq;
      }
      if (m_physicalTime >= m_finalPhysicalTime) {
        print       = true;
        computeFini = true;
      }
    }
    //Managing Sources evolutions
    for (unsigned int s = 0; s < m_sources.size(); s++) {
      m_sources[s]->sourceEvolution(m_physicalTime);
    }

#ifdef DEBUG
    //if (rankCpu == 0) std::cout << m_iteration << ", dt = " << m_dt << std::endl;
#endif // DEBUG

    //------------------------ OUTPUT FILES PRINTING -------------------------
    nbCellsTotalAMRMax = std::max(nbCellsTotalAMRMax, m_nbCellsTotalAMR);
    m_dtNext           = m_cfl * dtMax;
    if (Ncpu > 1) {
      parallel.computeDt(m_dtNext);
    }

    m_stat.startIOTime();
    //Printing cuts data
    for (unsigned int c = 0; c < m_cuts.size(); c++) {
      if (m_cuts[c]->getNextTime() <= m_physicalTime) m_cuts[c]->writeResults(m_mesh, m_cellsLvl);
    }

    //Printing probes data
    for (unsigned int p = 0; p < m_probes.size(); p++) {
      if ((m_probes[p]->possesses()) && m_probes[p]->getNextTime() <= m_physicalTime) m_probes[p]->writeResults(m_mesh, m_cellsLvl);
    }

    //Printing boundary data
    for (unsigned int b = 0; b < m_recordBoundaries.size(); b++) {
      if (m_recordBoundaries[b]->getNextTime() <= m_physicalTime) m_recordBoundaries[b]->writeResults(m_cellInterfacesLvl);
    }

    if (print) {
      if (rankCpu == 0) std::cout << "T" << m_numTest << " | Printing file number: " << m_outPut->getNumFile() << "... " << std::flush;
      //General printings
      //Only for few test case
      //-----
      // int phase(1); //Phase wanted for total volume computation
      // m_volumePhaseK = 0.;
      // for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->computeVolumePhaseK(m_volumePhaseK, phase); }
      // if (Ncpu > 1) { parallel.computeSum(m_volumePhaseK); }
      // // m_volumePhaseK = std::pow(m_volumePhaseK, 1./3.);
      //if (Ncpu > 1) { parallel.computePMax(m_pMax[0], m_pMaxWall[0]); }
      //m_pMax[0] = 0.;
      //m_pMaxWall[0] = 0.;
      // m_alphaWanted = 1.; //Volume fraction corresponding to the wanted mass
      // double mass(0.);
      // do {
      //   m_alphaWanted -= 0.01;
      //   mass = 0.;
      //   for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->computeMass(mass, m_alphaWanted); }
      //   if (Ncpu > 1) { parallel.computeSum(mass); }
      // } while (mass < 0.999*m_massWanted && m_alphaWanted > 0.001);
      // if (m_alphaWanted < 1.e-10) m_alphaWanted = 0.;
      //-----
      if (algorithmSaadOn) {
        this->saveInfosLES(m_cellsLvl);
      }
      if (m_mesh->getType() == AMR) {
        m_outPut->printTree(m_mesh, m_cellsLvl, m_resumeAMRsaveFreq);
      }
      for (unsigned int g = 0; g < m_globalQuantities.size(); g++) {
        m_globalQuantities[g]->writeResults(m_mesh, m_cellsLvl);
      }
      m_outPut->writeResults(m_mesh, m_cellsLvl);

      m_stat.endIOTime();
      m_stat.updateComputationTime();
      if (rankCpu == 0) m_outPut->writeInfos();
      print = false;
    }
    //-------------------------- TIME STEP UPDATING --------------------------
    m_dt = m_dtNext;

  } //time iterative loop end
  MPI_Barrier(MPI_COMM_WORLD);
  if (m_mesh->getType() == AMR) {
    double localLoad(0.);
    for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) {
      m_cellsLvl[0][i]->computeLoad(localLoad, 0);
    }
    std::cout << "T" << m_numTest << " | Final local load on CPU " << rankCpu << " : " << localLoad << std::endl;
  }
}

//***********************************************************************

void Run::integrationProcedure(double& dt, int lvl, double& dtMax, int& nbCellsTotalAMR)
{
  //1) AMR Level time step determination
  double dtLvl = dt * std::pow(2., -(double)lvl);

  //2) Refinement procedure
  if (m_lvlMax > 0) {
    m_stat.startAMRTime();
    m_mesh->procedureRaffinement(m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, lvl, m_addPhys, nbCellsTotalAMR, m_eos);
    if (Ncpu > 1) {
      if (lvl == 0) {
        if (m_iteration % (static_cast<int>(1. / m_cfl / 0.6) + 1) == 0) {
          m_mesh->parallelLoadBalancingAMR(
            m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, simulationOrder, m_addPhys, m_eos, nbCellsTotalAMR, m_solidDomains);
          for (unsigned int p = 0; p < m_probes.size(); p++) {
            m_probes[p]->locateProbeInMesh(m_cellsLvl[0], m_mesh->getNumberCells()); //Locate new probes CPU after Load Balancing
          }
        }
      }
    }
    m_stat.endAMRTime();
  }

  //3) Slopes determination for second order and gradients for additional physics
  //
  //Fait ici pour avoir une mise a jour d'effectuer lors de l'execution de la procedure de niveau lvl+1
  //(donc pour les slopes plus besoin de les faire au debut de resolHyperboliqueO2)
  if (simulationOrder == SimulationOrder::second) {
    if (m_mesh->getType() != TypeM::UNS) {
      for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) {
        if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
          m_cellInterfacesLvl[lvl][i]->computeSlopes();
        }
      }
    }
    else {
      //KS//Why is this loop done only here and not as well later similarly to the Cartesian case?
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
        if (!m_cellsLvl[lvl][i]->getSplit()) {
          m_cellsLvl[lvl][i]->computeGradientsO2();
          m_cellsLvl[lvl][i]->limitGradientsO2(*m_globalLimiter);
        }
      }
    }
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      parallel.communicationsSlopes(lvl);
      if (lvl > 0) {
        parallel.communicationsSlopes(lvl - 1);
      }
      m_stat.endCommunicationTime();
    }
  }
  if (lvl < m_lvlMax) {
    if (m_numberAddPhys) {
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
        if (!m_cellsLvl[lvl][i]->getSplit()) {
          m_cellsLvl[lvl][i]->prepareAddPhys();
        }
      }
      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
          m_addPhys[pa]->communicationsAddPhys(m_dimension, lvl);
        }
        m_stat.endCommunicationTime();
      }
    }
    //4) Recursive call for level up integration procedure
    this->integrationProcedure(dt, lvl + 1, dtMax, nbCellsTotalAMR);
  }

  //5) Advancement procedure
  this->advancingProcedure(dtLvl, lvl, dtMax);

  //6) Additional calculations for AMR levels > 0
  if (lvl > 0) {
    if (simulationOrder == SimulationOrder::second) {
      for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) {
        if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
          m_cellInterfacesLvl[lvl][i]->computeSlopes();
        }
      }
      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        parallel.communicationsSlopes(lvl);
        if (lvl > 0) {
          parallel.communicationsSlopes(lvl - 1);
        }
        m_stat.endCommunicationTime();
      }
    }
    if (lvl < m_lvlMax) {
      this->integrationProcedure(dt, lvl + 1, dtMax, nbCellsTotalAMR);
    }
    this->advancingProcedure(dtLvl, lvl, dtMax);
  }
}

//***********************************************************************

void Run::advancingProcedure(double& dt, int& lvl, double& dtMax)
{
  if (simulationOrder == SimulationOrder::first) {
    //1) Finite volume scheme for hyperbolic systems (Godunov)
    this->solveHyperbolic(dt, lvl, dtMax);

    //2) Finite volume scheme for additional physics
    if (m_numberAddPhys) this->solveAdditionalPhysics(dt, lvl);

    //3) Source terms integration before relaxations
    if (m_numberSources) this->solveSourceTerms(dt, lvl);
  }
  else {                              //Second order
    if (m_StrangSplitting == false) { //Godunov MUSCL
      //1) m_cons saves for AMR--second-order combination
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
        if (!m_cellsLvl[lvl][i]->getSplit()) {
          m_cellsLvl[lvl][i]->saveCons();
        }
      }
      //2) Finite volume scheme for hyperbolic systems (Godunov MUSCL)
      this->solveHyperbolicO2(dt, lvl, dtMax);

      //3) Finite volume scheme for additional physics
      if (m_numberAddPhys) this->solveAdditionalPhysics(dt, lvl);

      //4) Source terms integration before relaxations
      if (m_numberSources) this->solveSourceTerms(dt, lvl);
    }
    else { // Strang-splitting MUSCL
      //1) m_cons saves for AMR--second-order combination
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
        if (!m_cellsLvl[lvl][i]->getSplit()) {
          m_cellsLvl[lvl][i]->saveCons();
        }
      }

      double halfDt = 0.5 * dt; //Half dt only for additional physics and source terms

      //2) Finite volume scheme for additional physics
      if (m_numberAddPhys) this->solveAdditionalPhysics(halfDt, lvl);

      //3) Source terms integration before relaxations
      if (m_numberSources) this->solveSourceTerms(halfDt, lvl);

      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        parallel.communicationsPrimitives(m_eos, lvl);
        m_stat.endCommunicationTime();
      }

      //4) m_cons recovery for AMR--second-order combination (substitute to setToZeroCons)
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
        if (!m_cellsLvl[lvl][i]->getSplit()) {
          m_cellsLvl[lvl][i]->getBackCons();
        }
      }

      //5) Slopes are updated only for O2 on Cartesian grids
      if (m_mesh->getType() != TypeM::UNS) {
        for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) {
          if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
            m_cellInterfacesLvl[lvl][i]->computeSlopes();
          }
        }
        if (Ncpu > 1) {
          m_stat.startCommunicationTime();
          parallel.communicationsSlopes(lvl);
          if (lvl > 0) {
            parallel.communicationsSlopes(lvl - 1);
          }
          m_stat.endCommunicationTime();
        }
      }

      //6) Finite volume scheme for hyperbolic systems
      this->solveHyperbolicO2(dt, lvl, dtMax);

      //7) Finite volume scheme for additional physics
      if (m_numberAddPhys) this->solveAdditionalPhysics(halfDt, lvl); //QD// Do AddPhysics before or after SourceTerms ?

      //8) Source terms integration before relaxations
      if (m_numberSources) this->solveSourceTerms(halfDt, lvl);
    }
  }

  //9) Relaxations to equilibria
  if (m_model->getRelaxations()->size() > 0) this->solveRelaxations(dt, lvl);

  //10) Averaging childs cells in mother cell (if AMR)
  if (lvl < m_lvlMax) {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      m_cellsLvl[lvl][i]->averageChildrenInParent();
    }
  }

  //11) Final communications
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl);
    m_stat.endCommunicationTime();
  }

  // Only for lagrangian probes: to locate probes after lagrangian evolution
  // -----
  // for (unsigned int p = 0; p < m_probes.size(); p++) {
  //   m_probes[p]->probeDisplacement(dt);
  //   m_probes[p]->locateProbeInMesh(m_cellsLvl[0], m_mesh->getNumberCells());
  // }
  // -----

  // Only for few test case
  // -----
  // for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->lookForPmax(m_pMax, m_pMaxWall); }
  // Compute total energy of the simulation
  // double totalEnergy(0.), internalEnergies(0.), compactionEnergies(0.), elasticEnergies(0.), momentum(0.);
  // for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
  //   for (int k = 0; k < m_numberPhases; k++) {
  //     internalEnergies += m_cellsLvl[lvl][i]->getPhase(k)->getAlpha() * m_cellsLvl[lvl][i]->getPhase(k)->getDensity() * m_cellsLvl[lvl][i]->getPhase(k)->getEnergy();
  //     compactionEnergies += m_cellsLvl[lvl][i]->getPhase(k)->getAlpha() * m_cellsLvl[lvl][i]->getPhase(k)->getDensity() * m_cellsLvl[lvl][i]->getPhase(k)->getEnergyCompaction();
  //     elasticEnergies += m_cellsLvl[lvl][i]->getPhase(k)->getAlpha() * m_cellsLvl[lvl][i]->getPhase(k)->getDensity() * m_cellsLvl[lvl][i]->getPhase(k)->getEnergyElastic();
  //     // totalEnergy += m_cellsLvl[lvl][i]->getPhase(k)->getAlpha() * m_cellsLvl[lvl][i]->getPhase(k)->getDensity() * m_cellsLvl[lvl][i]->getPhase(k)->getTotalEnergy();
  //   }
  //   momentum += 0.5 * m_cellsLvl[lvl][i]->getMixture()->getDensity() * m_cellsLvl[lvl][i]->getMixture()->getVelocity().squaredNorm();
  //   totalEnergy += m_cellsLvl[lvl][i]->getMixture()->getDensity() * m_cellsLvl[lvl][i]->getMixture()->getTotalEnergy();
  // }
  // totalEnergy = internalEnergies + momentum;
  // // std::cout << std::fixed;
  // // std::cout << std::setprecision(3) << "Total energy = " << totalEnergy << std::endl;
  // std::cout << std::setprecision(15) << "Iter = " << m_iteration << " ; Total energy = " << totalEnergy
  //           << " / " << internalEnergies + compactionEnergies + elasticEnergies + momentum << std::endl;
  // std::cout << std::setprecision(15) << "Iter = " << m_iteration << " ; Total energy = " << totalEnergy << std::endl;
  // if (isnan(totalEnergy) || isnan(internalEnergies + momentum)) exit(0);
  // -----
}

//***********************************************************************

void Run::solveHyperbolicO2(double& dt, int& lvl, double& dtMax)
{
  //1) Spatial second order scheme
  //------------------------------
  //Fluxes are determined at each cells interfaces and stored in the m_cons variableof corresponding cells. Hyperbolic maximum time step determination
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) {
    if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
      m_cellInterfacesLvl[lvl][i]->computeFlux(
        dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, *m_interfaceVolumeFractionLimiter);
    }
  }

  //2) Prediction step using slopes
  //-------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->predictionOrdre2(dt, m_symmetry);
    }
  }

  //2b) Option: Activate relaxations and optional energy correction during prediction
  //KS//FP// To implement dynamically
  if (m_model->getRelaxations()->size() > 0) {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      if (!m_cellsLvl[lvl][i]->getSplit()) {
        m_model->relaxations(m_cellsLvl[lvl][i], 0.5 * dt, vecPhasesO2);
        m_model->correctionEnergy(m_cellsLvl[lvl][i], vecPhasesO2); //Correction of energies for PUEq model
        m_cellsLvl[lvl][i]->fulfillState(vecPhasesO2);
      }
    }
  }
  //2c) Option: Activate additional physics during prediction //KS//FP// To implement dynamically
  //2d) Option: Activate source terms during prediction //KS//FP// To implement dynamically
  // if (m_sources.size() > 0) {
  //   for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
  //   if (!m_cellsLvl[lvl][i]->getSplit()) {
  //       for (unsigned int s = 0; s < m_sources.size(); s++) {
  //         m_sources[s]->integrateSourceTerms(m_cellsLvl[lvl][i], 0.5 * dt, vecPhasesO2);
  //       }
  //       //m_cellsLvl[lvl][i]->setToZeroCons();
  //     }
  //   }
  // }

  //3) m_cons recovery for AMR/second order combination (substitute to setToZeroCons)
  //---------------------------------------------------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->getBackCons();
    }
  }

  //4) vecPhasesO2 communications
  //-----------------------------
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl, vecPhasesO2);
    m_stat.endCommunicationTime();
  }

  //5) Optional new slopes determination (improves code stability)
  //--------------------------------------------------------------
  // Slopes are updated only for O2 on cartesian grids
  if (m_mesh->getType() != TypeM::UNS) {
    for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) {
      if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
        m_cellInterfacesLvl[lvl][i]->computeSlopes(vecPhasesO2);
      }
    }
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      parallel.communicationsSlopes(lvl);
      if (lvl > 0) {
        parallel.communicationsSlopes(lvl - 1);
      }
      m_stat.endCommunicationTime();
    }
  }

  //6) Spatial scheme on predicted variables
  //----------------------------------------
  //Fluxes are determined at each cells interfaces and stored in the m_cons variableof corresponding cells.
  //Hyperbolic maximum time step determination
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) {
    if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
      m_cellInterfacesLvl[lvl][i]->computeFlux(
        dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, *m_interfaceVolumeFractionLimiter, vecPhasesO2);
    }
  }

  //7) Time evolution
  //-----------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolution(dt, m_symmetry); //Obtaining Cons for the scheme on (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->buildPrim();                   //It is possible to reconstruct Prim from m_cons
      m_cellsLvl[lvl][i]->setToZeroCons();               //Setting the Cons to zero for the dU/dt : no need for time step at this point
    }
  }
}

//***********************************************************************

void Run::solveHyperbolic(double& dt, int& lvl, double& dtMax)
{
  //1) Spatial scheme
  //-----------------
  //Fluxes are determined at each cells interfaces and stored in the m_cons variableof corresponding cells. Hyperbolic maximum time step determination
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) {
    if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
      m_cellInterfacesLvl[lvl][i]->computeFlux(
        dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, *m_interfaceVolumeFractionLimiter);
    }
  }

  //2) Time evolution
  //-----------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolution(dt, m_symmetry); //Obtaining Cons for the scheme on (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->buildPrim();                   //It is possible to reconstruct Prim from m_cons
      m_cellsLvl[lvl][i]->setToZeroCons();               //Setting the Cons to zero for the dU/dt : no need for time step at this point
    }
  }
}

//***********************************************************************

void Run::solveAdditionalPhysics(double& dt, int& lvl)
{
  //1) Preparation of variables for additional (gradients computations, etc) and communications
  //-------------------------------------------------------------------------------------------
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl);
    m_stat.endCommunicationTime();
  }
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->prepareAddPhys();
    }
  }
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
      m_addPhys[pa]->communicationsAddPhys(m_dimension, lvl);
    }
    m_stat.endCommunicationTime();
  }

  //2) Additional physics fluxes determination (Surface tensions, viscosity, conductivity, ...)
  //-------------------------------------------------------------------------------------------
  //Calculation of the sum of the fluxes of the additional physical elements stored in m_cons of each cell
  for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
    for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) {
      if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
        m_cellInterfacesLvl[lvl][i]->computeFluxAddPhys(*m_addPhys[pa]);
      }
    }
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      if (!m_cellsLvl[lvl][i]->getSplit()) {
        m_cellsLvl[lvl][i]->addNonConsAddPhys(*m_addPhys[pa], m_symmetry);
      }
    }
  }

  //3) Time evolution for additional physics
  //----------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolutionAddPhys(dt); //Obtaining Cons for the scheme on (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->buildPrim();              //It is possible to reconstruct Prim from m_cons
      m_cellsLvl[lvl][i]->setToZeroCons();          //Setting the Cons to zero for the dU/dt : no need for time step at this point
    }
  }
}

//***********************************************************************

void Run::solveSourceTerms(double& dt, int& lvl)
{
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      for (unsigned int s = 0; s < m_sources.size(); s++) {
        m_sources[s]->integrateSourceTerms(m_cellsLvl[lvl][i], dt);
      }
      m_cellsLvl[lvl][i]->setToZeroCons();
    }
  }
}

//***********************************************************************

void Run::solveRelaxations(double& dt, int& lvl)
{
  //Relaxations
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_model->relaxations(m_cellsLvl[lvl][i], dt);
    }
  }
  //Reset of colour function (transports) using volume fraction
  for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
    if (m_addPhys[pa]->reinitializationActivated()) {
      m_addPhys[pa]->reinitializeColorFunction(m_cellsLvl, lvl);
      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        parallel.communicationsTransports(lvl);
        m_stat.endCommunicationTime();
      }
    }
  }
  //Prepare additional physics terms for next primitive build and/or time-step
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->prepareAddPhys();
    }
  }
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
      m_addPhys[pa]->communicationsAddPhys(m_dimension, lvl);
    }
    m_stat.endCommunicationTime();
  }
  //Correction of energies for PUEq model
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_model->correctionEnergy(m_cellsLvl[lvl][i]);
      m_cellsLvl[lvl][i]->fulfillState();
    }
  }
}

//***********************************************************************

void Run::verifyErrors()
{
  //Warnings -> continue the run
  //----------------------------
  m_outPut->writeErrorInFile(warnings, WARNING, m_iteration, m_warn);
  warnings.clear();

  //Errors -> stop the run
  //----------------------
  try {
    m_outPut->writeErrorInFile(errors, ERROR, m_iteration, m_warn);
  }
  catch (ErrorECOGEN&) {
    throw;
  }
}

//***********************************************************************

void Run::finalize()
{
  if (m_warn && rankCpu == 0) {
    std::cout << ECOGEN_color_magenta;
    std::cout << "-------------------------------------------------------\n";
    std::cout << " /!\\  Warning raised:\n      Please check '" << m_outPut->getFolderOutput() << "errorsAndWarnings/' folder.\n";
    std::cout << "-------------------------------------------------------\n";
    std::cout << ECOGEN_color_reset;
  }
  //Global desallocations (some are recursives)
  if (m_cellInterfacesLvl != nullptr) {
    for (unsigned int i = 0; i < m_cellInterfacesLvl[0].size(); i++) {
      delete m_cellInterfacesLvl[0][i];
    }
  }
  if (m_cellsLvl != nullptr) {
    for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) {
      delete m_cellsLvl[0][i];
    }
  }
  if (m_cellsLvlGhost != nullptr) {
    for (unsigned int i = 0; i < m_cellsLvlGhost[0].size(); i++) {
      delete m_cellsLvlGhost[0][i];
    }
  }
  if (m_eos != nullptr) {
    for (int i = 0; i < m_numberEos; i++) {
      delete m_eos[i];
    }
  }
  delete[] m_eos;

  //Additional physics desallocations
  for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
    delete m_addPhys[pa];
  }
  for (unsigned int s = 0; s < m_sources.size(); s++) {
    delete m_sources[s];
  }

  //Second order desallocations
  if (simulationOrder == SimulationOrder::second) {
    if (slopesPhasesLocal1) {
      for (int k = 0; k < m_numberPhases; k++) {
        if (slopesPhasesLocal1[k]) delete slopesPhasesLocal1[k];
      }
    }
    if (slopesPhasesLocal2) {
      for (int k = 0; k < m_numberPhases; k++) {
        if (slopesPhasesLocal2[k]) delete slopesPhasesLocal2[k];
      }
    }

    // In the following, the destroy and destroy_array functions delete the pointers and nullify them to prevent
    // heap use after free in case of chained test cases (multiple tests in ECOGEN.xml) and failure before arrays allocation
    destroy_array(slopesPhasesLocal1);
    destroy_array(slopesPhasesLocal2);
    destroy(slopesMixtureLocal1);
    destroy(slopesMixtureLocal2);
    destroy_array(slopesTransportLocal1);
    destroy_array(slopesTransportLocal2);

    if (m_mesh && m_mesh->getType() == TypeM::UNS) {
      if (buffPhasesMin) {
        for (int k = 0; k < m_numberPhases; k++) {
          if (buffPhasesMin[k]) delete buffPhasesMin[k];
        }
      }
      if (buffPhasesMax) {
        for (int k = 0; k < m_numberPhases; k++) {
          if (buffPhasesMax[k]) delete buffPhasesMax[k];
        }
      }
      destroy_array(buffPhasesMin);
      destroy_array(buffPhasesMax);
      destroy(buffMixtureMin);
      destroy(buffMixtureMax);

      if (buffTransportMin != nullptr) {
        destroy_array(buffTransportMin);
        destroy_array(buffTransportMax);
      }
    }
  }

  //Parallel desaloccations
  if (m_mesh != nullptr) {
    m_mesh->finalizeParallele(m_lvlMax);
  }

  //Desallocations others
  destroy(TB);
  destroy(bufferCellLeft);
  destroy(bufferCellRight);
  destroy(m_mesh);
  destroy(m_model);
  destroy(m_gradient);
  destroy(m_symmetry);
  destroy(m_globalLimiter);
  destroy(m_interfaceLimiter);
  destroy(m_globalVolumeFractionLimiter);
  destroy(m_interfaceVolumeFractionLimiter);
  destroy(m_input);
  destroy(m_outPut);
  for (unsigned int d = 0; d < m_solidDomains.size(); d++) {
    delete m_solidDomains[d];
  }
  for (unsigned int c = 0; c < m_cuts.size(); c++) {
    delete m_cuts[c];
  }
  for (unsigned int p = 0; p < m_probes.size(); p++) {
    delete m_probes[p];
  }
  for (unsigned int g = 0; g < m_globalQuantities.size(); g++) {
    delete m_globalQuantities[g];
  }
  for (unsigned int b = 0; b < m_recordBoundaries.size(); b++) {
    delete m_recordBoundaries[b];
  }

  //Desallocations AMR
  destroy_array(m_cellsLvl);
  destroy_array(m_cellInterfacesLvl);
  destroy_array(m_cellsLvlGhost);
}

//***********************************************************************

void Run::saveInfosLES(std::vector<Cell*>* cellsLvl)
{
  double localVolume(0.), localUiui(0.), localEnstrophy(0.);
  double localT(0.), localTT(0.), localThetaTheta(0.);

  for (unsigned int i = 0; i < cellsLvl[0].size(); i++) {
    double vol   = cellsLvl[0][i]->getElement()->getVolume();
    localVolume += vol;

    double u   = cellsLvl[0][i]->getVelocity().getX();
    double v   = cellsLvl[0][i]->getVelocity().getY();
    double w   = cellsLvl[0][i]->getVelocity().getZ();
    localUiui += vol * ((u * u) + (v * v) + (w * w));

    double T  = cellsLvl[0][i]->getPhase(0)->getTemperature();
    localT   += vol * T;
    localTT  += vol * (T * T);

    Coord gradU = cellsLvl[0][i]->getQPA(turbulenceNumQPA)->getGrad(1);
    Coord gradV = cellsLvl[0][i]->getQPA(turbulenceNumQPA)->getGrad(2);
    Coord gradW = cellsLvl[0][i]->getQPA(turbulenceNumQPA)->getGrad(3);

    double wx       = gradW.getY() - gradV.getZ();
    double wy       = gradU.getZ() - gradW.getX();
    double wz       = gradV.getX() - gradU.getY();
    localEnstrophy += vol * ((wx * wx) + (wy * wy) + (wz * wz));

    double theta     = gradU.getX() + gradV.getY() + gradW.getZ();
    localThetaTheta += vol * (theta * theta);
  }

  // MPI Reduction
  MPI_Allreduce(&localVolume, &m_globalVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&localUiui, &m_globalUiui, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&localEnstrophy, &m_globalEnstrophy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&localTT, &m_globalTT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&localT, &m_globalT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&localThetaTheta, &m_globalThetaTheta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  m_meanTemp = m_globalT / m_globalVolume;
}

//***********************************************************************
