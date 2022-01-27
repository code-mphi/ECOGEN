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

using namespace tinyxml2;

//***********************************************************************

Run::Run(std::string nameCasTest, const int& number) : m_numTest(number), m_simulationName(nameCasTest), m_numberEos(0),
  m_numberTransports(0), m_MRF(-1), m_smoothCrossSection1d(false), m_dt(1.e-15), m_physicalTime(0.), m_iteration(0), 
  m_restartSimulation(0), m_restartAMRsaveFreq(0)
{
  m_mesh = nullptr;
  m_model = nullptr;
  m_gradient = nullptr;
  m_cellsLvl = nullptr;
  m_cellsLvlGhost = nullptr;
  m_cellInterfacesLvl = nullptr;
  m_eos = nullptr;
  m_symmetry = nullptr;
  m_globalLimiter = nullptr;
  m_interfaceLimiter = nullptr;
  m_globalVolumeFractionLimiter = nullptr;
  m_interfaceVolumeFractionLimiter = nullptr;
  m_input = nullptr;
  m_outPut = nullptr;

  m_stat.initialize();
}

//***********************************************************************

Run::~Run(){}

//***********************************************************************

void Run::initialize()
{

  //1) Reading input file (XML format)
  //----------------------------------
  std::vector<GeometricalDomain*> domains;
  std::vector<BoundCond*> boundCond;
  try {
    m_input = new Input(this);
    m_input->readInputXML(domains, boundCond);
  }
  catch (ErrorXML &) { throw; }
  TB = new Tools(m_numberPhases, m_numberTransports);

  //2) Initialization of parallel computing (also needed for 1 CPU)
  //---------------------------------------------------------------
  parallel.initialization();
  if (Ncpu > 1){
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0) std::cout << "T" << m_numTest << " | Number of CPU: " << Ncpu << std::endl;
  }

  //3) Mesh data initialization
  //---------------------------
  m_mesh->attributLimites(boundCond);
  m_cellsLvl = new TypeMeshContainer<Cell*>[m_lvlMax + 1];
  m_cellsLvlGhost = new TypeMeshContainer<Cell*>[m_lvlMax + 1];
  m_cellInterfacesLvl = new TypeMeshContainer<CellInterface*>[m_lvlMax + 1];
  try {
    if (m_restartSimulation > 0) {
      if (m_outPut->getType() == TypeOutput::XML) {
        if (rankCpu == 0) std::cout << "T" << m_numTest << " | Restarting simulation from result file number: " << m_restartSimulation << "..." << std::endl;
        m_outPut->readInfos();
        if (m_mesh->getType() == AMR) {
          if (m_restartSimulation % m_restartAMRsaveFreq == 0) {
            m_outPut->readDomainDecompostion(m_mesh);
          }
          else {
            Errors::errorMessage("Run::restartSimulation: Restart files not available");
          }
        }
      }
      else {
        Errors::errorMessage("Run::restartSimulation: Restart option only available for XML output");
      }
    }
    m_dimension = m_mesh->initializeGeometrie(m_cellsLvl[0], m_cellsLvlGhost[0], m_cellInterfacesLvl[0], m_restartSimulation, m_parallelPreTreatment, m_order);
  }
  catch (ErrorECOGEN &) { throw; }

  //4) Main array initialization using model and phase number
  //---------------------------------------------------------
  m_cellsLvl[0][0]->associateExtVar(m_model, m_gradient);    //Associate external variables (model, gradient method)
  for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) { m_cellsLvl[0][i]->allocate(m_addPhys); }
  for (unsigned int i = 0; i < m_cellsLvlGhost[0].size(); i++) { m_cellsLvlGhost[0][i]->allocate(m_addPhys); }

  //5) Physical data initialization: filling fluid states
  //-----------------------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) { m_cellsLvl[0][i]->fill(domains, m_lvlMax); }
  for (unsigned int i = 0; i < m_cellsLvlGhost[0].size(); i++) { m_cellsLvlGhost[0][i]->fill(domains, m_lvlMax); }
  //EOS filling
  m_cellsLvl[0][0]->allocateEos();
  //Complete fluid state with additional calculations (sound speed, energies, mixture variables, etc.)
  for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) { m_cellsLvl[0][i]->completeFulfillState(); }

  //6) Allocate Sloped and buffer Cells for Riemann problems
  //--------------------------------------------------------
  int allocateSlopeLocal = 0;
  for (int i = 0; i < m_mesh->getNumberFaces(); i++) { m_cellInterfacesLvl[0][i]->allocateSlopes(allocateSlopeLocal); }
  bufferCellLeft = new Cell; bufferCellRight = new Cell;
  bufferCellLeft->allocate(m_addPhys);
  bufferCellRight->allocate(m_addPhys);
  domains[0]->fillIn(bufferCellLeft);
  domains[0]->fillIn(bufferCellRight);

  //7) Intialization of persistant communications for parallel computing
  //--------------------------------------------------------------------
  m_mesh->initializePersistentCommunications(m_cellsLvl[0], m_order);
  if (Ncpu > 1) { parallel.communicationsPrimitives(m_eos, 0); }
  
  //8) AMR initialization
  //---------------------
  m_mesh->procedureRaffinementInitialization(m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, m_addPhys, m_nbCellsTotalAMR, domains, m_eos, m_restartSimulation, m_order);

  for (unsigned int d = 0; d < domains.size(); d++) { delete domains[d]; }

  //9) Output file preparation
  //--------------------------
  m_outPut->initializeOutput(*bufferCellLeft);
  for (unsigned int c = 0; c < m_cuts.size(); c++) m_cuts[c]->initializeOutput(*bufferCellLeft);
  for (unsigned int p = 0; p < m_probes.size(); p++) m_probes[p]->initializeOutput(*bufferCellLeft);
  for (unsigned int g = 0; g < m_globalQuantities.size(); g++) m_globalQuantities[g]->initializeOutput(*bufferCellLeft);
  for (unsigned int b = 0; b < m_recordBoundaries.size(); b++) m_recordBoundaries[b]->initializeOutput(m_cellInterfacesLvl);

  //10) Restart simulation
  //----------------------
  if (m_restartSimulation > 0) {
    try { this->restartSimulation(); }
    catch (ErrorECOGEN &) { throw; }
  }
  
  //11) Printing t0 solution
  //------------------------
  if (m_restartSimulation == 0) {
    try {
      //Only for few test cases
      //Additional output with purpose to track the radius of a bubble over time, the maximum pressures and else.
      //To comment if not needed. Be carefull when using it, integration for bubble radius and maximum pressure at the wall are not generalized.
      //-----
      // int phase(1); //Phase wanted for total volume computation
      // m_volumePhaseK = 0.;
      // for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->computeVolumePhaseK(m_volumePhaseK, phase); }
      // if (Ncpu > 1) { parallel.computeSum(m_volumePhaseK); }
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
      m_outPut->initializeOutputInfos();
      if (rankCpu == 0) m_outPut->writeInfos();
      m_outPut->saveInfosMailles();
      if (m_mesh->getType() == AMR) m_outPut->printTree(m_mesh, m_cellsLvl, m_restartAMRsaveFreq);
      for (unsigned int c = 0; c < m_cuts.size(); c++) m_cuts[c]->writeResults(m_mesh, m_cellsLvl);
      for (unsigned int p = 0; p < m_probes.size(); p++) { if (m_probes[p]->possesses()) m_probes[p]->writeResults(m_mesh, m_cellsLvl); }
      for (unsigned int g = 0; g < m_globalQuantities.size(); g++) { m_globalQuantities[g]->writeResults(m_mesh, m_cellsLvl); }
      m_outPut->writeResults(m_mesh, m_cellsLvl);
      Errors::prepareErrorFiles(m_outPut->getFolderOutput());
    }
    catch (ErrorXML &) { throw; }
    if (rankCpu == 0) std::cout << " OK" << std::endl;
  }
}

//***********************************************************************

void Run::restartSimulation()
{
  std::ifstream fileStream;

  //Reconstruct the AMR mesh if any and get physical data from restart point
  try {
    if (m_mesh->getType() == AMR) {
      if (m_restartSimulation % m_restartAMRsaveFreq == 0) {
        m_outPut->readTree(m_mesh, m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, m_addPhys, m_nbCellsTotalAMR);
      }
    }
    m_outPut->readResults(m_mesh, m_cellsLvl);
  }
  catch (ErrorECOGEN &) { fileStream.close(); throw; }
  fileStream.close();

  //Communicate physical data between processors and complete fluid state with additional calculations (sound speed, energies, mixture variables, etc.)
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->fulfillStateRestart(); }
  }
  if (Ncpu > 1) {
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      parallel.communicationsPrimitives(m_eos, lvl);
      parallel.communicationsTransports(lvl);
    }
  }
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->completeFulfillState(); }
  }
  if (m_mesh->getType() == AMR) {
    for (int lvl = 0; lvl < m_lvlMax; lvl++) {
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->averageChildrenInParent(); }
    }
  }
  if (Ncpu > 1) {
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) { parallel.communicationsPrimitives(m_eos, lvl); }
  }

  if (rankCpu == 0) std::cout << " OK" << std::endl;
}

//***********************************************************************

void Run::solver()
{
  int nbCellsTotalAMRMax = m_nbCellsTotalAMR;
  double dtMax;

  //-------------------
  //Time iterative loop
  //-------------------
  bool computeFini(false); bool print(false);
  double printSuivante(m_physicalTime+m_timeFreq);
  while (!computeFini) {
    //Errors checking
    try {
      this->verifyErrors();
    }
    catch (ErrorECOGEN &) { throw; }
		
    //------------------- INTEGRATION PROCEDURE -------------------

    //Setting cons variable to zero for spatial scheme on dU/dt: no need for time step at this point
    for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) { m_cellsLvl[0][i]->setToZeroConsGlobal(); }
    dtMax = 1.e10;
    int lvlDep = 0;
    this->integrationProcedure(m_dt, lvlDep, dtMax, m_nbCellsTotalAMR);
    
    //-------------------- CONTROL ITERATIONS/TIME ---------------------

    //Still alive...
    if (m_iteration !=0 && m_iteration % 1000 == 0 && rankCpu == 0) { 
      std::cout << "T" << m_numTest << " | Iteration " << m_iteration << " / Timestep " << m_dt << " / Progress " << m_physicalTime / m_finalPhysicalTime * 100. << "%" << std::endl;
    }

    m_physicalTime += m_dt;
    m_iteration++;
    //Managing output files printing / End of time iterative loop
    if (m_controleIterations) {
      if (m_iteration%m_freq == 0) { print = true; }
      if (m_iteration >= m_nbIte) { computeFini = true; }
    }
    else {
      if (m_physicalTime >= printSuivante) { print = true; printSuivante += m_timeFreq; }
      if (m_physicalTime >= m_finalPhysicalTime) { print = true; computeFini = true; }
    }
    //Managing Sources evolutions
    for (unsigned int s = 0; s < m_sources.size(); s++) { m_sources[s]->sourceEvolution(m_physicalTime); }

    #ifdef DEBUG
      if (rankCpu == 0) std::cout << m_iteration << ", dt = " << m_dt << std::endl;
    #endif // DEBUG


    //------------------------ OUTPUT FILES PRINTING -------------------------
    nbCellsTotalAMRMax = std::max(nbCellsTotalAMRMax, m_nbCellsTotalAMR);
    m_dtNext = m_cfl * dtMax;
    if (Ncpu > 1) { parallel.computeDt(m_dtNext); }
    if (print) {
      m_stat.updateComputationTime();
      //General printings
      //Only for few test case
      //-----
      // int phase(1); //Phase wanted for total volume computation
      // m_volumePhaseK = 0.;
      // for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->computeVolumePhaseK(m_volumePhaseK, phase); }
      // if (Ncpu > 1) { parallel.computeSum(m_volumePhaseK); }
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
      if (rankCpu == 0) m_outPut->writeInfos();
      m_outPut->saveInfosMailles();
      if (m_mesh->getType() == AMR) m_outPut->printTree(m_mesh, m_cellsLvl, m_restartAMRsaveFreq);
      for (unsigned int c = 0; c < m_cuts.size(); c++) { m_cuts[c]->writeResults(m_mesh, m_cellsLvl); }
      for (unsigned int g = 0; g < m_globalQuantities.size(); g++) { m_globalQuantities[g]->writeResults(m_mesh, m_cellsLvl); }
      m_outPut->writeResults(m_mesh, m_cellsLvl);
      if (rankCpu == 0) std::cout << "OK" << std::endl;
      print = false;
	}
    //Printing probes data
    for (unsigned int p = 0; p < m_probes.size(); p++) { 
      if((m_probes[p]->possesses()) && m_probes[p]->getNextTime()<=m_physicalTime) m_probes[p]->writeResults(m_mesh, m_cellsLvl);
    }

    //Printing boundary data
    for (unsigned int b = 0; b < m_recordBoundaries.size(); b++)
    {
      if (m_recordBoundaries[b]->getNextTime() <= m_physicalTime)
        m_recordBoundaries[b]->writeResults(m_cellInterfacesLvl);
    }

    //-------------------------- TIME STEP UPDATING --------------------------
    m_dt = m_dtNext;

  } //time iterative loop end
  if (rankCpu == 0) std::cout << "T" << m_numTest << " | -------------------------------------------" << std::endl;
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
        if (m_iteration % (static_cast<int>(1./m_cfl/0.6) + 1) == 0) {
          m_mesh->parallelLoadBalancingAMR(m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, m_order, m_addPhys, m_eos, nbCellsTotalAMR);
          for (unsigned int p = 0; p < m_probes.size(); p++) {
            m_probes[p]->locateProbeInMesh(m_cellsLvl[0], m_mesh->getNumberCells()); //Locate new probes CPU after Load Balancing
          }
        }
      }
    }
    m_stat.endAMRTime();
  }

  //3) Slopes determination for second order and gradients for additional physics
  //Fait ici pour avoir une mise a jour d'effectuer lors de l'execution de la procedure de niveau lvl+1 (donc pour les slopes plus besoin de les faire au debut de resolHyperboliqueO2)
  if (m_order == "SECONDORDER") {
    for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { m_cellInterfacesLvl[lvl][i]->computeSlopes(); } }
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      parallel.communicationsSlopes(lvl);
      if (lvl > 0) { parallel.communicationsSlopes(lvl - 1); }
      m_stat.endCommunicationTime();
    }
  }
  if (lvl < m_lvlMax) {
    if (m_numberAddPhys) {
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->prepareAddPhys(); } }
      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_dimension, lvl); }
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
    if (m_order == "SECONDORDER") {
      for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { m_cellInterfacesLvl[lvl][i]->computeSlopes(); } }
      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        parallel.communicationsSlopes(lvl);
        if (lvl > 0) { parallel.communicationsSlopes(lvl - 1); }
        m_stat.endCommunicationTime();
      }
    }
    if (lvl < m_lvlMax) { this->integrationProcedure(dt, lvl + 1, dtMax, nbCellsTotalAMR); }
    this->advancingProcedure(dtLvl, lvl, dtMax);
  }
}

//***********************************************************************

void Run::advancingProcedure(double& dt, int& lvl, double& dtMax)
{
  //1) Finite volume scheme for hyperbolic systems (Godunov or MUSCL)
  if (m_order == "FIRSTORDER") { this->solveHyperbolic(dt, lvl, dtMax); }
  else { this->solveHyperbolicO2(dt, lvl, dtMax); }
  //2) Finite volume scheme for additional physics
  if (m_numberAddPhys) this->solveAdditionalPhysics(dt, lvl);
  //3) Source terms integration before relaxations
  if (m_numberSources) this->solveSourceTerms(dt, lvl);
  //4) Relaxations to equilibria
  if (m_numberPhases > 1) this->solveRelaxations(dt, lvl);
  //5) Averaging childs cells in mother cell (if AMR)
  if (lvl < m_lvlMax) { for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->averageChildrenInParent(); } }
  //6) Final communications
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl);
    m_stat.endCommunicationTime();
  }
  //Only for few test case
  //-----
  //for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->lookForPmax(m_pMax, m_pMaxWall); }
  //Compute total energy of the simulation
  // double totalEnergy(0.), internalEnergies(0.), momentum(0.);
  // for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
  //   for (int k = 0; k < m_numberPhases; k++) {
  //     internalEnergies += m_cellsLvl[lvl][i]->getPhase(k)->getAlpha() * m_cellsLvl[lvl][i]->getPhase(k)->getDensity() * m_cellsLvl[lvl][i]->getPhase(k)->getEnergy();
  //   }
  //   momentum += 0.5 * m_cellsLvl[lvl][i]->getMixture()->getDensity() * m_cellsLvl[lvl][i]->getMixture()->getVelocity().squaredNorm();
  //   totalEnergy += m_cellsLvl[lvl][i]->getMixture()->getDensity() * m_cellsLvl[lvl][i]->getMixture()->getTotalEnergy();
  // }
  // // totalEnergy = internalEnergies + momentum;
  // // std::cout << std::fixed;
  // // std::cout << std::setprecision(3) << "Total energy = " << totalEnergy << std::endl;
  // std::cout << std::setprecision(15) << "Iter = " << m_iteration << " ; Total energy = " << totalEnergy << " / " << internalEnergies + momentum << std::endl;
  // if (isnan(totalEnergy) || isnan(internalEnergies + momentum)) exit(0);
  //-----
}

//***********************************************************************

void Run::solveHyperbolicO2(double& dt, int& lvl, double& dtMax)
{
  //1) m_cons saves for AMR/second order combination
  //------------------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->saveCons(); } }

  //2) Spatial second order scheme
  //------------------------------
  //Fluxes are determined at each cells interfaces and stored in the m_cons variableof corresponding cells. Hyperbolic maximum time step determination
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { m_cellInterfacesLvl[lvl][i]->computeFlux(dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, *m_interfaceVolumeFractionLimiter); } }

  //3)Prediction step using slopes
  //------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->predictionOrdre2(dt, m_symmetry); } }
  //3b) Option: Activate relaxation during prediction //KS//FP// To implement dynamically
  //3c) Option: Activate additional physics during prediction //KS//FP// To implement dynamically
  //3d) Option: Activate source terms during prediction //KS//FP// To implement dynamically

  //4) m_cons recovery for AMR/second order combination (substitute to setToZeroCons)
  //---------------------------------------------------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->recuperationCons(); } }

  //5) vecPhasesO2 communications
  //-----------------------------
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl, vecPhasesO2);
    m_stat.endCommunicationTime();
  }

  //6) Optional new slopes determination (improves code stability)
  //--------------------------------------------------------------
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { m_cellInterfacesLvl[lvl][i]->computeSlopes(vecPhasesO2); } }
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsSlopes(lvl);
    if (lvl > 0) { parallel.communicationsSlopes(lvl - 1); }
    m_stat.endCommunicationTime();
  }

  //7) Spatial scheme on predicted variables
  //----------------------------------------
  //Fluxes are determined at each cells interfaces and stored in the m_cons variableof corresponding cells. Hyperbolic maximum time step determination
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { m_cellInterfacesLvl[lvl][i]->computeFlux(dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, *m_interfaceVolumeFractionLimiter, vecPhasesO2); } }

  //8) Time evolution
  //-----------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolution(dt, m_symmetry);   //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->buildPrim();                     //On peut reconstruire Prim a partir de m_cons
      m_cellsLvl[lvl][i]->setToZeroCons();                 //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
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
      m_cellInterfacesLvl[lvl][i]->computeFlux(dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, *m_interfaceVolumeFractionLimiter); 
    }
  }

  //2) Time evolution
  //-----------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolution(dt, m_symmetry);   //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->buildPrim();                     //On peut reconstruire Prim a partir de m_cons
      m_cellsLvl[lvl][i]->setToZeroCons();                 //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
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
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->prepareAddPhys(); } }
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_dimension, lvl); }
    m_stat.endCommunicationTime();
  }

  //2) Additional physics fluxes determination (Surface tensions, viscosity, conductivity, ...)
  //-------------------------------------------------------------------------------------------
  //Calcul de la somme des flux des physiques additionnelles que l on stock dans m_cons de chaque cell
  for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
    for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { m_cellInterfacesLvl[lvl][i]->computeFluxAddPhys(*m_addPhys[pa]); } }
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->addNonConsAddPhys(*m_addPhys[pa], m_symmetry); } }
  }

  //3) Time evolution for additional physics
  //----------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolutionAddPhys(dt);  //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->buildPrim();               //On peut reconstruire Prim a partir de m_cons
      m_cellsLvl[lvl][i]->setToZeroCons();           //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
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
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->prepareAddPhys(); } }
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_dimension, lvl); }
    m_stat.endCommunicationTime();
  }
  //Optional energy corrections and other relaxations
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->correctionEnergy();               //Correction of energies
      m_cellsLvl[lvl][i]->fulfillState();
      //if (m_evaporation) m_cellsLvl[lvl][i]->relaxPTMu(); //Relaxation des pressures, temperatures et potentiels chimiques //KS//FP// To implement dynamically
    }
  }
}

//***********************************************************************

void Run::verifyErrors() const
{

  //Warnings -> continue the run
  //----------------------------
  for (unsigned int e = 0; e < warnings.size(); e++) {
    //errors[e].displayError(e);
    warnings[e].writeErrorInFile(e, m_outPut->getFolderOutput(), WARNING);
  }
  warnings.clear();

  //Errors -> stop the run
  //----------------------
  try {
    bool err(false);
    if (Ncpu > 1) {
      err = parallel.verifyStateCPUs();      
    }
    else{
      err = errors.size();
    }
    if(err){
      if (errors.size() != 0) {
        for (unsigned int e = 0; e < errors.size(); e++) {
          errors[e].displayError(e);
          errors[e].writeErrorInFile(e, m_outPut->getFolderOutput(), ERROR);
        }
      }
      throw ErrorECOGEN("Stop code after error... not managed");
    }
  }
  catch (ErrorECOGEN &) { throw; }


}

//***********************************************************************

void Run::finalize()
{
  //Global desallocations (some are recursives)
  if (m_cellInterfacesLvl != nullptr) { for (unsigned int i = 0; i < m_cellInterfacesLvl[0].size(); i++) { delete m_cellInterfacesLvl[0][i]; } }
  if (m_cellsLvl != nullptr) { for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) { delete m_cellsLvl[0][i]; } }
  if (m_cellsLvlGhost != nullptr) { for (unsigned int i = 0; i < m_cellsLvlGhost[0].size(); i++) { delete m_cellsLvlGhost[0][i]; } }
  if (m_eos != nullptr) { for (int i = 0; i < m_numberEos; i++) { delete m_eos[i]; } }
  delete[] m_eos;

  //Additional physics desallocations
  for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { delete m_addPhys[pa]; }
  for (unsigned int s = 0; s < m_sources.size(); s++) { delete m_sources[s]; }
  
  //Second order desallocations
  if (m_order == "SECONDORDER") {
    for (int k = 0; k < m_numberPhases; k++) { delete slopesPhasesLocal1[k]; }
    for (int k = 0; k < m_numberPhases; k++) { delete slopesPhasesLocal2[k]; }
    delete[] slopesPhasesLocal1;
    delete[] slopesPhasesLocal2;
    delete slopesMixtureLocal1;
    delete slopesMixtureLocal2;
    delete[] slopesTransportLocal1;
    delete[] slopesTransportLocal2;
  }

  //Parallel desaloccations
  if (m_mesh != nullptr) { m_mesh->finalizeParallele(m_lvlMax); }

  //Desallocations others
  delete TB;
  delete bufferCellLeft; delete bufferCellRight;
  delete m_mesh;
  delete m_model;
  delete m_gradient;
  delete m_symmetry;
  delete m_globalLimiter; delete m_interfaceLimiter; delete m_globalVolumeFractionLimiter; delete m_interfaceVolumeFractionLimiter;
  delete m_input;
  delete m_outPut;
  for (unsigned int c = 0; c < m_cuts.size(); c++) { delete m_cuts[c]; }
  for (unsigned int p = 0; p < m_probes.size(); p++) { delete m_probes[p]; }
  for (unsigned int g = 0; g < m_globalQuantities.size(); g++) { delete m_globalQuantities[g]; }
  for (unsigned int b = 0; b < m_recordBoundaries.size(); b++) { delete m_recordBoundaries[b]; }

  //Desallocations AMR
  delete[] m_cellsLvl;
  delete[] m_cellInterfacesLvl;
  delete[] m_cellsLvlGhost;
}

//***********************************************************************
