//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

//! \file      Cell.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      July 30 2018

#include "Cell.h"

using namespace std;

//***********************************************************************

Cell::Cell() : m_vecPhases(0), m_mixture(0), m_cons(0), m_vecTransports(0), m_consTransports(0), m_childrenCells(0),m_element(0)
{
  m_lvl = 0;
  m_xi = 0.;
	m_split = false;
}

//***********************************************************************

Cell::Cell(int lvl) : m_vecPhases(0), m_mixture(0), m_cons(0), m_vecTransports(0), m_consTransports(0), m_childrenCells(0), m_element(0)
{
  m_lvl = lvl;
  m_xi = 0.;
	m_split = false;
}

//***********************************************************************

Cell::~Cell()
{
  for (int k = 0; k < m_numberPhases; k++) {
    delete m_vecPhases[k];
  }
  delete[] m_vecPhases;
  delete m_mixture;
  delete[] m_vecTransports;
  delete m_cons;
  delete[] m_consTransports;
  for (unsigned int qpa = 0; qpa < m_vecQuantitiesAddPhys.size(); qpa++) {
    delete m_vecQuantitiesAddPhys[qpa];
  }
  for (unsigned int i = 0; i < m_childrenInternalBoundaries.size(); i++) {
    m_childrenInternalBoundaries[i]->finalizeFace();
    delete m_childrenInternalBoundaries[i];
  }
  m_childrenInternalBoundaries.clear();
  for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
    delete m_childrenCells[i];
  }
  m_childrenCells.clear();
}

//***********************************************************************

void Cell::addBoundary(CellInterface *bord)
{
  m_boundaries.push_back(bord);
}

//***********************************************************************

void Cell::deleteBoundary(CellInterface *bord)
{
  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    if (m_boundaries[b] == bord) { m_boundaries.erase(m_boundaries.begin() + b); }
  }
}

//***********************************************************************

void Cell::allocate(const int &numberPhases, const int &numberTransports, const std::vector<AddPhys*> &addPhys, Model *model)
{
  m_numberPhases = numberPhases;
  m_numberTransports = numberTransports;
  m_vecPhases = new Phase*[numberPhases];
  for (int k = 0; k < numberPhases; k++) {
    model->allocatePhase(&m_vecPhases[k]);
  }
  model->allocateMixture(&m_mixture);
  model->allocateCons(&m_cons,numberPhases);
  if (numberTransports > 0) {
    m_vecTransports = new Transport[numberTransports];
    m_consTransports = new Transport[numberTransports];
  }
  for (unsigned int k = 0; k < addPhys.size(); k++) {
    addPhys[k]->addQuantityAddPhys(this);
  }
  m_model = model;
}

//***********************************************************************

void Cell::allocateEos(const int &numberPhases, Model *model)
{
	model->allocateEos(*this, numberPhases);
}

//***********************************************************************

void Cell::fill(vector<GeometricalDomain*> &domains, const int &lvlMax)
{
  Coord coordinates;
  coordinates = m_element->getPosition();
  for (unsigned int geom = 0; geom < domains.size(); geom++) {
      domains[geom]->fillIn(this, m_numberPhases, m_numberTransports);
  }

  //Initial smearing of the interface. Uncomment only when needed.
  // double radius;
  // //radius = posElement.getX(); //1D
  // //radius = pow(pow(posElement.getX(), 2.) + pow(posElement.getY(), 2.), 0.5); //2D
  // radius = pow(pow(coordinates.getX(), 2.) + pow(coordinates.getY(), 2.) + pow(coordinates.getZ(), 2.), 0.5); //3D
  // double alphaAir(0.), alphaEau(0.);
  // double h(3200.e-6/50./pow(2., (double)lvlMax)); //0.04e-3
  // alphaAir = 1. / 2.*(1. - tanh((radius - 100.e-6)/1.5/h));
  // if (alphaAir > 1.) alphaAir = 1.;
  // if (alphaAir < 0.) alphaAir = 0.;
  // alphaEau = 1. - alphaAir;
  // m_vecPhases[1]->setAlpha(alphaAir);
  // m_vecPhases[0]->setAlpha(alphaEau);
  // double pressure(alphaAir*3.55e3+alphaEau*50.6625e5);
  // m_vecPhases[1]->setPressure(pressure);
  // m_vecPhases[0]->setPressure(pressure);
  // m_mixture->setPressure(pressure);
}

//***********************************************************************

void Cell::allocateAndCopyPhase(const int &phaseNumber, Phase *phase)
{
  phase->allocateAndCopyPhase(&m_vecPhases[phaseNumber]);
}

//***********************************************************************

void Cell::copyPhase(const int &phaseNumber, Phase *phase)
{
  m_vecPhases[phaseNumber]->copyPhase(*phase);
}

//***********************************************************************

void Cell::copyMixture(Mixture *mixture)
{
  m_mixture->copyMixture(*mixture);
}

//***********************************************************************

void Cell::setToZeroCons(const int &numberPhases, const int &numberTransports)
{
  m_cons->setToZero(numberPhases);
  for (int k = 0; k < numberTransports; k++) {
    m_consTransports[k].setValue(0.);
  }
}

//***********************************************************************

void Cell::setToZeroConsGlobal(const int &numberPhases, const int &numberTransports)
{
  if (!m_split) {
    m_cons->setToZero(numberPhases);
    for (int k = 0; k < numberTransports; k++) {
      m_consTransports[k].setValue(0.);
    }
  }
  else {
    for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
      m_childrenCells[i]->setToZeroConsGlobal(numberPhases, numberTransports);
    }
  }
}

//***********************************************************************

void Cell::setToZeroBufferFlux(const int &numberPhases)
{
  m_cons->setToZeroBufferFlux(numberPhases);
}

//***********************************************************************

void Cell::timeEvolution(const double &dt, const int &numberPhases, const int &numberTransports, Symmetry *symmetry, Prim type)
{
  m_cons->setBufferFlux(*this, numberPhases);                //fluxTempXXX receive conservative variables at time n : Un
  symmetry->addSymmetricTerms(this, numberPhases, type);     //m_cons is incremented by the symmetric terms from primitive variables at time n
  m_cons->multiply(dt, numberPhases);                        //m_cons is multiplied by dt
  m_cons->addFlux(1., numberPhases);                         //Adding the buffer fluxTempXXX to obtain Un+1 in m_cons
  m_cons->schemeCorrection(this, numberPhases);              //Specific correction for non conservative models

  //Same process for transport (Un construction not needed)
  for (int k = 0; k < numberTransports; k++) {
    m_consTransports[k].multiply(dt);
    m_vecTransports[k].add(m_consTransports[k].getValue());
  }
}

//***********************************************************************

void Cell::timeEvolutionAddPhys(const double &dt, const int &numberPhases, const int &numberTransports)
{
  m_cons->setBufferFlux(*this, numberPhases);                //fluxTempXXX receive conservative variables at time n : Un
  m_cons->multiply(dt, numberPhases);                        //m_cons is multiplied by dt
  m_cons->addFlux(1., numberPhases);                         //Adding the buffer fluxTempXXX to obtain Un+1 in m_cons
}

//***********************************************************************

void Cell::buildPrim(const int &numberPhases)
{
  m_cons->buildPrim(m_vecPhases, m_mixture, numberPhases); 
}

//***********************************************************************

void Cell::buildCons(const int &numberPhases)
{
  m_cons->buildCons(m_vecPhases, numberPhases, m_mixture);
}

//***********************************************************************

void Cell::correctionEnergy(const int &numberPhases)
{
  m_mixture->totalEnergyToInternalEnergy(m_vecQuantitiesAddPhys);    //Building specific internal energy form totale one
  m_cons->correctionEnergy(this, numberPhases);                      //Pressure correction
}

//***********************************************************************

void Cell::printPhasesMixture(const int &numberPhases, const int &numberTransports, ofstream &fileStream) const
{
  for (int k = 0; k < numberPhases; k++) { m_vecPhases[k]->printPhase(fileStream); }
  m_mixture->printMixture(fileStream);
  for (int k = 0; k < numberTransports; k++) { fileStream << m_vecTransports[k].getValue() << " "; }
}

//***********************************************************************

void Cell::completeFulfillState(Prim type)
{
  //Complete thermodynamical variables
  m_model->fulfillState(m_vecPhases, m_mixture, m_numberPhases, type);
  //Extended energies depending on additional physics
  this->prepareAddPhys();
  m_mixture->internalEnergyToTotalEnergy(m_vecQuantitiesAddPhys);
}

//***********************************************************************

void Cell::fulfillState(Prim type)
{
  //Complete thermodynamical variables
  m_model->fulfillState(m_vecPhases, m_mixture, m_numberPhases, type);
  //This routine is used in different configurations and a short note correspond to each one:
  //- Riemann solver: No need to reconstruct the total energy there because it isn't grabbed during the Riemann problem. The total energy is directly reconstruct there.
  //The reason is to avoid calculations on the gradients of additional physics which are not necessary and furthermore wrongly computed.
  //Note that the capillary energy is not required during the Riemann problem because the models are splitted.
  //- Parallel: No need to reconstruct the total energy there because it is already communicated.
  //Note that the gradients of additional physics would also be wrongly computed if done in the ghost cells.
  //- Relaxation or correction: The total energy doesn't have to be updated there.
}

//***********************************************************************

void Cell::localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal, const int &numberPhases, Prim type)
{
  for (int k = 0; k < numberPhases; k++) {
    m_vecPhases[k]->localProjection(normal, tangent, binormal);
  }
  m_mixture->localProjection(normal, tangent, binormal);
}

//***********************************************************************

void Cell::reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal, const int &numberPhases, Prim type)
{
  for (int k = 0; k < numberPhases; k++) {
    m_vecPhases[k]->reverseProjection(normal, tangent, binormal);
  }
  m_mixture->reverseProjection(normal, tangent, binormal);
}

//***********************************************************************

void Cell::copyVec(Phase **vecPhases, Mixture *mixture, Transport *vecTransports)
{
  for (int k = 0; k < m_numberPhases; k++) {
    m_vecPhases[k]->copyPhase(*vecPhases[k]);
  }
  m_mixture->copyMixture(*mixture);
  for (int k = 0; k < m_numberTransports; k++) {
    m_vecTransports[k] = vecTransports[k];
  }
}

//***********************************************************************

//void Cell::printCut1Dde2D(std::ofstream &fileStream, std::string variableConstanteCut, const double &valueCut, const double &dL)
//{
//  if (m_childrenCells.size() == 0) {
//    bool imprX(false), imprY(false), imprZ(false);
//    double dLsur2, position, epsilon;
//    dLsur2 = dL / pow(2., (double)m_lvl) / 2.;
//    epsilon = 1.e-3*dLsur2;
//    if (variableConstanteCut == "x") {
//      imprY = true;
//      position = m_element->getPosition().getX();
//    }
//    else if (variableConstanteCut == "y") {
//      imprX = true;
//      position = m_element->getPosition().getY();
//    }
//
//    //imprX = true;
//    //imprY = false;
//    //if (abs(m_element->getPosition().getX() - m_element->getPosition().getY()) < 1.e-6) {
//    if (abs(position - valueCut - epsilon) < dLsur2) {
//      m_element->ecritPos(fileStream, imprX, imprY, imprZ);
//      this->printPhasesMixture(m_numberPhases, m_numberTransports, fileStream);
//      fileStream << m_lvl << " " << m_xi << " ";
//      fileStream << endl;
//    }
//  }
//  else {
//    for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
//      m_childrenCells[i]->printCut1Dde2D(fileStream, variableConstanteCut, valueCut, dL);
//    }
//  }
//}
//
////***********************************************************************
//
//void Cell::printCut1Dde3D(std::ofstream &fileStream, std::string variableConstanteCut1, std::string variableConstanteCut2,
//  const double &valueCut1, const double &valueCut2, const double &dL1, const double &dL2)
//{
//  if (m_childrenCells.size() == 0) {
//    bool imprX(true), imprY(true), imprZ(true);
//    double dL1sur2, dL2sur2, position1, position2, epsilon1, epsilon2;
//    dL1sur2 = dL1 / pow(2., (double)m_lvl) / 2.;
//    dL2sur2 = dL2 / pow(2., (double)m_lvl) / 2.;
//    epsilon1 = 1.e-3*dL1sur2;
//    epsilon2 = 1.e-3*dL2sur2;
//
//    if (variableConstanteCut1 == "x") {
//      imprX = false;
//      position1 = m_element->getPosition().getX();
//    }
//    else if (variableConstanteCut1 == "y") {
//      imprY = false;
//      position1 = m_element->getPosition().getY();
//    }
//    else {
//      imprZ = false;
//      position1 = m_element->getPosition().getZ();
//    }
//
//    if (variableConstanteCut2 == "x") {
//      imprX = false;
//      position2 = m_element->getPosition().getX();
//    }
//    else if (variableConstanteCut2 == "y") {
//      imprY = false;
//      position2 = m_element->getPosition().getY();
//    }
//    else {
//      imprZ = false;
//      position2 = m_element->getPosition().getZ();
//    }
//
//    if ((abs(position1 - valueCut1 - epsilon1) < dL1sur2) && (abs(position2 - valueCut2 - epsilon2) < dL2sur2)) {
//      m_element->ecritPos(fileStream, imprX, imprY, imprZ);
//      this->printPhasesMixture(m_numberPhases, m_numberTransports, fileStream);
//      fileStream << m_lvl << " " << m_xi << " ";
//      fileStream << endl;
//    }
//  }
//  else {
//    for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
//      m_childrenCells[i]->printCut1Dde3D(fileStream, variableConstanteCut1, variableConstanteCut2, valueCut1, valueCut2, dL1, dL2);
//    }
//  }
//}

//****************************************************************************
//***************************Additional Physics*******************************
//****************************************************************************

void Cell::prepareAddPhys()
{
  for (unsigned int qpa = 0; qpa < m_vecQuantitiesAddPhys.size(); qpa++) {
    m_vecQuantitiesAddPhys[qpa]->computeQuantities(this);
  }
}

//***********************************************************************

double Cell::selectScalar(string nameVariable, int num) const
{
  //Selection scalar
  if (nameVariable == "TR") {
    return m_vecTransports[num].getValue();
    //double psi(0.), coeff(0.75);
    //psi = pow(m_vecTransports[num].getValue(), coeff) / (pow(m_vecTransports[num].getValue(), coeff) + pow((1 - m_vecTransports[num].getValue()), coeff));
    //return psi;
  }
  else if (nameVariable == "P") {
    if (m_numberPhases > 1) {
      return m_mixture->getPressure();
    }
    else {
      return m_vecPhases[num]->getPressure();
    }
  }
  else if (nameVariable == "RHO") {
    if (m_numberPhases > 1) {
      return m_mixture->getDensity();
    }
    else {
      return m_vecPhases[num]->getDensity();
    }
  }
  else if (nameVariable == "ALPHA") {
    if (m_numberPhases > 1) {
      return m_vecPhases[num]->getAlpha();
    }
    else { return 1.; }
  }
  else if (nameVariable == "u") {
    if (m_numberPhases > 1) {
      return m_mixture->getVelocity().getX();
    }
    else {
      return m_vecPhases[num]->getU();
    }
  }
  else if (nameVariable == "v") {
    if (m_numberPhases > 1) {
      return m_mixture->getVelocity().getY();
    }
    else {
      return m_vecPhases[num]->getV();
    }
  }
  else if (nameVariable == "w") {
    if (m_numberPhases > 1) {
      return m_mixture->getVelocity().getZ();
    }
    else {
      return m_vecPhases[num]->getW();
    }
  }
  else if (nameVariable == "T") {
    return m_vecPhases[num]->getTemperature();
  }
  else { Errors::errorMessage("nameVariable unknown in selectScalar (linked to QuantitiesAddPhys)"); return 0; }
}

//***********************************************************************

void Cell::setScalar(string nameVariable, const double &value, int num, int subscript)
{
  //Selection scalar
  if (nameVariable == "TR") { //transport
    m_vecTransports[num].setValue(value);
  }
  else { Errors::errorMessage("nameVariable unknown in setScalar (linked to QuantitiesAddPhys)"); }
}

//***********************************************************************

Coord Cell::selectVector(string nameVector, int num, int subscript) const
{
  //Selection vector
  if (nameVector == "QPA") { //additional physics
    return m_vecQuantitiesAddPhys[num]->getGrad(subscript);
  }
  else { Errors::errorMessage("nameVector unknown in selectVector (linked to QuantitiesAddPhys)"); return 0; }
}

//***********************************************************************

void Cell::setVector(std::string nameVector, const Coord &value, int num, int subscript)
{
  //Selection vector
  if (nameVector == "QPA") { //additional physics
    m_vecQuantitiesAddPhys[num]->setGrad(value, subscript);
  }
  else { Errors::errorMessage("nameVector unknown in setVector (linked to QuantitiesAddPhys)"); }
}

//***********************************************************************

Coord Cell::computeGradient(string nameVariable, int numPhase)
{
  int typeBord(0);
  double sommeDistanceX = 0.;
  double sommeDistanceY = 0.;
  double sommeDistanceZ = 0.;
  Coord grad(0.);               /*!< gradient vector for needed variable on the cell*/
  Coord gradProjectedFace(0.);  /*!< gradient vector for the needed variable on a boundary in the absolute system of coordinate*/             
  double gradBord(0.);          /*!< gradient for needed variable on the boundary in the face direction*/
 
  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    if (!m_boundaries[b]->getSplit()) {
      typeBord = m_boundaries[b]->whoAmI();
      if (typeBord == 0) //boundary type CellInterface/O2
      {
        // Extracting left and right variables values for each cell boundary
        // and calculus of the gradient normal to the face
        double cg = m_boundaries[b]->getCellGauche()->selectScalar(nameVariable, numPhase);
        double cd = m_boundaries[b]->getCellDroite()->selectScalar(nameVariable, numPhase);

        double distance(m_boundaries[b]->getCellGauche()->distance(m_boundaries[b]->getCellDroite()));
        gradBord = (cd - cg) / distance;

        // Projection in the absolute system of coordinate
        gradProjectedFace.setX(m_boundaries[b]->getFace()->getNormal().getX()*gradBord);
        gradProjectedFace.setY(m_boundaries[b]->getFace()->getNormal().getY()*gradBord);
        gradProjectedFace.setZ(m_boundaries[b]->getFace()->getNormal().getZ()*gradBord);

        // Summ for each boundary with ponderation using distance in each direction
        // then the cell gradient is normalized by summ of distances.
        double distanceX(m_boundaries[b]->getCellGauche()->distanceX(m_boundaries[b]->getCellDroite()));
        double distanceY(m_boundaries[b]->getCellGauche()->distanceY(m_boundaries[b]->getCellDroite()));
        double distanceZ(m_boundaries[b]->getCellGauche()->distanceZ(m_boundaries[b]->getCellDroite()));
        distanceX = abs(distanceX);
        distanceY = abs(distanceY);
        distanceZ = abs(distanceZ);

        gradProjectedFace.setXYZ(gradProjectedFace.getX()*distanceX, gradProjectedFace.getY()*distanceY, gradProjectedFace.getZ()*distanceZ);

        sommeDistanceX += distanceX;
        sommeDistanceY += distanceY;
        sommeDistanceZ += distanceZ;

        grad += gradProjectedFace;
      }
      else if (typeBord == 1) { //Boundary of type Abs
        double distanceX(this->distanceX(m_boundaries[b]));
        double distanceY(this->distanceY(m_boundaries[b]));
        double distanceZ(this->distanceZ(m_boundaries[b]));
        distanceX = abs(distanceX)*2.;
        distanceY = abs(distanceY)*2.;
        distanceZ = abs(distanceZ)*2.;
        sommeDistanceX += distanceX;
        sommeDistanceY += distanceY;
        sommeDistanceZ += distanceZ;
      }
      else if (typeBord == 6) { //Boundary of type Symetrie
        if (nameVariable == "u" || nameVariable == "v" || nameVariable == "w") {
          // Extracting left variables values
          // and calculus of the gradient normal to the face
          double cg = m_boundaries[b]->getCellGauche()->selectScalar(nameVariable, numPhase);

          double distance(this->distance(m_boundaries[b]));
          gradBord = cg / distance;

          // Multiplication of the gradient by the normal direction to guarantee symmetry
          if (nameVariable == "u") { gradBord = gradBord* m_boundaries[b]->getFace()->getNormal().getX(); }
          if (nameVariable == "v") { gradBord = gradBord* m_boundaries[b]->getFace()->getNormal().getY(); }
          if (nameVariable == "w") { gradBord = gradBord* m_boundaries[b]->getFace()->getNormal().getZ(); }

          // Projection in the absolute system of coordinate
          gradProjectedFace.setX(m_boundaries[b]->getFace()->getNormal().getX()*gradBord);
          gradProjectedFace.setY(m_boundaries[b]->getFace()->getNormal().getY()*gradBord);
          gradProjectedFace.setZ(m_boundaries[b]->getFace()->getNormal().getZ()*gradBord);

          // Summ for each boundary with ponderation using distance in each direction
          // then the cell gradient is normalized by summ of distances.
          double distanceX(this->distanceX(m_boundaries[b]));
          double distanceY(this->distanceY(m_boundaries[b]));
          double distanceZ(this->distanceZ(m_boundaries[b]));
          distanceX = abs(distanceX)*2.;
          distanceY = abs(distanceY)*2.;
          distanceZ = abs(distanceZ)*2.;

          gradProjectedFace.setXYZ(gradProjectedFace.getX()*distanceX, gradProjectedFace.getY()*distanceY, gradProjectedFace.getZ()*distanceZ);

          sommeDistanceX += distanceX;
          sommeDistanceY += distanceY;
          sommeDistanceZ += distanceZ;

          grad += gradProjectedFace;
        }
        else {
          double distanceX(this->distanceX(m_boundaries[b]));
          double distanceY(this->distanceY(m_boundaries[b]));
          double distanceZ(this->distanceZ(m_boundaries[b]));
          distanceX = abs(distanceX)*2.;
          distanceY = abs(distanceY)*2.;
          distanceZ = abs(distanceZ)*2.;
          sommeDistanceX += distanceX;
          sommeDistanceY += distanceY;
          sommeDistanceZ += distanceZ;
        }
      }
      else if (typeBord == 2) { //Boundary of type Wall
        if (nameVariable == "u" || nameVariable == "v" || nameVariable == "w") {
          // Extracting left variables values
          // and calculus of the gradient normal to the face
          double cg = m_boundaries[b]->getCellGauche()->selectScalar(nameVariable, numPhase);

          double distance(this->distance(m_boundaries[b]));
          gradBord = cg / distance;

          // Projection in the absolute system of coordinate
          gradProjectedFace.setX(m_boundaries[b]->getFace()->getNormal().getX()*gradBord);
          gradProjectedFace.setY(m_boundaries[b]->getFace()->getNormal().getY()*gradBord);
          gradProjectedFace.setZ(m_boundaries[b]->getFace()->getNormal().getZ()*gradBord);

          // Summ for each boundary with ponderation using distance in each direction
          // then the cell gradient is normalized by summ of distances.
          double distanceX(this->distanceX(m_boundaries[b]));
          double distanceY(this->distanceY(m_boundaries[b]));
          double distanceZ(this->distanceZ(m_boundaries[b]));
          distanceX = abs(distanceX)*2.;
          distanceY = abs(distanceY)*2.;
          distanceZ = abs(distanceZ)*2.;

          gradProjectedFace.setXYZ(gradProjectedFace.getX()*distanceX, gradProjectedFace.getY()*distanceY, gradProjectedFace.getZ()*distanceZ);

          sommeDistanceX += distanceX;
          sommeDistanceY += distanceY;
          sommeDistanceZ += distanceZ;

          grad += gradProjectedFace;
        }
        else {
          double distanceX(this->distanceX(m_boundaries[b]));
          double distanceY(this->distanceY(m_boundaries[b]));
          double distanceZ(this->distanceZ(m_boundaries[b]));
          distanceX = abs(distanceX)*2.;
          distanceY = abs(distanceY)*2.;
          distanceZ = abs(distanceZ)*2.;
          sommeDistanceX += distanceX;
          sommeDistanceY += distanceY;
          sommeDistanceZ += distanceZ;
        }
      }
    }
  }

  // Verifications in multiD
  if (sommeDistanceX <= 1.e-12) { sommeDistanceX = 1.; }
  if (sommeDistanceY <= 1.e-12) { sommeDistanceY = 1.; }
  if (sommeDistanceZ <= 1.e-12) { sommeDistanceZ = 1.; }

  // Final normalized gradient on the cell
  grad.setXYZ(grad.getX() / sommeDistanceX, grad.getY() / sommeDistanceY, grad.getZ() / sommeDistanceZ);

  return grad;
}

//***********************************************************************

QuantitiesAddPhys* Cell::getQPA(int &numGPA) const
{
  return m_vecQuantitiesAddPhys[numGPA];
}

//***********************************************************************

Coord Cell::getGradTk(int &numPhase, int &numAddPhys) const
{
  return m_vecQuantitiesAddPhys[numAddPhys]->getGradTk(numPhase);
}

//***********************************************************************

void Cell::setGradTk(int &numPhase, int &numAddPhys, double *buffer, int &counter)
{

  Coord grad(0.);
  grad.setX(buffer[++counter]);
  grad.setY(buffer[++counter]);
  grad.setZ(buffer[++counter]);
  m_vecQuantitiesAddPhys[numAddPhys]->setGradTk(numPhase, grad);
}

//***********************************************************************

void Cell::addNonConsAddPhys(const int &numberPhases, AddPhys &addPhys, Symmetry *symmetry)
{
  addPhys.addNonConsAddPhys(this, numberPhases);
  symmetry->addSymmetricTermsAddPhys(this, numberPhases, addPhys);
}

//***********************************************************************

void Cell::reinitializeColorFunction(const int &numTransport, const int &numPhase)
{
	m_vecTransports[numTransport].setValue(m_vecPhases[numPhase]->getAlpha());
}

//****************************************************************************
//*****************************Accessors**************************************
//****************************************************************************

int Cell::getBordsSize() const
{
  return m_boundaries.size();
}

//***********************************************************************

CellInterface* Cell::getBord(int &b)
{
  return m_boundaries[b];
}

//***********************************************************************

Phase* Cell::getPhase(const int &phaseNumber, Prim type) const
{
  return m_vecPhases[phaseNumber];
}

//***********************************************************************

Phase** Cell::getPhases(Prim type) const
{
  return m_vecPhases;
}

//***********************************************************************

Mixture* Cell::getMixture(Prim type) const
{
  return m_mixture;
}

//***********************************************************************

Flux* Cell::getCons() const
{
  return m_cons;
}

//***********************************************************************

void Cell::setCons(Flux *cons)
{
  m_cons->setCons(cons, m_numberPhases);
}

//***********************************************************************

Coord Cell::getPosition() const
{
  if (m_element != 0) { return m_element->getPosition(); }
  return 0.;
}

//***********************************************************************

Coord Cell::getSize() const
{
  return m_element->getSize();
}

//***********************************************************************

double Cell::getSizeX() const
{
  return m_element->getSizeX();
}

//***********************************************************************

double Cell::getSizeY() const
{
  return m_element->getSizeY();
}

//***********************************************************************

double Cell::getSizeZ() const
{
  return m_element->getSizeZ();
}

//***********************************************************************

void Cell::setElement(Element *element, const int &numCell)
{
  m_element = element;
  m_element->setCellAssociee(numCell);
}

//***********************************************************************

Element* Cell::getElement()
{
  return m_element;
}

//***********************************************************************

void Cell::setTransport(double value, int &numTransport, Prim type)
{
  m_vecTransports[numTransport].setValue(value);
}

//***********************************************************************

Transport& Cell::getTransport(const int &numTransport, Prim type) const
{
	return m_vecTransports[numTransport];
}

//***********************************************************************

Transport* Cell::getTransports(Prim type) const
{
	return m_vecTransports;
}

//***********************************************************************

Transport* Cell::getConsTransport(const int &numTransport) const
{
  return &m_consTransports[numTransport];
}

//***********************************************************************

void Cell::setConsTransport(double value, const int &numTransport)
{
  m_consTransports[numTransport].setValue(value);
}

//***********************************************************************

int Cell::getNumberPhases() const
{
  return m_numberPhases;
}

//***********************************************************************

int Cell::getNumberTransports() const
{
  return m_numberTransports;
}

//***********************************************************************

double Cell::getGradient()
{
  string nameVariable = "RHO";
  Coord grad(0.);
  int var = 0; //only for single phase
  grad = this->computeGradient(nameVariable, var);
  return grad.norm();
}

//***********************************************************************

Model* Cell::getModel()
{
  return m_model;
}

//***********************************************************************

Coord Cell::getVelocity()
{
  return m_model->getVelocity(this);
}

//***********************************************************************

vector<QuantitiesAddPhys*>& Cell::getVecQuantitiesAddPhys()
{
  return m_vecQuantitiesAddPhys;
}

//***********************************************************************

void Cell::printInfo() const
{
  m_element->printInfo();
}

//****************************************************************************
//******************************Distances*************************************
//****************************************************************************

double Cell::distance(Cell *c)
{
  return m_element->distance(c->getElement());
}

//***********************************************************************

double Cell::distanceX(Cell *c)
{
	return m_element->distanceX(c->getElement());
}

//***********************************************************************

double Cell::distanceY(Cell *c)
{
	return m_element->distanceY(c->getElement());
}

//***********************************************************************

double Cell::distanceZ(Cell *c)
{
	return m_element->distanceZ(c->getElement());
}

//***********************************************************************

double Cell::distance(CellInterface *b)
{
  return m_element->distance(b->getFace());
}

//***********************************************************************

double Cell::distanceX(CellInterface *b)
{
  return m_element->distanceX(b->getFace());
}

//***********************************************************************

double Cell::distanceY(CellInterface *b)
{
  return m_element->distanceY(b->getFace());
}

//***********************************************************************

double Cell::distanceZ(CellInterface *b)
{
  return m_element->distanceZ(b->getFace());
}

//***********************************************************************

bool Cell::traverseObjet(const GeometricObject &objet) const
{
  return m_element->traverseObjet(objet);
}

//****************************************************************************
//*****************************      AMR    **********************************
//****************************************************************************

void Cell::setToZeroXi()
{
  m_xi = 0.;
}

//***********************************************************************

void Cell::setToZeroConsXi()
{
  m_consXi = 0.;
}

//***********************************************************************

void Cell::timeEvolutionXi()
{
  m_xi += m_consXi;
}

//***********************************************************************

void Cell::chooseRefine(const double &xiSplit, const int &nbCellsY, const int &nbCellsZ,
  const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR)
{
  if (!m_split) {
    if (m_xi >= xiSplit) {
      if (!this->lvlNeighborTooLow()) {
        this->refineCellAndBoundaries(nbCellsY, nbCellsZ, addPhys, model);
        nbCellsTotalAMR += m_childrenCells.size() - 1;
      }
    }
  }
}

//***********************************************************************

void Cell::chooseUnrefine(const double &xiJoin, int &nbCellsTotalAMR)
{
  if (m_split) {
    bool deraffineGlobal(false);
    if (m_xi < xiJoin) {
      deraffineGlobal = true;
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //if one of my child possesses children, then unrefinement impossible
        if (m_childrenCells[i]->getNumberCellsChildren() > 0) { deraffineGlobal = false; }
      }
      //if one of my neighboor s level is too high, then unrefinement impossible
      if (deraffineGlobal) { if (this->lvlNeighborTooHigh()) { deraffineGlobal = false; }; }
    }
    //if all of my children can be unrefined, then unrefinement
    if (deraffineGlobal) {
      nbCellsTotalAMR -= m_childrenCells.size() - 1;
      this->unrefineCellAndBoundaries();
    }
  }
}

//***********************************************************************

void Cell::refineCellAndBoundaries(const int &nbCellsY, const int &nbCellsZ, const vector<AddPhys*> &addPhys, Model *model)
{
	m_split = true;

  //--------------------------------------
  //Initializations (children and dimension)
  //--------------------------------------

  double dimX(1.), dimY(0.), dimZ(0.);
  int numberCellsChildren(2);
  int dim(1);
  if (nbCellsZ != 1) {
    numberCellsChildren = 8;
    dimY = 1.;
    dimZ = 1.;
    dim = 3;
  }
  else if (nbCellsY != 1) {
    numberCellsChildren = 4;
    dimY = 1.;
    dim = 2;
  }
  CellInterface* bordRef(0);
  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    if (m_boundaries[b]->whoAmI() == 0) { bordRef = m_boundaries[b]; break; } //Boundary type CellInterface/O2
  }
  int allocateSlopeLocal = 1;
  
  //----------------
  //Cells refinement 
  //----------------

  //Mesh data initialization for children cells
  //-------------------------------------------
  double posXCellParent, posYCellParent, posZCellParent;
  double dXParent, dYParent, dZParent;
  double posXChild, posYChild, posZChild;
  posXCellParent = m_element->getPosition().getX();
  posYCellParent = m_element->getPosition().getY();
  posZCellParent = m_element->getPosition().getZ();
  dXParent = m_element->getSizeX();
  dYParent = m_element->getSizeY();
  dZParent = m_element->getSizeZ();
  double volumeCellParent, lCFLCellParent;
  volumeCellParent = m_element->getVolume();
  lCFLCellParent = m_element->getLCFL();

  for (int i = 0; i < numberCellsChildren; i++) {

    //Children cells creation
    //-----------------------
    this->createChildCell(i, m_lvl);
    m_element->creerElementChild();
    m_childrenCells[i]->setElement(m_element->getElementChild(i), i);
    m_childrenCells[i]->getElement()->setVolume(volumeCellParent / (double)numberCellsChildren);
    m_childrenCells[i]->getElement()->setLCFL(0.5*lCFLCellParent);
    m_childrenCells[i]->getElement()->setSize((1-dimX*0.5)*getSizeX(), (1 - dimY*0.5)*getSizeY(), (1 - dimZ*0.5)*getSizeZ());
    posXChild = posXCellParent + dimX*dXParent*(double)(-0.25 + 0.5 * (i % 2));
    posYChild = posYCellParent + dimY*dYParent*(double)(-0.25 + 0.5 * ((i / 2) % 2));
    posZChild = posZCellParent + dimZ*dZParent*(double)(-0.25 + 0.5 * ((i / 4) % 2));
    m_childrenCells[i]->getElement()->setPos(posXChild, posYChild, posZChild);

    //Initialization of main arrays according to model and number of phases
    //+ physical initialization: physical data for children cells
    //---------------------------------------------------------------------
    m_childrenCells[i]->allocate(m_numberPhases, m_numberTransports, addPhys, model);
    for (int k = 0; k < m_numberPhases; k++) {
      m_childrenCells[i]->copyPhase(k, m_vecPhases[k]);
    }
    m_childrenCells[i]->copyMixture(m_mixture);
    m_childrenCells[i]->getCons()->setToZero(m_numberPhases);
    for (int k = 0; k < m_numberTransports; k++) { m_childrenCells[i]->setTransport(m_vecTransports[k].getValue(), k); }
    for (int k = 0; k < m_numberTransports; k++) { m_childrenCells[i]->setConsTransport(0., k); }
    m_childrenCells[i]->setXi(m_xi);
  }

  //------------------------------
  //Internal boundaries refinement
  //------------------------------

  if (nbCellsZ == 1) {
    if (nbCellsY == 1) {

      //Case 1D
      //-------

      // |-------------|-------------| X
      // 1      0      0      1      2

      //Internal boundary child number 0 (face on X)
      //-------------------------------------------
      bordRef->creerBordChildInterne(m_lvl, &m_childrenInternalBoundaries);
      m_childrenInternalBoundaries[0]->creerFaceChild(bordRef);
      m_childrenInternalBoundaries[0]->getFace()->setNormal(1., 0., 0.);
      m_childrenInternalBoundaries[0]->getFace()->setTangent(0., 1., 0.);
      m_childrenInternalBoundaries[0]->getFace()->setBinormal(0., 0., 1.);
      m_childrenInternalBoundaries[0]->getFace()->setPos(posXCellParent, posYCellParent, posZCellParent);
      m_childrenInternalBoundaries[0]->getFace()->setSize(0., m_element->getSizeY(), m_element->getSizeZ());
      m_childrenInternalBoundaries[0]->getFace()->setSurface(m_element->getSizeY()*m_element->getSizeZ());
      m_childrenInternalBoundaries[0]->initializeGauche(m_childrenCells[0]);
      m_childrenInternalBoundaries[0]->initializeDroite(m_childrenCells[1]);
      m_childrenCells[0]->addBoundary(m_childrenInternalBoundaries[0]);
      m_childrenCells[1]->addBoundary(m_childrenInternalBoundaries[0]);

      //Attribution model and slopes
      //----------------------------
      m_childrenInternalBoundaries[0]->associeModel(model);
      m_childrenInternalBoundaries[0]->allocateSlopes(m_numberPhases, m_numberTransports, allocateSlopeLocal);
    }
    else {

      //Case 2D
      //-------
      // Y
      // |--------------|--------------|
      // |              |              |
      // |              |              |
      // |      2      1|      3       |
      // |              |              |
      // |              |              |
      // |--------------|--------------|
      // |      2       |      3       |
      // |              |              |
      // |      0      0|      1       |
      // |              |              |
      // |              |              |
      // |--------------|--------------|X

      for (int i = 0; i < 4; i++) {
        bordRef->creerBordChildInterne(m_lvl, &m_childrenInternalBoundaries);
        m_childrenInternalBoundaries[i]->creerFaceChild(bordRef);
        if (i < 2) {
          //Internal boundaries 0 and 1 (face on X)
          //-----------------------------------
          m_childrenInternalBoundaries[i]->getFace()->setNormal(1., 0., 0.);
          m_childrenInternalBoundaries[i]->getFace()->setTangent(0., 1., 0.);
          m_childrenInternalBoundaries[i]->getFace()->setBinormal(0., 0., 1.);
          m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent, posYCellParent + dYParent*(-0.25 + 0.5 * (double)i), posZCellParent);
          m_childrenInternalBoundaries[i]->getFace()->setSize(0., 0.5*m_element->getSizeY(), m_element->getSizeZ());
          m_childrenInternalBoundaries[i]->getFace()->setSurface(0.5*m_element->getSizeY()*m_element->getSizeZ());
          m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[2 * i]);
          m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[1 + 2 * i]);
          m_childrenCells[2 * i]->addBoundary(m_childrenInternalBoundaries[i]);
          m_childrenCells[1 + 2 * i]->addBoundary(m_childrenInternalBoundaries[i]);
        }
        else {
          //Internal boundaries 2 and 3 (face on Y)
          //-----------------------------------
          m_childrenInternalBoundaries[i]->getFace()->setNormal(0., 1., 0.);
          m_childrenInternalBoundaries[i]->getFace()->setTangent(-1., 0., 0.);
          m_childrenInternalBoundaries[i]->getFace()->setBinormal(0., 0., 1.);
          m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent + dXParent*(-0.25 + 0.5 * (double)(i % 2)), posYCellParent, posZCellParent);
          m_childrenInternalBoundaries[i]->getFace()->setSize(0.5*m_element->getSizeX(), 0., m_element->getSizeZ());
          m_childrenInternalBoundaries[i]->getFace()->setSurface(0.5*m_element->getSizeX()*m_element->getSizeZ());
          m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[i % 2]);
          m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[2 + i % 2]);
          m_childrenCells[i % 2]->addBoundary(m_childrenInternalBoundaries[i]);
          m_childrenCells[2 + i % 2]->addBoundary(m_childrenInternalBoundaries[i]);
        }
        //Attribution model and slopes
        //----------------------------
        m_childrenInternalBoundaries[i]->associeModel(model);
        m_childrenInternalBoundaries[i]->allocateSlopes(m_numberPhases, m_numberTransports, allocateSlopeLocal);
      }
    }
  }
  else {

    //Case 3D
    //-------

    //3 times the 2D plane on X, Y and Z with following numbering:
    //- Number au center correspond ici au number de la face (normal vers nous),
    //- Number dans le coin bas gauche pour cells devant le plan (cell droite),
    //- Number dans le coin haut droite pour cells derriere le plan (cell gauche).

    //             Selon X :                          Selon Y:                           Selon Z :
    //
    // |--------------|--------------|    |--------------|--------------|    |--------------|--------------|
    // |             6|             2|    |             4|             0|    |             2|             3|
    // |              |              |    |              |              |    |              |              |
    // |      2       |      3       |    |      2       |      3       |    |      2       |      3       |
    // |              |              |    |              |              |    |              |              |
    // |7             |3             |    |6             |2             |    |6             |7             |
    // |--------------|--------------|    |--------------|--------------|    |--------------|--------------|
    // |             4|             0|    |             5|             1|    |             0|             1|
    // |              |              |    |              |              |    |              |              |
    // |      0       |      1       |    |      0       |      1       |    |      0       |      1       |
    // |              |              |    |              |              |    |              |              |
    // |5             |1             |    |7             |3             |    |4             |5             |
    // |--------------|--------------|    |--------------|--------------|    |--------------|--------------|
    //
    //                y                              z---|                                  y
    //                |                                  |                                  |
    //            z---|                                  x                                  |---x

    //Face on X
    for (int i = 0; i < 4; i++) {
      bordRef->creerBordChildInterne(m_lvl, &m_childrenInternalBoundaries);
      m_childrenInternalBoundaries[i]->creerFaceChild(bordRef);
      m_childrenInternalBoundaries[i]->getFace()->setNormal(1., 0., 0.);
      m_childrenInternalBoundaries[i]->getFace()->setTangent(0., 1., 0.);
      m_childrenInternalBoundaries[i]->getFace()->setBinormal(0., 0., 1.);
      if (i == 0) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent, posYCellParent - 0.25*dYParent, posZCellParent + 0.25*dZParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[4]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[5]);
        m_childrenCells[4]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[5]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else if (i == 1) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent, posYCellParent - 0.25*dYParent, posZCellParent - 0.25*dZParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[0]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[1]);
        m_childrenCells[0]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[1]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else if (i == 2) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent, posYCellParent + 0.25*dYParent, posZCellParent + 0.25*dZParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[6]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[7]);
        m_childrenCells[6]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[7]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent, posYCellParent + 0.25*dYParent, posZCellParent - 0.25*dZParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[2]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[3]);
        m_childrenCells[2]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[3]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      m_childrenInternalBoundaries[i]->getFace()->setSize(0., 0.5*m_element->getSizeY(), 0.5*m_element->getSizeZ());
      m_childrenInternalBoundaries[i]->getFace()->setSurface(0.5*m_element->getSizeY()*0.5*m_element->getSizeZ());
      //Attribution model and slopes
      m_childrenInternalBoundaries[i]->associeModel(model);
      m_childrenInternalBoundaries[i]->allocateSlopes(m_numberPhases, m_numberTransports, allocateSlopeLocal);
    }

    //Face on Y
    for (int i = 4; i < 8; i++) {
      bordRef->creerBordChildInterne(m_lvl, &m_childrenInternalBoundaries);
      m_childrenInternalBoundaries[i]->creerFaceChild(bordRef);
      m_childrenInternalBoundaries[i]->getFace()->setNormal(0., 1., 0.);
      m_childrenInternalBoundaries[i]->getFace()->setTangent(-1., 0., 0.);
      m_childrenInternalBoundaries[i]->getFace()->setBinormal(0., 0., 1.);
      if (i == 4) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent + 0.25*dXParent, posYCellParent, posZCellParent + 0.25*dZParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[5]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[7]);
        m_childrenCells[5]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[7]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else if (i == 5) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent + 0.25*dXParent, posYCellParent, posZCellParent - 0.25*dZParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[1]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[3]);
        m_childrenCells[1]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[3]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else if (i == 6) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent - 0.25*dXParent, posYCellParent, posZCellParent + 0.25*dZParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[4]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[6]);
        m_childrenCells[4]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[6]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent - 0.25*dXParent, posYCellParent, posZCellParent - 0.25*dZParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[0]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[2]);
        m_childrenCells[0]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[2]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      m_childrenInternalBoundaries[i]->getFace()->setSize(0.5*m_element->getSizeX(), 0., 0.5*m_element->getSizeZ());
      m_childrenInternalBoundaries[i]->getFace()->setSurface(0.5*m_element->getSizeX()*0.5*m_element->getSizeZ());
      //Attribution model and slopes
      m_childrenInternalBoundaries[i]->associeModel(model);
      m_childrenInternalBoundaries[i]->allocateSlopes(m_numberPhases, m_numberTransports, allocateSlopeLocal);
    }

    //Face on Z
    for (int i = 8; i < 12; i++) {
      bordRef->creerBordChildInterne(m_lvl, &m_childrenInternalBoundaries);
      m_childrenInternalBoundaries[i]->creerFaceChild(bordRef);
      m_childrenInternalBoundaries[i]->getFace()->setNormal(0., 0., 1.);
      m_childrenInternalBoundaries[i]->getFace()->setTangent(1., 0., 0.);
      m_childrenInternalBoundaries[i]->getFace()->setBinormal(0., 1., 0.);
      if (i == 8) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent - 0.25*dXParent, posYCellParent - 0.25*dYParent, posZCellParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[0]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[4]);
        m_childrenCells[0]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[4]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else if (i == 9) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent + 0.25*dXParent, posYCellParent - 0.25*dYParent, posZCellParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[1]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[5]);
        m_childrenCells[1]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[5]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else if (i == 10) {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent - 0.25*dXParent, posYCellParent + 0.25*dYParent, posZCellParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[2]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[6]);
        m_childrenCells[2]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[6]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      else {
        m_childrenInternalBoundaries[i]->getFace()->setPos(posXCellParent + 0.25*dXParent, posYCellParent + 0.25*dYParent, posZCellParent);
        m_childrenInternalBoundaries[i]->initializeGauche(m_childrenCells[3]);
        m_childrenInternalBoundaries[i]->initializeDroite(m_childrenCells[7]);
        m_childrenCells[3]->addBoundary(m_childrenInternalBoundaries[i]);
        m_childrenCells[7]->addBoundary(m_childrenInternalBoundaries[i]);
      }
      m_childrenInternalBoundaries[i]->getFace()->setSize(0.5*m_element->getSizeX(), 0.5*m_element->getSizeY(), 0.);
      m_childrenInternalBoundaries[i]->getFace()->setSurface(0.5*m_element->getSizeX()*0.5*m_element->getSizeY());
      //Attribution model and slopes
      m_childrenInternalBoundaries[i]->associeModel(model);
      m_childrenInternalBoundaries[i]->allocateSlopes(m_numberPhases, m_numberTransports, allocateSlopeLocal);
    }
  }

  //------------------------------
  //External boundaries refinement
  //------------------------------

  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    if (!m_boundaries[b]->getSplit()) { m_boundaries[b]->raffineBordExterne(nbCellsY, nbCellsZ, dXParent, dYParent, dZParent, this, dim); }
  }
}

//***********************************************************************

void Cell::createChildCell(const int &num, const int &lvl)
{
  m_childrenCells.push_back(new Cell(lvl + 1));
}

//***********************************************************************

void Cell::unrefineCellAndBoundaries()
{
  //---------------------
  //Update of parent cell
  //---------------------

  this->averageChildrenInParent();

  //----------------------------------------
  //Internal children boundaries destruction
  //----------------------------------------

  for (unsigned int i = 0; i < m_childrenInternalBoundaries.size(); i++) {
    m_childrenInternalBoundaries[i]->finalizeFace();
    delete m_childrenInternalBoundaries[i];
  }
  m_childrenInternalBoundaries.clear();

  //--------------------------------------
  //External children boundaries destruction
  //--------------------------------------

  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    m_boundaries[b]->deraffineBordExterne(this);
  }

  //--------------------------
  //Children cells destruction
  //--------------------------

  for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
    delete m_childrenCells[i];
  }
  m_childrenCells.clear();
  m_element->finalizeElementsChildren();

	m_split = false;
}

//***********************************************************************

void Cell::averageChildrenInParent()
{
  int numberCellsChildren(m_childrenCells.size());
  if (numberCellsChildren > 0) {
    //Thermodynamical reconstruction of parent cell
    //Averaging conservative variables
    m_cons->setToZero(m_numberPhases);
    for (int i = 0; i < numberCellsChildren; i++) {
      m_cons->setBufferFlux(*m_childrenCells[i], m_numberPhases);              //Building children cells conservative variables in fluxTempXXX
      m_cons->addFlux(1., m_numberPhases);                                     //Adding them in m_cons
    }
    m_cons->multiply(1. / (double)numberCellsChildren, m_numberPhases);        //Division of m_cons by the children number to obatin average
    //Primitive variables reconstruction using relaxation
    m_cons->buildPrim(m_vecPhases, m_mixture, m_numberPhases);
	m_model->relaxations(this, m_numberPhases);

    //Transport
    for (int k = 0; k < m_numberTransports; k++) {
      double transport(0.);
      for (int i = 0; i < numberCellsChildren; i++) {
        transport += m_childrenCells[i]->getTransport(k).getValue();
      }
      transport = transport / (double)numberCellsChildren;
      m_vecTransports[k].setValue(transport);
    }

    //setting m_cons to zero for next
    m_cons->setToZero(m_numberPhases);
    for (int k = 0; k < m_numberTransports; k++) {
      m_consTransports[k].setValue(0.);
    }
  }
}

//***********************************************************************

bool Cell::lvlNeighborTooHigh()
{
  bool criteria(false);
  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    if (m_boundaries[b]->getLvl() == m_lvl) {
      for (int bChild = 0; bChild < m_boundaries[b]->getNumberBordsChildren(); bChild++) {
        if (m_boundaries[b]->getBordChild(bChild)->getSplit()) { criteria = true; return criteria; }
      }
    }
    else {
      if (m_boundaries[b]->getSplit()) { criteria = true; return criteria; }
    }
  }
  return criteria;
}

//***********************************************************************

bool Cell::lvlNeighborTooLow()
{
  bool criteria(false);
  for (unsigned int b = 0; b < m_boundaries.size(); b++) {
    if (!m_boundaries[b]->getSplit()) {
      if (m_boundaries[b]->whoAmI() == 0) //Boundary type CellInterface/O2
      {
        //Getting left and right AMR levels for each boundary
        int lvlg = m_boundaries[b]->getCellGauche()->getLvl();
        int lvld = m_boundaries[b]->getCellDroite()->getLvl();

        //Looking for neighbor level to be not too low
        if (lvlg < m_lvl) { criteria = true; return criteria; }
        if (lvld < m_lvl) { criteria = true; return criteria; }
      }
      else
      {
        //Getting Left AMR levels for boundary conditions
        int lvlg = m_boundaries[b]->getCellGauche()->getLvl();

        //Looking for neighbor level to be not too low
        if (lvlg < m_lvl) { criteria = true; return criteria; }
      }
    }
  }
  return criteria;
}

//***********************************************************************

void Cell::buildLvlCellsAndLvlInternalBoundariesArrays(vector<Cell *> *cellsLvl, vector<CellInterface *> *boundariesLvl)
{
  for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
    cellsLvl[m_lvl + 1].push_back(m_childrenCells[i]);
  }
  for (unsigned int i = 0; i < m_childrenInternalBoundaries.size(); i++) {
    boundariesLvl[m_lvl + 1].push_back(m_childrenInternalBoundaries[i]);
  }
}

//***********************************************************************

bool Cell::printGnuplotAMR(std::ofstream &fileStream, const int &dim, GeometricObject *objet)
{
  bool ecrit(true);
  int dimension(dim);
  Coord position = m_element->getPosition();
  //for cut printing
  if (objet != 0) 
  {
    if (objet->getType() != 0) { //For non probes objects
      ecrit = m_element->traverseObjet(*objet);
      position = objet->projectionPoint(position);
      dimension = objet->getType();
    }
  }
  if (ecrit) {
    if (!m_split) {
      //CAUTION: the order has to be kept for conformity reasons with gnuplot scripts
      if (dimension >= 1) fileStream << position.getX() << " ";
      if (dimension >= 2) fileStream << position.getY() << " ";
      if (dimension == 3)fileStream << position.getZ() << " ";
      this->printPhasesMixture(m_numberPhases, m_numberTransports, fileStream);
      fileStream << m_lvl << " " << m_xi << " ";
      fileStream << endl;
      if (objet != 0) { if (objet->getType() == 0) return true; } //probe specificity, unique.
    }
    else {
      
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        m_childrenCells[i]->printGnuplotAMR(fileStream, dim, objet);
      }
    }
  } 
  return false;
}

//***********************************************************************

void Cell::computeIntegration(double &integration)
{
  if (!m_split) {
    integration += m_element->getVolume()*m_vecPhases[1]->getAlpha();
  }
  else {
    for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
      m_childrenCells[i]->computeIntegration(integration);
    }
  }
}

//***********************************************************************

void Cell::lookForPmax(double *pMax, double*pMaxWall)
{
  if (!m_split) {
    if (m_mixture->getPressure() > pMax[0]) {
      pMax[0] = m_mixture->getPressure();
      pMax[1] = m_element->getPosition().getX();
      pMax[2] = m_element->getPosition().getY();
      pMax[3] = m_element->getPosition().getZ();
    }
    if (m_mixture->getPressure() > pMaxWall[0] && m_element->getPosition().getX() < 0.005) {
      pMaxWall[0] = m_mixture->getPressure();
      pMaxWall[1] = m_element->getPosition().getX();
      pMaxWall[2] = m_element->getPosition().getY();
      pMaxWall[3] = m_element->getPosition().getZ();
    }
  }
  else {
    for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
      m_childrenCells[i]->lookForPmax(pMax, pMaxWall);
    }
  }
}

//***********************************************************************

int Cell::getLvl()
{
  return m_lvl;
}

//***********************************************************************

bool Cell::getSplit()
{
	return m_split;
}

//***********************************************************************

double Cell::getXi()
{
  return m_xi;
}

//***********************************************************************

void Cell::setXi(double value)
{
  m_xi = value;
}

//***********************************************************************

void Cell::addFluxXi(double value)
{
  m_consXi += value;
}

//***********************************************************************

void Cell::subtractFluxXi(double value)
{
  m_consXi -= value;
}

//***********************************************************************

int Cell::getNumberCellsChildren()
{
  return m_childrenCells.size();
}

//***********************************************************************

Cell* Cell::getCellChild(const int &num)
{
  return m_childrenCells[num];
}

//***********************************************************************

std::vector<Cell *>* Cell::getChildVector()
{
  return &m_childrenCells;
}

//****************************************************************************
//************************** Parallel non-AMR *******************************
//****************************************************************************

void Cell::fillBufferPrimitives(double *buffer, int &counter, Prim type) const
{
	for (int k = 0; k < m_numberPhases; k++) {
		this->getPhase(k, type)->fillBuffer(buffer, counter);
	}
	this->getMixture(type)->fillBuffer(buffer, counter);
	for (int k = 0; k < m_numberTransports; k++) {
		buffer[++counter] = this->getTransport(k, type).getValue();
	}
}

//***********************************************************************

void Cell::getBufferPrimitives(double *buffer, int &counter, Eos **eos, Prim type)
{
	for (int k = 0; k < m_numberPhases; k++) {
		this->getPhase(k, type)->getBuffer(buffer, counter, eos);
	}
	this->getMixture(type)->getBuffer(buffer, counter);
	for (int k = 0; k < m_numberTransports; k++) {
		this->setTransport(buffer[++counter], k, type);
	}
  this->fulfillState(type);
}

//***********************************************************************

void Cell::fillBufferVector(double *buffer, int &counter, const int &dim, std::string nameVector, int num, int index) const
{
	buffer[++counter] = this->selectVector(nameVector, num, index).getX();
	if (dim > 1) buffer[++counter] = this->selectVector(nameVector, num, index).getY();
	if (dim > 2) buffer[++counter] = this->selectVector(nameVector, num, index).getZ();
}

//***********************************************************************

void Cell::getBufferVector(double *buffer, int &counter, const int &dim, std::string nameVector, int num, int index)
{
	Coord temp;
	temp.setX(buffer[++counter]);
	if (dim > 1) temp.setY(buffer[++counter]);
	if (dim > 2) temp.setZ(buffer[++counter]);
	this->setVector(nameVector, temp, num, index);
}

//***********************************************************************

void Cell::fillBufferTransports(double *buffer, int &counter) const
{
  for (int k = 0; k < m_numberTransports; k++) {
    buffer[++counter] = this->getTransport(k).getValue();
  }
}

//***********************************************************************

void Cell::getBufferTransports(double *buffer, int &counter)
{
  for (int k = 0; k < m_numberTransports; k++) {
    this->setTransport(buffer[++counter], k);
  }
}

//****************************************************************************
//**************************** AMR Parallel **********************************
//****************************************************************************

void Cell::fillBufferPrimitivesAMR(double *buffer, int &counter, const int &lvl, string whichCpuAmIForNeighbour, Prim type) const
{
  if (m_lvl == lvl) {
    for (int k = 0; k < m_numberPhases; k++) {
      this->getPhase(k, type)->fillBuffer(buffer, counter);
    }
    this->getMixture(type)->fillBuffer(buffer, counter);
    for (int k = 0; k < m_numberTransports; k++) {
      buffer[++counter] = this->getTransport(k, type).getValue();
    }
  }
  else {
    if (whichCpuAmIForNeighbour == "LEFT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am left CPU, I send all right cells
        if ((i % 2) == 1) { m_childrenCells[i]->fillBufferPrimitivesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, type); }
      }
    }
    else if (whichCpuAmIForNeighbour == "RIGHT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am right CPU, I send all left cells
        if ((i % 2) == 0) { m_childrenCells[i]->fillBufferPrimitivesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, type); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BOTTOM") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am bottom CPU, I send all top cells
        if ((i % 4) > 1) { m_childrenCells[i]->fillBufferPrimitivesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, type); }
      }
    }
    else if (whichCpuAmIForNeighbour == "TOP") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am top CPU, I send all bottom cells
        if ((i % 4) <= 1) { m_childrenCells[i]->fillBufferPrimitivesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, type); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BACK") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am back CPU, I send all front cells
        if (i > 3) { m_childrenCells[i]->fillBufferPrimitivesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, type); }
      }
    }
    else if (whichCpuAmIForNeighbour == "FRONT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am front CPU, I send all back cells
        if (i <= 3) { m_childrenCells[i]->fillBufferPrimitivesAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, type); }
      }
    }
  }
}

//***********************************************************************

void Cell::getBufferPrimitivesAMR(double *buffer, int &counter, const int &lvl, Eos **eos, Prim type)
{
	if (m_lvl == lvl) {
		for (int k = 0; k < m_numberPhases; k++) {
			this->getPhase(k, type)->getBuffer(buffer, counter, eos);
		}
		this->getMixture(type)->getBuffer(buffer, counter);
		for (int k = 0; k < m_numberTransports; k++) {
			this->setTransport(buffer[++counter], k, type);
		}
    this->fulfillState(type);
	}
	else {
		for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
			m_childrenCells[i]->getBufferPrimitivesAMR(buffer, counter, lvl, eos, type);
		}
	}
}

//***********************************************************************

void Cell::fillBufferVectorAMR(double *buffer, int &counter, const int &lvl, string whichCpuAmIForNeighbour, const int &dim, std::string nameVector, int num, int index) const
{
	if (m_lvl == lvl) {
		buffer[++counter] = this->selectVector(nameVector, num, index).getX();
		if (dim > 1) buffer[++counter] = this->selectVector(nameVector, num, index).getY();
		if (dim > 2) buffer[++counter] = this->selectVector(nameVector, num, index).getZ();
	}
	else {
    if (whichCpuAmIForNeighbour == "LEFT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am left CPU, I send all right cells
        if ((i % 2) == 1) { m_childrenCells[i]->fillBufferVectorAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, dim, nameVector, num, index); }
      }
    }
    else if (whichCpuAmIForNeighbour == "RIGHT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am right CPU, I send all left cells
        if ((i % 2) == 0) { m_childrenCells[i]->fillBufferVectorAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, dim, nameVector, num, index); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BOTTOM") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am bottom CPU, I send all top cells
        if ((i % 4) > 1) { m_childrenCells[i]->fillBufferVectorAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, dim, nameVector, num, index); }
      }
    }
    else if (whichCpuAmIForNeighbour == "TOP") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am top CPU, I send all bottom cells
        if ((i % 4) <= 1) { m_childrenCells[i]->fillBufferVectorAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, dim, nameVector, num, index); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BACK") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am back CPU, I send all front cells
        if (i > 3) { m_childrenCells[i]->fillBufferVectorAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, dim, nameVector, num, index); }
      }
    }
    else if (whichCpuAmIForNeighbour == "FRONT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am front CPU, I send all back cells
        if (i <= 3) { m_childrenCells[i]->fillBufferVectorAMR(buffer, counter, lvl, whichCpuAmIForNeighbour, dim, nameVector, num, index); }
      }
    }
	}
}

//***********************************************************************

void Cell::getBufferVectorAMR(double *buffer, int &counter, const int &lvl, const int &dim, std::string nameVector, int num, int index)
{
	if (m_lvl == lvl) {
		Coord temp;
		temp.setX(buffer[++counter]);
		if (dim > 1) temp.setY(buffer[++counter]);
		if (dim > 2) temp.setZ(buffer[++counter]);
		this->setVector(nameVector, temp, num, index);
	}
	else {
		for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
			m_childrenCells[i]->getBufferVectorAMR(buffer, counter, lvl, dim, nameVector, num, index);
		}
	}
}

//***********************************************************************

void Cell::fillBufferTransportsAMR(double *buffer, int &counter, const int &lvl, string whichCpuAmIForNeighbour) const
{
  if (m_lvl == lvl) {
    for (int k = 0; k < m_numberTransports; k++) {
      buffer[++counter] = this->getTransport(k).getValue();
    }
  }
  else {
    if (whichCpuAmIForNeighbour == "LEFT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am left CPU, I send all right cells
        if ((i % 2) == 1) { m_childrenCells[i]->fillBufferTransportsAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "RIGHT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am right CPU, I send all left cells
        if ((i % 2) == 0) { m_childrenCells[i]->fillBufferTransportsAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BOTTOM") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am bottom CPU, I send all top cells
        if ((i % 4) > 1) { m_childrenCells[i]->fillBufferTransportsAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "TOP") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am top CPU, I send all bottom cells
        if ((i % 4) <= 1) { m_childrenCells[i]->fillBufferTransportsAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BACK") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am back CPU, I send all front cells
        if (i > 3) { m_childrenCells[i]->fillBufferTransportsAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "FRONT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am front CPU, I send all back cells
        if (i <= 3) { m_childrenCells[i]->fillBufferTransportsAMR(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
  }
}

//***********************************************************************

void Cell::getBufferTransportsAMR(double *buffer, int &counter, const int &lvl)
{
  if (m_lvl == lvl) {
    for (int k = 0; k < m_numberTransports; k++) {
      this->setTransport(buffer[++counter], k);
    }
  }
  else {
    for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
      m_childrenCells[i]->getBufferTransportsAMR(buffer, counter, lvl);
    }
  }
}

//***********************************************************************

void Cell::chooseRefineDeraffineGhost(const int &nbCellsY, const int &nbCellsZ,	const vector<AddPhys*> &addPhys, Model *model, vector<Cell *> *cellsLvlGhost)
{
  if (m_split) {
    //cout << "cpu " << rankCpu << " split" << endl;
    if (m_childrenCells.size() == 0) { this->refineCellAndBoundariesGhost(nbCellsY, nbCellsZ, addPhys, model); }
  }
  else {
    //cout << "cpu " << rankCpu << " non-split" << endl;
    if (m_childrenCells.size() > 0) { this->unrefineCellAndBoundariesGhost(); }
  }
  //cout << "cpu " << rankCpu << " push-back" << endl;
  for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
    cellsLvlGhost[m_lvl + 1].push_back(m_childrenCells[i]);
  }
}

//***********************************************************************

void Cell::refineCellAndBoundariesGhost(const int &nbCellsY, const int &nbCellsZ, const vector<AddPhys*> &addPhys, Model *model)
{
	//--------------------------------------
	//Initializations (children and dimension)
	//--------------------------------------

  //Notice that the children number of ghost cells is different than for internal cells
	double dimX(1.), dimY(0.), dimZ(0.);
	int numberCellsChildren(1);
	int dim(1);
	if (nbCellsZ != 1) {
		numberCellsChildren = 4;
		dimY = 1.;
		dimZ = 1.;
		dim = 3;
	}
	else if (nbCellsY != 1) {
		numberCellsChildren = 2;
		dimY = 1.;
		dim = 2;
	}
	CellInterface* bordRef(0); //bordRef always is the one with the level corresponding to the current ghost-cell level
	for (unsigned int b = 0; b < m_boundaries.size(); b++) {
		if (m_boundaries[b]->whoAmI() == 0) { bordRef = m_boundaries[b]; break; } //Boundary type CellInterface/O2
	}
	int allocateSlopeLocal = 1;

	//---------------
	//Cell refinement
	//---------------

	//Initialization of mesh data for children cells
	//----------------------------------------------
	double posXCellParent, posYCellParent, posZCellParent;
	double dXParent, dYParent, dZParent;
	double posXChild, posYChild, posZChild;
	posXCellParent = m_element->getPosition().getX();
	posYCellParent = m_element->getPosition().getY();
	posZCellParent = m_element->getPosition().getZ();
  dXParent = m_element->getSizeX();
  dYParent = m_element->getSizeY();
  dZParent = m_element->getSizeZ();
	double volumeCellParent, lCFLCellParent;
	volumeCellParent = m_element->getVolume();
	lCFLCellParent = m_element->getLCFL();

	for (int i = 0; i < numberCellsChildren; i++) {

		//Children cells creation
		//-----------------------
		this->createChildCell(i, m_lvl);
		m_element->creerElementChild();
		m_childrenCells[i]->setElement(m_element->getElementChild(i), i);
		m_childrenCells[i]->getElement()->setVolume(volumeCellParent / (double)numberCellsChildren / 2);
		m_childrenCells[i]->getElement()->setLCFL(0.5*lCFLCellParent);
    m_childrenCells[i]->getElement()->setSize((1 - dimX*0.5)*dXParent, (1 - dimY*0.5)*dYParent, (1 - dimZ*0.5)*dZParent);
    //Face in the x-direction
    if (abs(bordRef->getFace()->getNormal().getX()) > 0.99) {
      //Ghost cells are on the right according to internal cells, the positions of the children are thus those on the left side
      if (bordRef->getFace()->getPos().getX() < m_element->getPosition().getX()) { posXChild = posXCellParent - dimX*dXParent*0.25; }
      //Ghost cells are on the left according to internal cells, the positions of the children are thus those on the right side
      else { posXChild = posXCellParent + dimX*dXParent*0.25; }
      posYChild = posYCellParent + dimY*dYParent*(double)(-0.25 + 0.5 * (i % 2));
      posZChild = posZCellParent + dimZ*dZParent*(double)(-0.25 + 0.5 * ((i / 2) % 2));
    }
    //Face in the y-direction
    else if (abs(bordRef->getFace()->getNormal().getY()) > 0.99) {
      //Ghost cells are on the top according to internal cells, the positions of the children are thus those on the bottom side
      if (bordRef->getFace()->getPos().getY() < m_element->getPosition().getY()) { posYChild = posYCellParent - dimY*dYParent*0.25; }
      //Ghost cells are on the bottom according to internal cells, the positions of the children are thus those on the top side
      else { posYChild = posYCellParent + dimY*dYParent*0.25; }
      posXChild = posXCellParent + dimX*dXParent*(double)(-0.25 + 0.5 * (i % 2));
      posZChild = posZCellParent + dimZ*dZParent*(double)(-0.25 + 0.5 * ((i / 2) % 2));
    }
    //Face in the z-direction
    else if (abs(bordRef->getFace()->getNormal().getZ()) > 0.99) {
      //Ghost cells are on the front according to internal cells, the positions of the children are thus those on the back side
      if (bordRef->getFace()->getPos().getZ() < m_element->getPosition().getZ()) { posZChild = posZCellParent - dimZ*dZParent*0.25; }
      //Ghost cells are on the back according to internal cells, the positions of the children are thus those on the front side
      else { posZChild = posZCellParent + dimZ*dZParent*0.25; }
      posXChild = posXCellParent + dimX*dXParent*(double)(-0.25 + 0.5 * (i % 2));
      posYChild = posYCellParent + dimY*dYParent*(double)(-0.25 + 0.5 * ((i / 2) % 2));
    }
		m_childrenCells[i]->getElement()->setPos(posXChild, posYChild, posZChild);

    //Initialization of main arrays according to model and number of phases
    //+ physical initialization: physical data for children cells
		//-------------------------------------------------------------------------------------
		m_childrenCells[i]->allocate(m_numberPhases, m_numberTransports, addPhys, model);
		for (int k = 0; k < m_numberPhases; k++) {
			m_childrenCells[i]->copyPhase(k, m_vecPhases[k]);
		}
		m_childrenCells[i]->copyMixture(m_mixture);
		m_childrenCells[i]->getCons()->setToZero(m_numberPhases);
		for (int k = 0; k < m_numberTransports; k++) { m_childrenCells[i]->setTransport(m_vecTransports[k].getValue(), k); }
		for (int k = 0; k < m_numberTransports; k++) { m_childrenCells[i]->setConsTransport(0., k); }
		m_childrenCells[i]->setXi(m_xi);
	}

	//------------------------------
	//External boundaries refinement
	//------------------------------

	for (unsigned int b = 0; b < m_boundaries.size(); b++) {
		if (!m_boundaries[b]->getSplit()) { m_boundaries[b]->raffineBordExterneGhost(nbCellsY, nbCellsZ, dXParent, dYParent, dZParent, this, dim); }
	}
}

//***********************************************************************

void Cell::unrefineCellAndBoundariesGhost()
{
	//----------------------------------------
	//External children boundaries destruction
	//----------------------------------------

	for (unsigned int b = 0; b < m_boundaries.size(); b++) {
		m_boundaries[b]->deraffineBordExterne(this);
	}

	//--------------------------
	//Children cells destruction
	//--------------------------

	for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
		delete m_childrenCells[i];
	}
	m_childrenCells.clear();
	m_element->finalizeElementsChildren();
}

//***********************************************************************

void Cell::fillBufferXi(double *buffer, int &counter, const int &lvl, string whichCpuAmIForNeighbour) const
{
	if (m_lvl == lvl) {
		buffer[++counter] = m_xi;
	}
	else {
    if (whichCpuAmIForNeighbour == "LEFT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am left CPU, I send all right cells
        if ((i % 2) == 1) { m_childrenCells[i]->fillBufferXi(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "RIGHT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am right CPU, I send all left cells
        if ((i % 2) == 0) { m_childrenCells[i]->fillBufferXi(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BOTTOM") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am bottom CPU, I send all top cells
        if ((i % 4) > 1) { m_childrenCells[i]->fillBufferXi(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "TOP") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am top CPU, I send all bottom cells
        if ((i % 4) <= 1) { m_childrenCells[i]->fillBufferXi(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BACK") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am back CPU, I send all front cells
        if (i > 3) { m_childrenCells[i]->fillBufferXi(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "FRONT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am front CPU, I send all back cells
        if (i <= 3) { m_childrenCells[i]->fillBufferXi(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
	}
}

//***********************************************************************

void Cell::getBufferXi(double *buffer, int &counter, const int &lvl)
{
	if (m_lvl == lvl) {
		m_xi = buffer[++counter];
	}
	else {
		for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
			m_childrenCells[i]->getBufferXi(buffer, counter, lvl);
		}
	}
}

//***********************************************************************

void Cell::fillBufferSplit(bool *buffer, int &counter, const int &lvl, string whichCpuAmIForNeighbour) const
{
	if (m_lvl == lvl) {
		buffer[++counter] = m_split;
	}
	else {
    if (whichCpuAmIForNeighbour == "LEFT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am left CPU, I send all right cells
        if ((i % 2) == 1) { m_childrenCells[i]->fillBufferSplit(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "RIGHT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am right CPU, I send all left cells
        if ((i % 2) == 0) { m_childrenCells[i]->fillBufferSplit(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BOTTOM") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am bottom CPU, I send all top cells
        if ((i % 4) > 1) { m_childrenCells[i]->fillBufferSplit(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "TOP") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am top CPU, I send all bottom cells
        if ((i % 4) <= 1) { m_childrenCells[i]->fillBufferSplit(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BACK") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am back CPU, I send all front cells
        if (i > 3) { m_childrenCells[i]->fillBufferSplit(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "FRONT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am front CPU, I send all back cells
        if (i <= 3) { m_childrenCells[i]->fillBufferSplit(buffer, counter, lvl, whichCpuAmIForNeighbour); }
      }
    }
	}
}

//***********************************************************************

void Cell::getBufferSplit(bool *buffer, int &counter, const int &lvl)
{
	if (m_lvl == lvl) {
		m_split = buffer[++counter];
	}
	else {
		for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
			m_childrenCells[i]->getBufferSplit(buffer, counter, lvl);
		}
	}
}

//***********************************************************************

void Cell::fillNumberElementsToSendToNeighbour(int &numberElementsToSendToNeighbor, const int &lvl, string whichCpuAmIForNeighbour)
{
	if (m_lvl == lvl) {
		numberElementsToSendToNeighbor++;
	}
	else {
    if (whichCpuAmIForNeighbour == "LEFT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am left CPU, I send all right cells
        if ((i % 2) == 1) { m_childrenCells[i]->fillNumberElementsToSendToNeighbour(numberElementsToSendToNeighbor, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "RIGHT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am right CPU, I send all left cells
        if ((i % 2) == 0) { m_childrenCells[i]->fillNumberElementsToSendToNeighbour(numberElementsToSendToNeighbor, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BOTTOM") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am bottom CPU, I send all top cells
        if ((i % 4) > 1) { m_childrenCells[i]->fillNumberElementsToSendToNeighbour(numberElementsToSendToNeighbor, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "TOP") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am top CPU, I send all bottom cells
        if ((i % 4) <= 1) { m_childrenCells[i]->fillNumberElementsToSendToNeighbour(numberElementsToSendToNeighbor, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "BACK") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am back CPU, I send all front cells
        if (i > 3) { m_childrenCells[i]->fillNumberElementsToSendToNeighbour(numberElementsToSendToNeighbor, lvl, whichCpuAmIForNeighbour); }
      }
    }
    else if (whichCpuAmIForNeighbour == "FRONT") {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        //I am front CPU, I send all back cells
        if (i <= 3) { m_childrenCells[i]->fillNumberElementsToSendToNeighbour(numberElementsToSendToNeighbor, lvl, whichCpuAmIForNeighbour); }
      }
    }
	}
}

//***************************************************************************