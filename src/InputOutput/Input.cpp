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

//! \file      Input.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      July 20 2018

#include "Input.h"
#include "HeaderInputOutput.h"
#include <iostream>
#include "../Meshes/stretchZone.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

Input::Input(Run *run) : m_run(run)
{
	//Attribution des numbers de version pour la lecture
	m_vMain = 5;
	m_vMesh = 5;
	m_vCI = 4;
	m_vModel = 4;

	m_nameMain = "mainV" + IO::toString(m_vMain) + ".xml";
	m_nameMesh = "meshV" + IO::toString(m_vMesh) + ".xml";
	m_nameCI = "initialConditionsV" + IO::toString(m_vCI) + ".xml";
	m_nameModel = "modelV" + IO::toString(m_vModel) + ".xml";
}

//***********************************************************************

Input::~Input(){}

//***********************************************************************

void Input::lectureInputXML(vector<GeometricalDomain*> &domains, vector<BoundCond*> &boundCond)
{
  try{
    //1) Parametres generaux du compute
    entreeMain(m_run->m_simulationName);
    //2) Donnees de Mesh
    entreeMesh(m_run->m_simulationName);
    //3) Donnees models et fluides
    entreeModel(m_run->m_simulationName);
    //4) Lecture des conditions initiales
    entreeConditionsInitiales(m_run->m_simulationName, domains, boundCond);
  }
  catch (ErrorXML &e){
    if(rankCpu==0) cerr << e.infoError() << endl;
    throw;
  }
}

//***********************************************************************

void Input::entreeMain(string casTest)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
	  stringstream fileName(casTest + m_nameMain);
    XMLDocument xmlMain;
    XMLError error(xmlMain.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(),__FILE__, __LINE__);
    
    //2) Recuperation des donnees principales du compute
    //-------------------------------------------------
    //Recuperation racine du document XML
    XMLNode *computationParam = xmlMain.FirstChildElement("computationParam");
    if (computationParam == NULL) throw ErrorXMLRacine("computationParam", fileName.str(), __FILE__, __LINE__);

    XMLElement *element, *sousElement;

    //Recuperation name du run
    element = computationParam->FirstChildElement("run");
    if (element == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);
    XMLNode* xmlNode2 = element->FirstChild();
    if (xmlNode2 == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);
    XMLText* xmlText = xmlNode2->ToText();
    if (xmlText == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);

    //Lecture des informations de sorties/prints
    element = computationParam->FirstChildElement("outputMode");
    if (element == NULL) throw ErrorXMLElement("outputMode", fileName.str(), __FILE__, __LINE__);
    //Lecture format sortie
    string format(element->Attribute("format"));
    if (format == "") throw ErrorXMLAttribut("format", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(format);
    if (format == "XML") { m_run->m_outPut = new OutputXML(casTest, xmlText->Value(), element, fileName.str(), this); }
    else if (format == "GNU") { m_run->m_outPut = new OutputGNU(casTest, xmlText->Value(), element, fileName.str(), this); }
    else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }

    //Lecture des cuts 1D
    element = computationParam->FirstChildElement("cut1D");
    while (element != NULL)
    {
      m_run->m_cuts.push_back(new OutputCutGNU(casTest, xmlText->Value(), element, fileName.str(), LINE, this));
      element = element->NextSiblingElement("cut1D");
    }
    //Lecture des cuts 2D
    element = computationParam->FirstChildElement("cut2D");
    while (element != NULL)
    {
      m_run->m_cuts.push_back(new OutputCutGNU(casTest, xmlText->Value(), element, fileName.str(), PLAN, this));
      element = element->NextSiblingElement("cut2D");
    }
    //Reading probes
    element = computationParam->FirstChildElement("probe");
    while (element != NULL)
    {
      m_run->m_probes.push_back(new OutputProbeGNU(casTest, xmlText->Value(), element, fileName.str(), this));
      element = element->NextSiblingElement("probe");
    }
    
    //Recuperation Iteration / temps Physique
    element = computationParam->FirstChildElement("timeControlMode");
    if (element == NULL) throw ErrorXMLElement("timeControlMode", fileName.str(), __FILE__, __LINE__);
    error = element->QueryBoolAttribute("iterations", &m_run->m_controleIterations);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("iterations", fileName.str(), __FILE__, __LINE__);
    if (m_run->m_controleIterations)
    {
      //Recuperation Iterations / Frequence
      sousElement = element->FirstChildElement("iterations");
      if (sousElement == NULL) throw ErrorXMLElement("iterations", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryIntAttribute("number", &m_run->m_nbIte);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryIntAttribute("iterFreq", &m_run->m_freq);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("iterFreq", fileName.str(), __FILE__, __LINE__);
    }
    else
    {
      //Recuperation Temps / Frequence
      sousElement = element->FirstChildElement("physicalTime");
      if (sousElement == NULL) throw ErrorXMLElement("physicalTime", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryFloatAttribute("totalTime", &m_run->m_finalPhysicalTime);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("totalTime", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryFloatAttribute("timeFreq", &m_run->m_timeFreq);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("timeFreq", fileName.str(), __FILE__, __LINE__);
    }

    //Recuperation CFL
    element = computationParam->FirstChildElement("computationControl");
    if (element == NULL) throw ErrorXMLElement("computationControl", fileName.str(), __FILE__, __LINE__);
    error = element->QueryDoubleAttribute("CFL", &m_run->m_cfl);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("CFL", fileName.str(), __FILE__, __LINE__);

    //Lecture Ordre2
    element = computationParam->FirstChildElement("secondOrder");
    if(element == NULL) {
      m_run->m_order = "FIRSTORDER";
      m_run->m_globalLimiter = new Limiter;
      m_run->m_interfaceLimiter = new Limiter;
      m_run->m_globalVolumeFractionLimiter = new Limiter;
      m_run->m_interfaceVolumeFractionLimiter = new Limiter;
    }
    else{
      m_run->m_order = "SECONDORDER";
      XMLNode *contenu;
      //Recuperation global limiter
      sousElement = element->FirstChildElement("globalLimiter");
      if (sousElement == NULL) throw ErrorXMLElement("globalLimiter", fileName.str(), __FILE__, __LINE__);
      contenu = sousElement->FirstChild();
      if (contenu == NULL) throw ErrorXMLElement("globalLimiter", fileName.str(), __FILE__, __LINE__);
      string globalLimiter = contenu->ToText()->Value();
      Tools::uppercase(globalLimiter);
      if (globalLimiter == "MINMOD") { m_run->m_globalLimiter = new LimiterMinmod; }
      else if (globalLimiter == "VANLEER") { m_run->m_globalLimiter = new LimiterVanLeer; }
      else if (globalLimiter == "VANALBADA") { m_run->m_globalLimiter = new LimiterVanAlbada; }
      else if (globalLimiter == "SUPERBEE") { m_run->m_globalLimiter = new LimiterSuperBee; }
      else if (globalLimiter == "MC") { m_run->m_globalLimiter = new LimiterMC; }
      else if (globalLimiter == "THINC") { throw ErrorXMLAttribut("THINC can only be a volume-fraction limiter", fileName.str(), __FILE__, __LINE__); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      //Recuperation interface limiter
      string interfaceLimiter = globalLimiter;
      sousElement = element->FirstChildElement("interfaceLimiter");
      if (sousElement != NULL) {
        contenu = sousElement->FirstChild();
        if (contenu == NULL) throw ErrorXMLElement("interfaceLimiter", fileName.str(), __FILE__, __LINE__);
        interfaceLimiter = contenu->ToText()->Value();
      }
      Tools::uppercase(interfaceLimiter);
      if (interfaceLimiter == "MINMOD") { m_run->m_interfaceLimiter = new LimiterMinmod; }
      else if (interfaceLimiter == "VANLEER") { m_run->m_interfaceLimiter = new LimiterVanLeer; }
      else if (interfaceLimiter == "VANALBADA") { m_run->m_interfaceLimiter = new LimiterVanAlbada; }
      else if (interfaceLimiter == "SUPERBEE") { m_run->m_interfaceLimiter = new LimiterSuperBee; }
      else if (interfaceLimiter == "MC") { m_run->m_interfaceLimiter = new LimiterMC; }
      else if (interfaceLimiter == "THINC") { throw ErrorXMLAttribut("THINC can only be a volume-fraction limiter", fileName.str(), __FILE__, __LINE__); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      //Recuperation global volume-fraction limiter
      string globalVolumeFractionLimiter = globalLimiter;
      sousElement = element->FirstChildElement("globalVolumeFractionLimiter");
      if (sousElement != NULL) {
        contenu = sousElement->FirstChild();
        if (contenu == NULL) throw ErrorXMLElement("globalVolumeFractionLimiter", fileName.str(), __FILE__, __LINE__);
        globalVolumeFractionLimiter = contenu->ToText()->Value();
      }
      Tools::uppercase(globalVolumeFractionLimiter);
      if (globalVolumeFractionLimiter == "MINMOD") { m_run->m_globalVolumeFractionLimiter = new LimiterMinmod; }
      else if (globalVolumeFractionLimiter == "VANLEER") { m_run->m_globalVolumeFractionLimiter = new LimiterVanLeer; }
      else if (globalVolumeFractionLimiter == "VANALBADA") { m_run->m_globalVolumeFractionLimiter = new LimiterVanAlbada; }
      else if (globalVolumeFractionLimiter == "SUPERBEE") { m_run->m_globalVolumeFractionLimiter = new LimiterSuperBee; }
      else if (globalVolumeFractionLimiter == "MC") { m_run->m_globalVolumeFractionLimiter = new LimiterMC; }
      else if (globalVolumeFractionLimiter == "THINC") { m_run->m_globalVolumeFractionLimiter = new LimiterTHINC; }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      //Recuperation interface volume-fraction limiter
      string interfaceVolumeFractionLimiter = interfaceLimiter;
      sousElement = element->FirstChildElement("interfaceVolumeFractionLimiter");
      if (sousElement != NULL) {
        contenu = sousElement->FirstChild();
        if (contenu == NULL) throw ErrorXMLElement("interfaceVolumeFractionLimiter", fileName.str(), __FILE__, __LINE__);
        interfaceVolumeFractionLimiter = contenu->ToText()->Value();
      }
      Tools::uppercase(interfaceVolumeFractionLimiter);
      if (interfaceVolumeFractionLimiter == "MINMOD") { m_run->m_interfaceVolumeFractionLimiter = new LimiterMinmod; }
      else if (interfaceVolumeFractionLimiter == "VANLEER") { m_run->m_interfaceVolumeFractionLimiter = new LimiterVanLeer; }
      else if (interfaceVolumeFractionLimiter == "VANALBADA") { m_run->m_interfaceVolumeFractionLimiter = new LimiterVanAlbada; }
      else if (interfaceVolumeFractionLimiter == "SUPERBEE") { m_run->m_interfaceVolumeFractionLimiter = new LimiterSuperBee; }
      else if (interfaceVolumeFractionLimiter == "MC") { m_run->m_interfaceVolumeFractionLimiter = new LimiterMC; }
      else if (interfaceVolumeFractionLimiter == "THINC") { m_run->m_interfaceVolumeFractionLimiter = new LimiterTHINC; }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }

    }

    //Reprise de Calcul depuis file resultat
    element = computationParam->FirstChildElement("resumeSimulation");
    if (element != NULL) { 
      error = element->QueryIntAttribute("fileNumber", &m_run->m_resumeSimulation);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName.str(), __FILE__, __LINE__);
    }

  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Input::entreeMesh(string casTest)
{
  try{
    //Methode AMR: initialization des variables
    m_run->m_lvlMax = 0;
		double criteriaVar(1.e10);
		bool varRho(false), varP(false), varU(false), varAlpha(false);
		double xiSplit(1.), xiJoin(1.);

    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    stringstream fileName(casTest + m_nameMesh);
    XMLDocument xmlMesh;
    XMLError error(xmlMesh.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees principales du compute
    //-------------------------------------------------
    //Recuperation racine du document XML
    XMLNode *mesh = xmlMesh.FirstChildElement("mesh");
    if (mesh == NULL) throw ErrorXMLRacine("mesh", fileName.str(), __FILE__, __LINE__);

    XMLElement *element;
    //Recuperation du type de mesh
    element = mesh->FirstChildElement("type");
    if (element == NULL) throw ErrorXMLElement("type", fileName.str(), __FILE__, __LINE__);
    string structureMesh(element->Attribute("structure"));
    if (structureMesh == "") throw ErrorXMLAttribut("structure", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(structureMesh);
    if (structureMesh == "UNSTRUCTURED")
    {
      //----------------MESH NON STRUCTURE ---------------------
      XMLElement *meshNS;
      meshNS = mesh->FirstChildElement("unstructuredMesh");
      if (meshNS == NULL) throw ErrorXMLElement("unstructuredMesh", fileName.str(), __FILE__, __LINE__);
      //Lecture name du file contenant les informations de mesh (.msh)
      element = meshNS->FirstChildElement("file");
      if (element == NULL) throw ErrorXMLElement("file", fileName.str(), __FILE__, __LINE__);
      string fichierMesh(element->Attribute("name"));
      if (fichierMesh == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      m_run->m_mesh = new MeshUnStruct(fichierMesh);
      //Recuperation pretraitement parallele
      element = meshNS->FirstChildElement("parallel");
      if (element != NULL) {
        error = element->QueryBoolAttribute("GMSHPretraitement", &m_run->m_parallelPreTreatment);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("GMSHPretraitement", fileName.str(), __FILE__, __LINE__);
      }
      //Methode AMR non possible avec mesh non structure
      element = meshNS->FirstChildElement("AMR");
      if (element != NULL) { throw ErrorXMLAttribut("Methode AMR non possible avec mesh non structure", fileName.str(), __FILE__, __LINE__); }
    }
    else if (structureMesh == "CARTESIAN")
    {
      //----------------MESH CARTESIAN ---------------------
      XMLElement *cartesianMesh;
      cartesianMesh = mesh->FirstChildElement("cartesianMesh");
      if (cartesianMesh == NULL) throw ErrorXMLElement("cartesianMesh", fileName.str(), __FILE__, __LINE__);
      //Recuperation des dimensions
      double lX, lY, lZ;
      element = cartesianMesh->FirstChildElement("dimensions");
      if (element == NULL) throw ErrorXMLElement("dimensions", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("x", &lX);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("y", &lY);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("z", &lZ);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName.str(), __FILE__, __LINE__);
      //Recuperation des numbers de mailles
      int nbX, nbY, nbZ;
      element = cartesianMesh->FirstChildElement("numberCells");
      if (element == NULL) throw ErrorXMLElement("numberCells", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("x", &nbX);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName.str(), __FILE__, __LINE__);
      if (nbX <= 1)  throw ErrorXMLAttribut("Number of cells in the x-direction has to be superior to 1", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("y", &nbY);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("z", &nbZ);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName.str(), __FILE__, __LINE__);
      if (nbZ > 1 && nbY <= 1) throw ErrorXMLAttribut("Number of cells in the z-direction can't be superior to 1 if number of cells in the y-direction is equal to 1", fileName.str(), __FILE__, __LINE__);
      
      //Stretching data
      vector<stretchZone> stretchX, stretchY, stretchZ;
      element = cartesianMesh->FirstChildElement("meshStretching");
      if (element != NULL) {
        double tempStart, tempEnd, tempFactor;
        int tempNumberCells;
        XMLElement *sousElement;
        //X stretching
        if (nbX > 1) {
          sousElement = element->FirstChildElement("XStretching");
          if (sousElement != NULL) {
            XMLElement *stretchElement(sousElement->FirstChildElement("stretch"));
            while (stretchElement != NULL) {
              error = stretchElement->QueryDoubleAttribute("startAt", &tempStart);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("startAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("endAt", &tempEnd);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("endAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("factor", &tempFactor);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("factor", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryIntAttribute("numberCells", &tempNumberCells);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberCells", fileName.str(), __FILE__, __LINE__);
              stretchX.push_back(stretchZone(tempStart, tempEnd, tempFactor, tempNumberCells));
              stretchElement = stretchElement->NextSiblingElement("stretch");
            }
          }
        }
        //Y stretching
        if (nbY > 1) {
          sousElement = element->FirstChildElement("YStretching");
          if (sousElement != NULL) {
            XMLElement *stretchElement(sousElement->FirstChildElement("stretch"));
            while (stretchElement != NULL) {
              error = stretchElement->QueryDoubleAttribute("startAt", &tempStart);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("startAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("endAt", &tempEnd);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("endAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("factor", &tempFactor);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("factor", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryIntAttribute("numberCells", &tempNumberCells);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberCells", fileName.str(), __FILE__, __LINE__);
              stretchY.push_back(stretchZone(tempStart, tempEnd, tempFactor, tempNumberCells));
              stretchElement = stretchElement->NextSiblingElement("stretch");
            }
          }
        }
        //Z stretching
        if (nbZ > 1) {
          sousElement = element->FirstChildElement("ZStretching");
          if (sousElement != NULL) {
            XMLElement *stretchElement(sousElement->FirstChildElement("stretch"));
            while (stretchElement != NULL) {
              error = stretchElement->QueryDoubleAttribute("startAt", &tempStart);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("startAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("endAt", &tempEnd);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("endAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("factor", &tempFactor);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("factor", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryIntAttribute("numberCells", &tempNumberCells);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberCells", fileName.str(), __FILE__, __LINE__);
              stretchZ.push_back(stretchZone(tempStart, tempEnd, tempFactor, tempNumberCells));
              stretchElement = stretchElement->NextSiblingElement("stretch");
            }
          }
        }

        //Stretching errors verifications
        stretchZone::verifyStretching(stretchX, lX, fileName.str());
        stretchZone::verifyStretching(stretchY, lY, fileName.str());
        stretchZone::verifyStretching(stretchZ, lZ, fileName.str());

      } //End stretching

      //Recuperation des variables pour methode AMR
      element = cartesianMesh->FirstChildElement("AMR");
      if (element != NULL) {
        //Interdiction parallele
        error = element->QueryIntAttribute("lvlMax", &m_run->m_lvlMax);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lvlMax", fileName.str(), __FILE__, __LINE__);
        error = element->QueryDoubleAttribute("criteriaVar", &criteriaVar);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("criteriaVar", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varRho", &varRho);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varRho", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varP", &varP);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varP", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varU", &varU);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varU", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varAlpha", &varAlpha);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varAlpha", fileName.str(), __FILE__, __LINE__);
        error = element->QueryDoubleAttribute("xiSplit", &xiSplit);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("xiSplit", fileName.str(), __FILE__, __LINE__);
        error = element->QueryDoubleAttribute("xiJoin", &xiJoin);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("xiJoin", fileName.str(), __FILE__, __LINE__);
        m_run->m_mesh = new MeshCartesianAMR(lX, nbX, lY, nbY, lZ, nbZ, stretchX, stretchY, stretchZ, m_run->m_lvlMax, criteriaVar, varRho, varP, varU, varAlpha, xiSplit, xiJoin);
      }
      else {
        m_run->m_mesh = new MeshCartesian(lX, nbX, lY, nbY, lZ, nbZ, stretchX, stretchY, stretchZ);
      }

    }
    else{ throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Input::entreeModel(string casTest)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    stringstream fileName(casTest + m_nameModel);
    XMLDocument xmlModel;
    XMLError error(xmlModel.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees du Model
    //-------------------------------------
    //Recuperation racine du document XML
    XMLNode *xmlNode = xmlModel.FirstChildElement("model");
    if (xmlNode == NULL) throw ErrorXMLRacine("model", fileName.str(), __FILE__, __LINE__);
    //Recuperation name du model resolu
    XMLElement *element(xmlNode->FirstChildElement("flowModel"));
    if (element == NULL) throw ErrorXMLElement("flowModel", fileName.str(), __FILE__, __LINE__);
    string model(element->Attribute("name"));
    if (model == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
    m_run->m_numberTransports = 0;
    error = element->QueryIntAttribute("numberTransports", &m_run->m_numberTransports);
    //if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberTransports", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(model);
    //Switch selon model
    bool alphaNull(false);
    if (model == "EULER"){ m_run->m_model = new ModEuler(m_run->m_numberTransports); m_run->m_numberPhases = 1; }
    else if (model == "EULERHOMOGENEOUS") { 
      int liquid, vapor;
      error = element->QueryIntAttribute("liquid", &liquid);
      error = element->QueryIntAttribute("vapor", &vapor);
      m_run->m_model = new ModEulerHomogeneous(m_run->m_numberTransports, liquid, vapor); m_run->m_numberPhases = 2; 
    }
    else if (model == "KAPILA")
    { 
      error = element->QueryIntAttribute("numberPhases", &m_run->m_numberPhases);
      m_run->m_model = new ModKapila(m_run->m_numberTransports, m_run->m_numberPhases);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberPhases", fileName.str(), __FILE__, __LINE__);
      error = element->QueryBoolAttribute("alphaNull", &alphaNull);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alphaNull", fileName.str(), __FILE__, __LINE__);
    }
    else if (model == "MULTIP")
    {
      error = element->QueryIntAttribute("numberPhases", &m_run->m_numberPhases);
      m_run->m_model = new ModMultiP(m_run->m_numberTransports, m_run->m_numberPhases);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberPhases", fileName.str(), __FILE__, __LINE__);
    }
    else if (model == "THERMALEQ")
    {
      error = element->QueryIntAttribute("numberPhases", &m_run->m_numberPhases);
      m_run->m_model = new ModThermalEq(m_run->m_numberTransports, m_run->m_numberPhases);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberPhases", fileName.str(), __FILE__, __LINE__);
    }
    else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
    if (m_run->m_globalVolumeFractionLimiter->AmITHINC() && (m_run->m_numberPhases != 2)) { throw ErrorXMLAttribut("Limiter THINC only works for 2 phases", fileName.str(), __FILE__, __LINE__); }
    if (m_run->m_interfaceVolumeFractionLimiter->AmITHINC() && (m_run->m_numberPhases != 2)) { throw ErrorXMLAttribut("Limiter THINC only works for 2 phases", fileName.str(), __FILE__, __LINE__); }
    
    //Thermodynamique
    vector<string> nameEOS;
    int EOSTrouvee(0);
    XMLElement *sousElement(xmlNode->FirstChildElement("EOS"));
    while (sousElement != NULL)
    {
      //Lecture Eos
      nameEOS.push_back(sousElement->Attribute("name"));
      if (nameEOS[EOSTrouvee] == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      EOSTrouvee++;
      sousElement = sousElement->NextSiblingElement("EOS");
    }
    m_run->m_numberEos = nameEOS.size();
    if (m_run->m_numberEos == 0) throw ErrorXMLEOS(fileName.str(), __FILE__, __LINE__);
    m_run->m_eos = new Eos*[m_run->m_numberPhases];
    int numberEos(0);
    for (int i = 0; i < m_run->m_numberEos; i++){ m_run->m_eos[i] = entreeEOS(nameEOS[i], numberEos); }
    m_run->m_eos[0]->assignEpsilonForAlphaNull(alphaNull, fileName.str()); //initialization of epsilon value for alphaNull option

    //Transports
    int transportsTrouvee(0);
    XMLElement *elementTransport(xmlNode->FirstChildElement("transport"));
    while (elementTransport != NULL)
    {
      //Lecture name transport
      m_run->m_nameGTR.push_back(elementTransport->Attribute("name"));
      if (m_run->m_nameGTR[transportsTrouvee] == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      transportsTrouvee++;
      elementTransport = elementTransport->NextSiblingElement("transport");
    }
    //Verifications
    if(transportsTrouvee < m_run->m_numberTransports) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //Reading of the additional physics
    int numberGPA(0);
    int physiqueAddTrouvee(0);
    element = xmlNode->FirstChildElement("additionalPhysic");
    while (element != NULL) {
      physiqueAddTrouvee++;
      //Lecture physique additionnelle
      string typeAddPhys(element->Attribute("type"));
      if (typeAddPhys == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeAddPhys);
      //switch sur le type de physique additionelle
      if (typeAddPhys == "SURFACETENSION") { 
        //switch model d ecoulement
        if (model == "KAPILA") { m_run->m_addPhys.push_back(new APKSurfaceTension(element, numberGPA, m_run->m_nameGTR, nameEOS, fileName.str())); }
        else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      }
      else if (typeAddPhys == "VISCOSITY") {
        //switch model d ecoulement
        if (model == "KAPILA") { m_run->m_addPhys.push_back(new APKViscosity(numberGPA, m_run->m_eos, m_run->m_numberPhases, fileName.str())); }
        else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      }
      else if (typeAddPhys == "CONDUCTIVITY") {
        //switch model d ecoulement
        if (model == "KAPILA") { m_run->m_addPhys.push_back(new APKConductivity(numberGPA, m_run->m_eos, m_run->m_numberPhases, fileName.str())); }
        else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("additionalPhysic");
    }
    m_run->m_numberAddPhys = physiqueAddTrouvee;

    //Symmetry terms reading
    m_run->m_symmetryAddPhys = new Symmetry();
    element = xmlNode->FirstChildElement("symmetryTerm");
    if (element == NULL) {
      m_run->m_symmetry = new Symmetry();
    }
    else {
      //Symmetry reading
      string typeSym(element->Attribute("type"));
      if (typeSym == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeSym);
      //Switch on the type of symmetry
      if (typeSym == "CYLINDRICAL") { m_run->m_symmetry = new SymCylindrical(element, fileName.str()); }
      else if (typeSym == "SPHERICAL") { m_run->m_symmetry = new SymSpherical(element, fileName.str()); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
    }

    //Source terms reading
    int sourceFound(0);
    element = xmlNode->FirstChildElement("sourceTerms");
    while (element != NULL) {
      sourceFound++;
      //Source reading
      string typeSource(element->Attribute("type"));
      if (typeSource == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeSource);
      //Switch on the type of source
      if (typeSource == "GRAVITY") { m_run->m_sources.push_back(new SourceGravity(element, fileName.str())); }
      else if (typeSource == "HEATING") { m_run->m_sources.push_back(new SourceHeating(element, fileName.str())); }
      else if (typeSource == "MRF") { m_run->m_MRF = m_run->m_sources.size(); m_run->m_sources.push_back(new SourceMRF(element, fileName.str())); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("sourceTerms");
    }
    m_run->m_numberSources = sourceFound;

    //Relaxations reading
    element = xmlNode->FirstChildElement("relaxation");
    while (element != NULL) {
      //Lecture source
      string typeRelax(element->Attribute("type"));
      if (typeRelax == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeRelax);
      //switch sur le type de relaxation
      if (typeRelax == "PTMU") { m_run->m_model->getRelaxations()->push_back(new RelaxationPTMu(element, fileName.str())); }
	  else if (typeRelax == "PT") {m_run->m_model->getRelaxations()->push_back(new RelaxationPT()); }
      else if (typeRelax == "P") { m_run->m_model->getRelaxations()->push_back(new RelaxationP()); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("relaxation");
    }

  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

Eos* Input::entreeEOS(string EOS, int &numberEOS)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    stringstream fileName("./libEOS/" + EOS);
    XMLDocument xmlEOS;
    XMLError error(xmlEOS.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees de l'EOS
    //------------------------------------
    //Recuperation racine du document XML
    XMLNode *xmlNode = xmlEOS.FirstChildElement("parametersEOS");
    if (xmlNode == NULL) throw ErrorXMLRacine("parametersEOS", fileName.str(), __FILE__, __LINE__);
    //Recuperation type d'EOS
    XMLElement *element;
    element = xmlNode->FirstChildElement("EOS");
    if (element == NULL) throw ErrorXMLElement("EOS", fileName.str(), __FILE__, __LINE__);
    string typeEOS(element->Attribute("type"));
    if (typeEOS == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(typeEOS);
    //Switch selon EOS
    Eos *eos;
    vector<string> NamesParametresEos;
    if      (typeEOS == "IG"){ eos = new EosIG(NamesParametresEos, numberEOS); }
    else if (typeEOS == "SG"){ eos = new EosSG(NamesParametresEos, numberEOS); }
    else{ throw ErrorXMLEOSInconnue(typeEOS, fileName.str(), __FILE__, __LINE__); } //Cas ou la loi state est inconnue

    //Recuperation des parametres de l'EOS
    element = xmlNode->FirstChildElement("parameters");
    if (element == NULL) throw ErrorXMLElement("parameters", fileName.str(), __FILE__, __LINE__);
    vector<double> parametresEos(NamesParametresEos.size());
    for (unsigned int p = 0; p < NamesParametresEos.size(); p++)
    {
      error = element->QueryDoubleAttribute(NamesParametresEos[p].c_str(), &parametresEos[p]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut(NamesParametresEos[p].c_str(), fileName.str(), __FILE__, __LINE__);
    }
    eos->assignParametersEos(EOS.c_str(), parametresEos);
    //Lecture des parametres physiques (viscosite, conductivite, etc.)
    eos->readPhysicalParameter(xmlNode, fileName.str());

    return eos;
  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Input::entreeConditionsInitiales(string casTest, vector<GeometricalDomain*> &domains, vector<BoundCond*> &boundCond)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    stringstream fileName(casTest + m_nameCI);
    XMLDocument xmlModel;
    XMLError error(xmlModel.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //------------------------------ CONDITIONS INITIALES -------------------------------

    //2) Recuperation racine CI du document XML
    //-----------------------------------------
    XMLNode *xmlNode = xmlModel.FirstChildElement("CI");
    if (xmlNode == NULL) throw ErrorXMLRacine("CI", fileName.str(), __FILE__, __LINE__);
    XMLElement *elementDomaine(xmlNode->FirstChildElement("physicalDomains"));
    if(elementDomaine == NULL) throw ErrorXMLElement("physicalDomains", fileName.str(), __FILE__, __LINE__);

    //3) Recuperation des domains definis
    //------------------------------------
    string nameDomaine, stateDomaine;
    int domaineTrouve(0);
    XMLElement *element(elementDomaine->FirstChildElement("domain"));
    while (element != NULL)
    {
      //A)Lecture name domain
      //*********************
      nameDomaine = element->Attribute("name");
      if (nameDomaine == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(nameDomaine);

      //B)Lecture state domain
      //**********************
      stateDomaine = element->Attribute("state");
      if (stateDomaine == "") throw ErrorXMLAttribut("state", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(stateDomaine);
      //Recherche de l state associe au domain
      XMLElement *state(xmlNode->FirstChildElement("state"));
      bool trouve(false); string nameEtat;
      while (state != NULL)
      {
        nameEtat = state->Attribute("name");
        Tools::uppercase(nameEtat);
        if (nameEtat == stateDomaine){ trouve = true; break; }
        state = state->NextSiblingElement("state");
      }
      if (!trouve){ throw ErrorXMLEtat(stateDomaine, fileName.str(), __FILE__, __LINE__); }

      //Reading mixture state first
      Mixture *stateMixture(0);
      if (m_run->m_model->whoAmI() == "EULER") { stateMixture = new MixEuler(); }
      else if (m_run->m_model->whoAmI() == "KAPILA") { stateMixture = new MixKapila(state, fileName.str()); }
      else if (m_run->m_model->whoAmI() == "MULTIP") { stateMixture = new MixMultiP(state, fileName.str()); }
      else if (m_run->m_model->whoAmI() == "THERMALEQ") { stateMixture = new MixThermalEq(state, fileName.str()); }
      else  if (m_run->m_model->whoAmI() == "EULERHOMOGENEOUS") { stateMixture = new MixEulerHomogeneous(state, fileName.str()); }
      else { throw ErrorXMLElement("Couplage Model-Fluide non valide", fileName.str(), __FILE__, __LINE__); }

      //Then Reading phases states
      vector<Phase*> statesPhases;
      string typeMateriau; string nameEOS;
      XMLElement* material(state->FirstChildElement("material"));
      int nbMateriauxEtat(0);
      while (material != NULL)
      {
        typeMateriau = material->Attribute("type");
        Tools::uppercase(typeMateriau);

        //LECTURE FLUIDE
        if (typeMateriau == "FLUIDE") { 
          nbMateriauxEtat++;
          //Recuperation de l EOS
          nameEOS = material->Attribute("EOS");
          int e(0);//  bool eosTrouvee(false);
          for (e = nbMateriauxEtat - 1; e < m_run->m_numberPhases; e++) {
            if (nameEOS == m_run->m_eos[e]->getName()) { break; }
          }
          if (e == m_run->m_numberPhases) { throw ErrorXMLEOSInconnue(nameEOS, fileName.str(), __FILE__, __LINE__); }
          if (e != nbMateriauxEtat-1) { throw ErrorXMLEtat("Materials in used states should be ordered in reference to modelXX.xml", fileName.str(), __FILE__, __LINE__); }
          if (m_run->m_model->whoAmI() == "EULER") { statesPhases.push_back(new PhaseEuler(material, m_run->m_eos[e], fileName.str())); }
          else if (m_run->m_model->whoAmI() == "KAPILA") { statesPhases.push_back(new PhaseKapila(material, m_run->m_eos[e], stateMixture->getPressure(), fileName.str())); }
          else if (m_run->m_model->whoAmI() == "MULTIP") { statesPhases.push_back(new PhaseMultiP(material, m_run->m_eos[e], fileName.str())); }
          else if (m_run->m_model->whoAmI() == "THERMALEQ") { statesPhases.push_back(new PhaseThermalEq(material, m_run->m_eos[e], fileName.str())); }
          else if (m_run->m_model->whoAmI() == "EULERHOMOGENEOUS") { statesPhases.push_back(new PhaseEulerHomogeneous(material, m_run->m_eos[e], fileName.str())); }
          else { throw ErrorXMLElement("Not valid Model-Fluide coupling", fileName.str(), __FILE__, __LINE__); }
        }
        
        //MATERIAU INCONNU
        else { throw ErrorXMLMateriauInconnu(typeMateriau, fileName.str(), __FILE__, __LINE__); } //Cas ou le type de material n a pas ete implemente
        material = material->NextSiblingElement("material");
      }
      if(nbMateriauxEtat!=m_run->m_numberPhases) throw ErrorXMLEtat(stateDomaine, fileName.str(), __FILE__, __LINE__);

      //Lecture des variables transportees pour l'state trouve
      vector<Transport> statesTransport(m_run->m_numberTransports);
      string nameTransport; double valueTransport(0.);
      XMLElement *elementTransport(state->FirstChildElement("transport"));
      int nbTransports(0);
      while (elementTransport != NULL) {
        nameTransport = elementTransport->Attribute("name");
        elementTransport->QueryDoubleAttribute("value", &valueTransport);
        int e(0);
        for (e = 0; e < m_run->m_numberTransports; e++) {
          if (nameTransport == m_run->m_nameGTR[e]) { break; }
        }
        if (e != m_run->m_numberTransports) {
          statesTransport[e].setValue(valueTransport);
          nbTransports++;
        }
        elementTransport = elementTransport->NextSiblingElement("transport");
      }

      //C)Reading optional physical entity
      //**********************************
      int physicalEntity(element->IntAttribute("physicalEntity")); //Default value: 0      

      //C)Reading domain type
      //*********************
      string typeDomaine(element->Attribute("type"));
      Tools::uppercase(typeDomaine);
      vector<string> NamesParametresDomaine;
      if      (typeDomaine == "ENTIREDOMAIN")                    { domains.push_back(new GDEntireDomain(nameDomaine, statesPhases, stateMixture, statesTransport, physicalEntity)); }
      else if (typeDomaine == "ENTIREDOMAINWITHPARTICULARITIES") { domains.push_back(new GDEntireDomainWithParticularities(nameDomaine, statesPhases, stateMixture, statesTransport, physicalEntity)); }
      else if (typeDomaine == "HALFSPACE")                       { domains.push_back(new GDHalfSpace(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "DISC")                            { domains.push_back(new GDDisc(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "ELLIPSE")                         { domains.push_back(new GDEllipse(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "RECTANGLE")                       { domains.push_back(new GDRectangle(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "PAVEMENT")                        { domains.push_back(new GDPavement(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "SPHERE")                          { domains.push_back(new GDSphere(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else{ throw ErrorXMLDomaineInconnu(typeDomaine, fileName.str(), __FILE__, __LINE__); } //Cas ou le domain n a pas ete implemente
      domaineTrouve++;
      //Domaine suivant
      element = element->NextSiblingElement("domain");

      for (int k = 0; k < m_run->m_numberPhases; k++) { delete statesPhases[k]; }
      delete stateMixture;
    }

    //------------------------------ CONDITIONS AUX LIMITES -------------------------------

    //4) Recuperation racine CL du document XML
    //-----------------------------------------
    XMLElement *elementLimites(xmlNode->FirstChildElement("boundaryConditions"));
    if (elementLimites == NULL) throw ErrorXMLElement("boundaryConditions", fileName.str(), __FILE__, __LINE__);

    //5) Recuperation des boundCond definies
    //------------------------------------
    string nameBoundCond, stateLimite, typeBoundCond;
    int numBoundCond;
    int boundCondTrouve(0);
    element = elementLimites->FirstChildElement("boundCond");
    while (element != NULL)
    {
      //A)Lecture name boundCond
      //************************
      nameBoundCond = element->Attribute("name");
      if (nameBoundCond == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(nameBoundCond);

      //B)Lecture type boundCond
      //************************
      typeBoundCond = element->Attribute("type");
      if (typeBoundCond == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeBoundCond);

      //C)Number de la boundCond associe a la geometrie mesh
      //****************************************************
      error = element->QueryIntAttribute("number", &numBoundCond);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName.str(), __FILE__, __LINE__);

      //D)Reading conLim specific data
      //******************************
      if (typeBoundCond == "ABS") { boundCond.push_back(new BoundCondAbs(numBoundCond)); }
      else if (typeBoundCond == "SYMMETRY") { if (m_run->m_order == "FIRSTORDER") { boundCond.push_back(new BoundCondWall(numBoundCond)); } else { boundCond.push_back(new BoundCondSymmetryO2(numBoundCond)); } }
      else if (typeBoundCond == "WALL") { if (m_run->m_order == "FIRSTORDER") { boundCond.push_back(new BoundCondWall(numBoundCond)); } else { boundCond.push_back(new BoundCondWallO2(numBoundCond)); } }
      else if (typeBoundCond == "INJECTION") { boundCond.push_back(new BoundCondInj(numBoundCond, element, m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_nameGTR, m_run->m_eos, fileName.str())); }
      else if (typeBoundCond == "TANK") { boundCond.push_back(new BoundCondTank(numBoundCond, element, m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_nameGTR, m_run->m_eos, fileName.str())); }
      else if (typeBoundCond == "OUTFLOW") { boundCond.push_back(new BoundCondOutflow(numBoundCond, element, m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_nameGTR, fileName.str())); }
      else { throw ErrorXMLBoundCondInconnue(typeBoundCond, fileName.str(), __FILE__, __LINE__); } //Cas ou la limite n a pas ete implemente
      boundCondTrouve++;
      //Domaine suivant
      element = element->NextSiblingElement("boundCond");
    }

  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************
