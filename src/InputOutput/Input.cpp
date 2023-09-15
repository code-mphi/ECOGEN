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

#include "Input.h"
#include "HeaderInputOutput.h"
#include "../Meshes/stretchZone.h"
#include "../Config.h"

using namespace tinyxml2;

//***********************************************************************

Input::Input(Run *run) : m_run(run)
{
	m_nameMain = "main.xml";
	m_nameMesh = "mesh.xml";
	m_nameCI = "initialConditions.xml";
	m_nameModel = "model.xml";
}

//***********************************************************************

Input::~Input(){}

//***********************************************************************

void Input::readInputXML(std::vector<GeometricalDomain*>& domains, std::vector<BoundCond*>& boundCond, std::vector<GeometricalDomain*>& solidDomains)
{
  try{
    //1) General computation parameters
    inputMain(m_run->m_simulationName);
    //2) Mesh data
    inputMesh(m_run->m_simulationName);
    //3) Model and fluid data
    inputModel(m_run->m_simulationName);
    //4) Read initial conditions
    inputInitialConditions(m_run->m_simulationName, domains, boundCond, solidDomains);
    
    // Check compatibility of some specific input options
    verifyCompatibilityInput(m_run->m_simulationName);
  }
  catch (ErrorInput &e){
    if(rankCpu==0) std::cerr << e.infoError() << std::endl;
    throw;
  }
}

//***********************************************************************

void Input::inputMain(std::string casTest)
{
  try{
    //1) Parsing XML file with the library tinyxml2
    //---------------------------------------------
	  std::stringstream fileName(casTest + m_nameMain);
    XMLDocument xmlMain;
    XMLError error(xmlMain.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(),__FILE__, __LINE__);
    
    //2) Get main compute data
    //------------------------
    //Get root of XML document
    XMLNode *computationParam = xmlMain.FirstChildElement("computationParam");
    if (computationParam == NULL) throw ErrorXMLRacine("computationParam", fileName.str(), __FILE__, __LINE__);

    XMLElement* element, *sousElement;

    //Get run name
    element = computationParam->FirstChildElement("run");
    if (element == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);
    XMLNode* xmlNode2 = element->FirstChild();
    if (xmlNode2 == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);
    XMLText* xmlText = xmlNode2->ToText();
    if (xmlText == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);

    //Read des informations de sorties/prints
    element = computationParam->FirstChildElement("outputMode");
    if (element == NULL) throw ErrorXMLElement("outputMode", fileName.str(), __FILE__, __LINE__);
    //Read format sortie
    std::string format(element->Attribute("format"));
    if (format == "") throw ErrorXMLAttribut("format", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(format);
    if (format == "XML") { m_run->m_outPut = new OutputXML(casTest, xmlText->Value(), element, fileName.str(), this); }
    else if (format == "GNU") { m_run->m_outPut = new OutputGNU(casTest, xmlText->Value(), element, fileName.str(), this); }
    else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }

    //Read des cuts 1D
    element = computationParam->FirstChildElement("cut1D");
    while (element != NULL)
    {
      m_run->m_cuts.push_back(new OutputCutGNU(casTest, xmlText->Value(), element, fileName.str(), LINE, this));
      element = element->NextSiblingElement("cut1D");
    }
    //Read des cuts 2D
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

    //Record global quantity on whole domain (mass, totalEnergy)
    element = computationParam->FirstChildElement("globalQuantity");
    while (element != NULL) {
      std::string quantity(element->Attribute("quantity"));
      if (quantity == "") throw ErrorXMLAttribut("quantity", fileName.str(), __FILE__, __LINE__);
      Tools::lowercase(quantity);
      if (quantity == "mass" || quantity == "totalenergy") { 
        m_run->m_globalQuantities.push_back(new OutputGlobalGNU(casTest, xmlText->Value(), element, this, quantity)); 
      }
      else { throw ErrorXMLAttribut("quantity", fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("globalQuantity");
    }

    //Record saturation pressure
    //Be careful when using it, not fully generalized
    element = computationParam->FirstChildElement("psat");
    if (element != NULL) {
      bool recordPsat(false);
      error = element->QueryBoolAttribute("record", &recordPsat);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("record", fileName.str(), __FILE__, __LINE__);
      m_run->m_recordPsat = recordPsat;
    }

    //Record massflow on a given boundary
    element = computationParam->FirstChildElement("boundary");
    while (element != NULL)
    {
      std::string recordBoundType(element->Attribute("record"));
      if (recordBoundType == "") throw ErrorXMLAttribut("record", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(recordBoundType);
      if (recordBoundType == "FLUX") {
        m_run->m_recordBoundaries.push_back(new OutputBoundaryFluxGNU(casTest, xmlText->Value(), element, fileName.str(), this));
      }
      else if (recordBoundType == "ALL") {
        m_run->m_recordBoundaries.push_back(new OutputBoundaryAllGNU(casTest, xmlText->Value(), element, fileName.str(), this));
      }
      else {
        throw ErrorXMLAttribut("record", fileName.str(), __FILE__, __LINE__);
      }
      element = element->NextSiblingElement("boundary");
    }
    
    //Get Iteration / temps Physique
    element = computationParam->FirstChildElement("timeControlMode");
    if (element == NULL) throw ErrorXMLElement("timeControlMode", fileName.str(), __FILE__, __LINE__);
    error = element->QueryBoolAttribute("iterations", &m_run->m_timeControlIterations);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("iterations", fileName.str(), __FILE__, __LINE__);
    if (m_run->m_timeControlIterations)
    {
      //Get Iterations / Frequence
      sousElement = element->FirstChildElement("iterations");
      if (sousElement == NULL) throw ErrorXMLElement("iterations", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryIntAttribute("number", &m_run->m_nbIte);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryIntAttribute("iterFreq", &m_run->m_freq);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("iterFreq", fileName.str(), __FILE__, __LINE__);
    }
    else
    {
      //Get Temps / Frequence
      sousElement = element->FirstChildElement("physicalTime");
      if (sousElement == NULL) throw ErrorXMLElement("physicalTime", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryDoubleAttribute("totalTime", &m_run->m_finalPhysicalTime);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("totalTime", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryDoubleAttribute("timeFreq", &m_run->m_timeFreq);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("timeFreq", fileName.str(), __FILE__, __LINE__);
    }

    //Get CFL
    element = computationParam->FirstChildElement("computationControl");
    if (element == NULL) throw ErrorXMLElement("computationControl", fileName.str(), __FILE__, __LINE__);
    error = element->QueryDoubleAttribute("CFL", &m_run->m_cfl);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("CFL", fileName.str(), __FILE__, __LINE__);

    //Reading gradient method
    element = computationParam->FirstChildElement("gradient");
    if (element != NULL) {
      sousElement = element->FirstChildElement("method");
      if (sousElement == NULL) throw ErrorXMLElement("method", fileName.str(), __FILE__, __LINE__);
      XMLNode *method = sousElement->FirstChild();
      if (method == NULL) throw ErrorXMLElement("method", fileName.str(), __FILE__, __LINE__);
      std::string methodName = method->ToText()->Value();
      Tools::uppercase(methodName);
      if (methodName == "FINITE-DIFFERENCE") { m_run->m_gradient = new GradientFiniteDifference(); }
      else if (methodName == "GREEN-GAUSS") { m_run->m_gradient = new GradientGreenGauss(); }
      else { throw ErrorXMLElement("method", fileName.str(), __FILE__, __LINE__); }
    }
    else {
      m_run->m_gradient = new GradientFiniteDifference();
    }

    //Read Ordre2
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
      //Get global limiter
      sousElement = element->FirstChildElement("globalLimiter");
      if (sousElement == NULL) throw ErrorXMLElement("globalLimiter", fileName.str(), __FILE__, __LINE__);
      contenu = sousElement->FirstChild();
      if (contenu == NULL) throw ErrorXMLElement("globalLimiter", fileName.str(), __FILE__, __LINE__);
      std::string globalLimiter = contenu->ToText()->Value();
      Tools::uppercase(globalLimiter);
      if (globalLimiter == "MINMOD") { m_run->m_globalLimiter = new LimiterMinmod; }
      else if (globalLimiter == "VANLEER") { m_run->m_globalLimiter = new LimiterVanLeer; }
      else if (globalLimiter == "VANALBADA") { m_run->m_globalLimiter = new LimiterVanAlbada; }
      else if (globalLimiter == "SUPERBEE") { m_run->m_globalLimiter = new LimiterSuperBee; }
      else if (globalLimiter == "MC") { m_run->m_globalLimiter = new LimiterMC; }
      else if (globalLimiter == "THINC") { throw ErrorXMLAttribut("THINC can only be a volume-fraction limiter", fileName.str(), __FILE__, __LINE__); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      //Get interface limiter
      std::string interfaceLimiter = globalLimiter;
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
      //Get global volume-fraction limiter
      std::string globalVolumeFractionLimiter = globalLimiter;
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
      //Get interface volume-fraction limiter
      std::string interfaceVolumeFractionLimiter = interfaceLimiter;
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

    //Reprise de calcul depuis file resultat
    element = computationParam->FirstChildElement("restartSimulation");
    if (element != NULL) { 
      error = element->QueryIntAttribute("restartFileNumber", &m_run->m_restartSimulation);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("restartFileNumber", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("AMRsaveFreq", &m_run->m_restartAMRsaveFreq);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("AMRsaveFreq", fileName.str(), __FILE__, __LINE__);
    }

  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Input::inputMesh(std::string casTest)
{
  try{
    //Methode AMR: initialization des variables
    m_run->m_lvlMax = 0;
		double criteriaVar(1.e10);
		bool varRho(false), varP(false), varU(false), varAlpha(false);
		double xiSplit(1.), xiJoin(1.);

    //1) Parsing XML file with the library tinyxml2
    //---------------------------------------------
    std::stringstream fileName(casTest + m_nameMesh);
    XMLDocument xmlMesh;
    XMLError error(xmlMesh.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Get main mesh data
    //---------------------
    //Get root of XML document
    XMLNode *mesh = xmlMesh.FirstChildElement("mesh");
    if (mesh == NULL) throw ErrorXMLRacine("mesh", fileName.str(), __FILE__, __LINE__);

    XMLElement* element;
    //Get mesh type
    element = mesh->FirstChildElement("type");
    if (element == NULL) throw ErrorXMLElement("type", fileName.str(), __FILE__, __LINE__);
    std::string structureMesh(element->Attribute("structure"));
    if (structureMesh == "") throw ErrorXMLAttribut("structure", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(structureMesh);
    if (structureMesh == "UNSTRUCTURED")
    {
      if (m_run->m_order == "SECONDORDER") {
        Errors::errorMessage("Second-order method for unstructured meshes is not released yet because of bugs");
        exit(0);
      }

      //----------------MESH NON STRUCTURE ---------------------
      XMLElement* meshNS;
      meshNS = mesh->FirstChildElement("unstructuredMesh");
      if (meshNS == NULL) throw ErrorXMLElement("unstructuredMesh", fileName.str(), __FILE__, __LINE__);
      //Read name du file contenant les informations de mesh
      element = meshNS->FirstChildElement("file");
      if (element == NULL) throw ErrorXMLElement("file", fileName.str(), __FILE__, __LINE__);
      std::string meshFile(element->Attribute("name"));
      if (meshFile == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      //Read format of the mesh file (.msh, .mesh ...)
      meshFile = config.getWorkFolder() + meshFile;

      std::string meshExtension(MeshUnStruct::readMeshFileExtension(meshFile));
      Tools::uppercase(meshExtension);
      if (meshExtension == "MSH") // Gmsh format
      {
        std::string version(MUSGmsh::readVersion(meshFile));
        if (version == "2.2") {
			    bool switchTags(false);
			    element = meshNS->FirstChildElement("tag");
			    if (element != NULL) {
				    error = element->QueryBoolAttribute("GMSHSwitchTags", &switchTags);
				    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("GMSHSwitchTags", fileName.str(), __FILE__, __LINE__);
			    }
			    m_run->m_mesh = new MUSGmshV2(meshFile, meshExtension, switchTags);
		    }
        else if (version == "4.1") m_run->m_mesh = new MUSGmshV4(meshFile, meshExtension);
        else throw ErrorXML("mesh version not found for file : " + meshFile, __FILE__, __LINE__);
      }
      else if (meshExtension == "MESH") throw ErrorXML("MESH format is not supported for file : " + meshFile, __FILE__, __LINE__);
      else { throw ErrorXML("mesh extension not supported for file : " + meshFile, __FILE__, __LINE__); }
      
      //Get pretraitement parallele
      element = meshNS->FirstChildElement("parallel");
      if (element != NULL) {
        error = element->QueryBoolAttribute("GMSHPretraitement", &m_run->m_parallelPreTreatment);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("GMSHPretraitement", fileName.str(), __FILE__, __LINE__);
      }
      //AMR method not possible with unstructured mesh
      element = meshNS->FirstChildElement("AMR");
      if (element != NULL) { throw ErrorXMLAttribut("Methode AMR non possible avec mesh non structure", fileName.str(), __FILE__, __LINE__); }

      //Mesh analysis option
      element = meshNS->FirstChildElement("extractReferenceLength");
      if (element != NULL) {
        error = element->QueryBoolAttribute("state", &m_run->m_extractRefLength);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("extractReferenceLength", fileName.str(), __FILE__, __LINE__);
      }

      //Mesh mapping restart option
      element = meshNS->FirstChildElement("meshMappingRestart");
      if (element != NULL) {
        // Read attributes
        std::string resultFolderMapped = element->Attribute("resultFolder");
        if (resultFolderMapped == "") throw ErrorXMLAttribut("resultFolder", fileName.str(), __FILE__, __LINE__);
        int restartFileNumber(0);
        error = element->QueryIntAttribute("restartFileNumber", &restartFileNumber);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("restartFileNumber", fileName.str(), __FILE__, __LINE__);
        m_run->m_meshFileMapped = element->Attribute("meshFile");
        if (m_run->m_meshFileMapped == "") throw ErrorXMLAttribut("meshFile", fileName.str(), __FILE__, __LINE__);
        m_run->m_meshFileMapped = config.getWorkFolder() + m_run->m_meshFileMapped;
        // Check if mapped file is conform
        std::string meshMappedExtension(MeshUnStruct::readMeshFileExtension(m_run->m_meshFileMapped));
        Tools::uppercase(meshMappedExtension);
        if (meshMappedExtension == "MSH" && meshMappedExtension == meshExtension) // Gmsh format & same on both simulation
        {
          std::string versionMapped(MUSGmsh::readVersion(m_run->m_meshFileMapped));
          if (versionMapped == "2.2") {
            m_run->m_restartMeshMapping = true;
            m_run->m_outputMeshMapping = new OutputXML(resultFolderMapped, restartFileNumber, this);
          }
          else throw ErrorXML("mesh version not supported for restart with mesh mapping, see file: " + m_run->m_meshFileMapped, __FILE__, __LINE__);
        }
        else { throw ErrorXML("mesh format not supported for restart with mesh mapping, see file: " + m_run->m_meshFileMapped, __FILE__, __LINE__); }
      }
    }
    else if (structureMesh == "CARTESIAN")
    {
      //----------------MESH CARTESIAN ---------------------
      XMLElement* cartesianMesh;
      cartesianMesh = mesh->FirstChildElement("cartesianMesh");
      if (cartesianMesh == NULL) throw ErrorXMLElement("cartesianMesh", fileName.str(), __FILE__, __LINE__);
      //Get dimensions
      double lX, lY, lZ;
      element = cartesianMesh->FirstChildElement("dimensions");
      if (element == NULL) throw ErrorXMLElement("dimensions", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("x", &lX);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("y", &lY);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("z", &lZ);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName.str(), __FILE__, __LINE__);
      //Get number of cells
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
      std::vector<stretchZone> stretchX, stretchY, stretchZ;
      element = cartesianMesh->FirstChildElement("meshStretching");
      if (element != NULL) {
        double tempStart, tempEnd, tempFactor;
        int tempNumberCells;
        XMLElement* sousElement;
        //X stretching
        if (nbX > 1) {
          sousElement = element->FirstChildElement("XStretching");
          if (sousElement != NULL) {
            XMLElement* stretchElement(sousElement->FirstChildElement("stretch"));
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
            XMLElement* stretchElement(sousElement->FirstChildElement("stretch"));
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
            XMLElement* stretchElement(sousElement->FirstChildElement("stretch"));
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

      //Get des variables pour methode AMR
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

void Input::inputModel(std::string casTest)
{
  try{
    //1) Parsing XML file with the library tinyxml2
    //---------------------------------------------
    std::stringstream fileName(casTest + m_nameModel);
    XMLDocument xmlModel;
    XMLError error(xmlModel.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Get model data
    //-----------------
    //Get root of XML document
    XMLNode *xmlNode = xmlModel.FirstChildElement("model");
    if (xmlNode == NULL) throw ErrorXMLRacine("model", fileName.str(), __FILE__, __LINE__);
    //Get name of solved model
    XMLElement* element(xmlNode->FirstChildElement("flowModel"));
    if (element == NULL) throw ErrorXMLElement("flowModel", fileName.str(), __FILE__, __LINE__);
    std::string model(element->Attribute("name"));
    if (model == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
    m_run->m_numberTransports = 0;
    error = element->QueryIntAttribute("numberTransports", &m_run->m_numberTransports);
    //if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberTransports", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(model);
    //Switch selon model
    bool alphaNull(false);
    if (model == "EULER"){ 
      m_run->m_model = new ModEuler(m_run->m_numberTransports); m_run->m_numberPhases = 1; 
    }
    else if (model == "EULERHOMOGENEOUS") {
      int liquid, vapor;
      error = element->QueryIntAttribute("liquid", &liquid);
      error = element->QueryIntAttribute("vapor", &vapor);
      m_run->m_model = new ModEulerHomogeneous(m_run->m_numberTransports, liquid, vapor); m_run->m_numberPhases = 2; 
    }
    else if (model == "PRESSUREVELOCITYEQ")
    { 
      error = element->QueryIntAttribute("numberPhases", &m_run->m_numberPhases);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberPhases", fileName.str(), __FILE__, __LINE__);
      m_run->m_model = new ModPUEq(m_run->m_numberTransports, m_run->m_numberPhases);
      error = element->QueryBoolAttribute("alphaNull", &alphaNull);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alphaNull", fileName.str(), __FILE__, __LINE__);
      if (m_run->m_outPut->getReducedOutput()) { 
        numberScalarsMixture = 2; //Density and pressure
        numberScalarsPhase = 2; //Volume fraction and density
      }
      else {
        numberScalarsMixture = 3; //Density, pressure and total energy
        numberScalarsPhase = 4; //Volume fraction, density, temperature and mass fraction
      }
    }
    else if (model == "VELOCITYEQ")
    {
      error = element->QueryIntAttribute("numberPhases", &m_run->m_numberPhases);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberPhases", fileName.str(), __FILE__, __LINE__);
      m_run->m_model = new ModUEq(m_run->m_numberTransports, m_run->m_numberPhases);
      error = element->QueryBoolAttribute("alphaNull", &alphaNull);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alphaNull", fileName.str(), __FILE__, __LINE__);
      if (m_run->m_outPut->getReducedOutput()) { 
        numberScalarsMixture = 2; //Density and pressure
        numberScalarsPhase = 3; //Volume fraction, density and pressure
      }
      else {
        numberScalarsMixture = 3; //Density, pressure and total energy
        numberScalarsPhase = 5; //Volume fraction, density, pressure, temperature and mass fraction
      }
    }
    else if (model == "VELOCITYEQTOTENERGY") 
    {
      error = element->QueryIntAttribute("numberPhases", &m_run->m_numberPhases);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberPhases", fileName.str(), __FILE__, __LINE__);
      m_run->m_model = new ModUEqTotEnergy(m_run->m_numberTransports, m_run->m_numberPhases);
      error = element->QueryBoolAttribute("alphaNull", &alphaNull);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alphaNull", fileName.str(), __FILE__, __LINE__);
      numberScalarsMixture = 2; //Density and pressure
    }
    else if (model == "TEMPERATUREPRESSUREVELOCITYEQ")
    {
      error = element->QueryIntAttribute("numberPhases", &m_run->m_numberPhases);
      m_run->m_model = new ModPTUEq(m_run->m_numberTransports, m_run->m_numberPhases);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberPhases", fileName.str(), __FILE__, __LINE__);
    }
    else if (model == "NONLINEARSCHRODINGER")
    {
      double alpha(0.), beta(0.);
      error = element->QueryDoubleAttribute("alpha", &alpha);
      error = element->QueryDoubleAttribute("beta", &beta);
      m_run->m_model = new ModNonLinearSchrodinger(m_run->m_numberTransports, alpha, beta);
      m_run->m_numberPhases = 1;
    }
    else if (model == "EULERKORTEWEG")
    {
      double alpha(0.), beta(0.), temperature(0.), kappa(0.);
      error = element->QueryDoubleAttribute("alpha", &alpha);
      error = element->QueryDoubleAttribute("beta", &beta);
      error = element->QueryDoubleAttribute("temperature", &temperature);
      error = element->QueryDoubleAttribute("kappa", &kappa);
      m_run->m_model = new ModEulerKorteweg(m_run->m_numberTransports, alpha, beta, temperature, kappa);
      m_run->m_numberPhases = 1;
    }
    else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
    if (m_run->m_globalVolumeFractionLimiter->AmITHINC() && (m_run->m_numberPhases != 2)) { throw ErrorXMLAttribut("Limiter THINC only works for 2 phases", fileName.str(), __FILE__, __LINE__); }
    if (m_run->m_interfaceVolumeFractionLimiter->AmITHINC() && (m_run->m_numberPhases != 2)) { throw ErrorXMLAttribut("Limiter THINC only works for 2 phases", fileName.str(), __FILE__, __LINE__); }
    
    //Low-Mach preconditioning
    element = xmlNode->FirstChildElement("lowMach");
    if (element != NULL) {
      bool lowMach(false);
      error = element->QueryBoolAttribute("state", &lowMach);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("state", fileName.str(), __FILE__, __LINE__);
      m_run->m_model->setLowMach(lowMach);
      if (lowMach) {
        double machRefMin(1.e-2);
        error = element->QueryDoubleAttribute("machRefMin", &machRefMin);
        if (error == XML_NO_ERROR) {
          m_run->m_model->setMachRefMin(machRefMin);
        }
      }
    }

    //Thermodynamique
    std::vector<std::string> nameEOS;
    if (model != "NONLINEARSCHRODINGER") {
      int EOSTrouvee(0);
      XMLElement* sousElement(xmlNode->FirstChildElement("EOS"));
      while (sousElement != NULL)
      {
        //Read Eos
        nameEOS.push_back(sousElement->Attribute("name"));
        if (nameEOS[EOSTrouvee] == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
        EOSTrouvee++;
        sousElement = sousElement->NextSiblingElement("EOS");
      }
      m_run->m_numberEos = nameEOS.size();
      if (m_run->m_numberEos == 0) throw ErrorXMLEOS(fileName.str(), __FILE__, __LINE__);
      m_run->m_eos = new Eos*[m_run->m_numberPhases];
      int numberEos(0);
      for (int i = 0; i < m_run->m_numberEos; i++){ m_run->m_eos[i] = inputEOS(nameEOS[i], numberEos); }
      m_run->m_eos[0]->assignEpsilonForAlphaNull(alphaNull); //initialization of epsilon value for alphaNull option
    }

    //Transports
    int transportsTrouvee(0);
    XMLElement* elementTransport(xmlNode->FirstChildElement("transport"));
    while (elementTransport != NULL)
    {
      //Read name transport
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
      //Read physique additionnelle
      std::string typeAddPhys(element->Attribute("type"));
      if (typeAddPhys == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeAddPhys);
      //switch sur le type de physique additionelle
      if (typeAddPhys == "SURFACETENSION") { 
        //switch model d ecoulement
        if (model == "PRESSUREVELOCITYEQ" || model == "VELOCITYEQ") {
          m_run->m_addPhys.push_back(new APUEqSurfaceTension(element, numberGPA, m_run->m_nameGTR, nameEOS, fileName.str()));
        }
        else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      }
      else if (typeAddPhys == "VISCOSITY") {
        //switch model d ecoulement
        if (model == "PRESSUREVELOCITYEQ" || model == "VELOCITYEQ") {
          m_run->m_addPhys.push_back(new APUEqViscosity(numberGPA, m_run->m_eos, m_run->m_numberPhases));
        }
        else if (model == "EULER") { m_run->m_addPhys.push_back(new APEViscosity(numberGPA, m_run->m_eos)); }
        else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
        m_run->m_viscous = true;
      }
      else if (typeAddPhys == "CONDUCTIVITY") {
        //switch model d ecoulement
        if (model == "PRESSUREVELOCITYEQ" || model == "VELOCITYEQ") {
          m_run->m_addPhys.push_back(new APUEqConductivity(numberGPA, m_run->m_eos, m_run->m_numberPhases));
        }
        else if (model == "EULER") { m_run->m_addPhys.push_back(new APEConductivity(numberGPA, m_run->m_eos)); }
        else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("additionalPhysic");
    }
    m_run->m_numberAddPhys = physiqueAddTrouvee;

    //Symmetry terms reading
    element = xmlNode->FirstChildElement("symmetryTerm");
    if (element == NULL) {
      m_run->m_symmetry = new Symmetry();
    }
    else {
      //Symmetry reading
      std::string typeSym(element->Attribute("type"));
      if (typeSym == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeSym);
      //Switch on the type of symmetry
      if (typeSym == "CYLINDRICAL") { m_run->m_symmetry = new SymCylindrical(element, fileName.str()); }
      else if (typeSym == "SPHERICAL") { m_run->m_symmetry = new SymSpherical(element, fileName.str()); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
    }

    //1D geometry with smooth cross section variation
    element = xmlNode->FirstChildElement("geometry");
    if (element != NULL) {
      bool isSmoothCrossSection1d(false);
      error = element->QueryBoolAttribute("smoothCrossSection1d", &isSmoothCrossSection1d);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("smoothCrossSection1d", fileName.str(), __FILE__, __LINE__);
      m_run->m_model->setSmoothCrossSection1d(isSmoothCrossSection1d);
    }

    //Source terms reading
    int sourceFound(0);
    element = xmlNode->FirstChildElement("sourceTerms");
    while (element != NULL) {
      sourceFound++;
      //Source reading
      std::string typeSource(element->Attribute("type"));
      if (typeSource == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeSource);
      std::string orderSource(element->Attribute("order"));
      Tools::uppercase(orderSource);
      //Order
      int order = 1;
      if (orderSource == "EULER") { order = 1; }
      else if (orderSource == "RK2") { order = 2; }
      else if (orderSource == "RK4") { order = 4; }
      int physicalEntity(0);
      error = element->QueryIntAttribute("physicalEntity", &physicalEntity);
      if (error != XML_NO_ERROR) { physicalEntity = 0; } //Default = all cells have this source term
      //Switch on the type of source
      if (typeSource == "GRAVITY") { m_run->m_sources.push_back(new SourceNumGravity(element, order, physicalEntity, fileName.str())); }
      else if (typeSource == "HEATING") { m_run->m_sources.push_back(new SourceNumHeating(element, order, physicalEntity, fileName.str())); }
      else if (typeSource == "MRF") { m_run->m_MRF = m_run->m_sources.size(); m_run->m_sources.push_back(new SourceNumMRF(element, order, physicalEntity, fileName.str())); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("sourceTerms");
    }
    if (model == "NONLINEARSCHRODINGER" || model == "EULERKORTEWEG") {
      sourceFound++;
      int physicalEntity(0); //Default = all cells have this source term
      m_run->m_sources.push_back(new SourceExactEulerKorteweg(physicalEntity));
    }
    m_run->m_numberSources = sourceFound;

    //Relaxations reading
    element = xmlNode->FirstChildElement("relaxation");
    while (element != NULL) {
      //Read source
      std::string typeRelax(element->Attribute("type"));
      if (typeRelax == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeRelax);
      //Switch on the relaxation type
      if (typeRelax == "PTMU") {
        //Verify if relaxation not already added
        for (unsigned int r = 0; r < m_run->m_model->getRelaxations()->size(); r++) {
          if ((*m_run->m_model->getRelaxations())[r]->getType() == TypeRelax::PTMU) { throw ErrorXMLRelaxation(typeRelax, fileName.str(), __FILE__, __LINE__); }
        }
        m_run->m_model->getRelaxations()->push_back(new RelaxationPTMu(element, nameEOS, fileName.str()));
      }
	    else if (typeRelax == "PT") {
        //Verify if relaxation not already added
        for (unsigned int r = 0; r < m_run->m_model->getRelaxations()->size(); r++) {
          if ((*m_run->m_model->getRelaxations())[r]->getType() == TypeRelax::PT) { throw ErrorXMLRelaxation(typeRelax, fileName.str(), __FILE__, __LINE__); }
        }
        m_run->m_model->getRelaxations()->push_back(new RelaxationPT());
      }
      else if (typeRelax == "P") {
        //Verify if relaxation not already added
        for (unsigned int r = 0; r < m_run->m_model->getRelaxations()->size(); r++) {
          if ((*m_run->m_model->getRelaxations())[r]->getType() == TypeRelax::P) { throw ErrorXMLRelaxation(typeRelax, fileName.str(), __FILE__, __LINE__); }
        }
        std::string speedRelax(element->Attribute("speed"));
        Tools::uppercase(speedRelax);
        if (speedRelax == "FINITE") {
          m_run->m_model->getRelaxations()->push_back(new RelaxationPFinite(element, fileName.str()));
        }
        else { //Infinite
          m_run->m_model->getRelaxations()->push_back(new RelaxationPInfinite());
        }
      }
      else { throw ErrorXMLRelaxation(typeRelax, fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("relaxation");
    }

  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

Eos* Input::inputEOS(std::string EOS, int& numberEOS)
{
  try{
    //1) Parsing XML file with the library tinyxml2
    //---------------------------------------------
    std::stringstream fileName(config.getWorkFolder() + "./libEOS/" + EOS);
    XMLDocument xmlEOS;
    XMLError error(xmlEOS.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Get EOS data
    //---------------
    //Get root of XML document
    XMLNode *xmlNode = xmlEOS.FirstChildElement("parametersEOS");
    if (xmlNode == NULL) throw ErrorXMLRacine("parametersEOS", fileName.str(), __FILE__, __LINE__);
    //Get EOS type
    XMLElement* element;
    element = xmlNode->FirstChildElement("EOS");
    if (element == NULL) throw ErrorXMLElement("EOS", fileName.str(), __FILE__, __LINE__);
    std::string typeEOS(element->Attribute("type"));
    if (typeEOS == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(typeEOS);
    //Switch selon EOS
    Eos* eos;
    std::vector<std::string> NamesParametresEos;
    if      (typeEOS == "IG"){ eos = new EosIG(NamesParametresEos, numberEOS); }
    else if (typeEOS == "SG"){ eos = new EosSG(NamesParametresEos, numberEOS); }
    else if (typeEOS == "NASG"){ eos = new EosNASG(NamesParametresEos, numberEOS); }
    else if (typeEOS == "VDW"){ eos = new EosVDW(NamesParametresEos, numberEOS); }
    else if (typeEOS == "POLYNOMIAL"){ eos = new EosPolynomial(NamesParametresEos, numberEOS); }
    else{ throw ErrorXMLEOSUnknown(typeEOS, fileName.str(), __FILE__, __LINE__); } //Cas ou la loi state est unknown

    //Get of EOS parameters
    element = xmlNode->FirstChildElement("parameters");
    if (element == NULL) throw ErrorXMLElement("parameters", fileName.str(), __FILE__, __LINE__);
    std::vector<double> parametresEos(NamesParametresEos.size());
    for (unsigned int p = 0; p < NamesParametresEos.size(); p++)
    {
      error = element->QueryDoubleAttribute(NamesParametresEos[p].c_str(), &parametresEos[p]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut(NamesParametresEos[p].c_str(), fileName.str(), __FILE__, __LINE__);
    }
    eos->assignParametersEos(EOS.c_str(), parametresEos);
    //Read physical parameters (viscosity, thermal conductivity, etc.)
    eos->readPhysicalParameters(xmlNode);

    return eos;
  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Input::inputInitialConditions(std::string casTest, std::vector<GeometricalDomain*>& domains, std::vector<BoundCond*>& boundCond, std::vector<GeometricalDomain*>& solidDomains)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    std::stringstream fileName(casTest + m_nameCI);
    XMLDocument xmlModel;
    XMLError error(xmlModel.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //------------------------------ CONDITIONS INITIALES -------------------------------

    //2) Get root CI of XML document
    //------------------------------
    XMLNode *xmlNode = xmlModel.FirstChildElement("CI");
    if (xmlNode == NULL) throw ErrorXMLRacine("CI", fileName.str(), __FILE__, __LINE__);
    XMLElement* elementDomaine(xmlNode->FirstChildElement("physicalDomains"));
    if(elementDomaine == NULL) throw ErrorXMLElement("physicalDomains", fileName.str(), __FILE__, __LINE__);

    //3) Get defined domains
    //----------------------
    std::string nameDomaine, stateDomaine;
    XMLElement* element(elementDomaine->FirstChildElement("domain"));
    while (element != NULL)
    {
      //A)Read name domain
      //*********************
      nameDomaine = element->Attribute("name");
      if (nameDomaine == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(nameDomaine);

      //B)Read state domain
      //**********************
      stateDomaine = element->Attribute("state");
      if (stateDomaine == "") throw ErrorXMLAttribut("state", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(stateDomaine);
      //Recherche de l state associe au domain
      XMLElement* state(xmlNode->FirstChildElement("state"));
      bool trouve(false); std::string nameEtat;
      while (state != NULL)
      {
        nameEtat = state->Attribute("name");
        Tools::uppercase(nameEtat);
        if (nameEtat == stateDomaine){ trouve = true; break; }
        state = state->NextSiblingElement("state");
      }
      if (!trouve){ throw ErrorXMLEtat(stateDomaine, fileName.str(), __FILE__, __LINE__); }

      //Reading mixture state first
      Mixture* stateMixture(0);
      if (m_run->m_model->whoAmI() == "EULER") { stateMixture = new MixEuler(); }
      else if (m_run->m_model->whoAmI() == "PRESSUREVELOCITYEQ") { stateMixture = new MixPUEq(state, fileName.str()); }
      else if (m_run->m_model->whoAmI() == "VELOCITYEQ") { stateMixture = new MixUEq(state, fileName.str()); }
      else if (m_run->m_model->whoAmI() == "VELOCITYEQTOTENERGY") { stateMixture = new MixUEqTotEnergy(state, fileName.str()); }
      else if (m_run->m_model->whoAmI() == "TEMPERATUREPRESSUREVELOCITYEQ") { stateMixture = new MixPTUEq(state, fileName.str()); }
      else if (m_run->m_model->whoAmI() == "EULERHOMOGENEOUS") { stateMixture = new MixEulerHomogeneous(state, fileName.str()); }
      else if (m_run->m_model->whoAmI() == "NONLINEARSCHRODINGER") { stateMixture = new MixNonLinearSchrodinger(); }
      else if (m_run->m_model->whoAmI() == "EULERKORTEWEG") { stateMixture = new MixEulerKorteweg(); }
      else { throw ErrorXMLElement("Not valid Model-Fluid coupling", fileName.str(), __FILE__, __LINE__); }

      //Then Reading phases states
      std::vector<Phase*> statesPhases;
      std::string typeMateriau; std::string nameEOS;
      XMLElement* material(state->FirstChildElement("material"));
      int nbMateriauxState(0);
      int nbSolidState(0);
      while (material != NULL)
      {
        typeMateriau = material->Attribute("type");
        Tools::uppercase(typeMateriau);

        //Verifying that solids are setup first
        if (nbSolidState < m_run->m_numberSolids) {
          if (typeMateriau != "SOLID") { throw ErrorXMLEtat("Solid materials should be setup first", fileName.str(), __FILE__, __LINE__); }
        }

        //READING FLUID AND SOLID
        if (typeMateriau == "FLUID" || typeMateriau == "SOLID") {
          nbMateriauxState++;
          //Find EOS
          nameEOS = material->Attribute("EOS");
          int e(0);
          for (e = nbMateriauxState - 1; e < m_run->m_numberPhases; e++) {
            if (nameEOS == m_run->m_eos[e]->getName()) { break; }
          }
          if (e == m_run->m_numberPhases) { throw ErrorXMLEOSUnknown(nameEOS, fileName.str(), __FILE__, __LINE__); }
          if (e != nbMateriauxState - 1) { throw ErrorXMLEtat("Materials in used states should be ordered in reference to model.xml", fileName.str(), __FILE__, __LINE__); }
          //Fluid phase
          if (typeMateriau == "FLUID") {
            if (m_run->m_model->whoAmI() == "EULER") { statesPhases.push_back(new PhaseEuler(material, m_run->m_eos[e], fileName.str())); }
            else if (m_run->m_model->whoAmI() == "PRESSUREVELOCITYEQ") { statesPhases.push_back(new PhasePUEq(material, m_run->m_eos[e], stateMixture->getPressure(), fileName.str())); }
            else if (m_run->m_model->whoAmI() == "VELOCITYEQ") { statesPhases.push_back(new PhaseUEq(material, m_run->m_eos[e], fileName.str())); }
            else if (m_run->m_model->whoAmI() == "VELOCITYEQTOTENERGY") { statesPhases.push_back(new PhaseUEqTotEnergy(material, m_run->m_eos[e], fileName.str())); }
            else if (m_run->m_model->whoAmI() == "TEMPERATUREPRESSUREVELOCITYEQ") { statesPhases.push_back(new PhasePTUEq(material, m_run->m_eos[e], fileName.str())); }
            else if (m_run->m_model->whoAmI() == "EULERHOMOGENEOUS") { statesPhases.push_back(new PhaseEulerHomogeneous(material, m_run->m_eos[e], fileName.str())); }
            else if (m_run->m_model->whoAmI() == "EULERKORTEWEG") { statesPhases.push_back(new PhaseEulerKorteweg(material, m_run->m_eos[e], fileName.str())); }
            else { throw ErrorXMLElement("Not valid Model-Fluid coupling", fileName.str(), __FILE__, __LINE__); }
          }
          //Solid phase
          else {
            nbSolidState++;
            throw ErrorXMLElement("Not valid Model-Solid coupling", fileName.str(), __FILE__, __LINE__);
          }
        }
        //READING "NONE" (not a fluid or solid)
        else if (typeMateriau == "NONE") {
          nbMateriauxState++;
          if (m_run->m_model->whoAmI() == "NONLINEARSCHRODINGER") { statesPhases.push_back(new PhaseNonLinearSchrodinger(material, nullptr, fileName.str())); }
        }
        //UNKNOWN MATERIAL
        else { throw ErrorXMLMaterialUnknown(typeMateriau, fileName.str(), __FILE__, __LINE__); } //Cas ou le type de material n a pas ete implemente

        material = material->NextSiblingElement("material");
      }
      if (nbMateriauxState != m_run->m_numberPhases) throw ErrorXMLEtat(stateDomaine, fileName.str(), __FILE__, __LINE__);
      if (nbSolidState != m_run->m_numberSolids) throw ErrorXMLEtat(stateDomaine, fileName.str(), __FILE__, __LINE__);

      //Read des variables transportees pour l'state trouve
      std::vector<Transport> statesTransport(m_run->m_numberTransports);
      std::string nameTransport; double valueTransport(0.);
      XMLElement* elementTransport(state->FirstChildElement("transport"));
      while (elementTransport != NULL) {
        nameTransport = elementTransport->Attribute("name");
        elementTransport->QueryDoubleAttribute("value", &valueTransport);
        int e(0);
        for (e = 0; e < m_run->m_numberTransports; e++) {
          if (nameTransport == m_run->m_nameGTR[e]) { break; }
        }
        if (e != m_run->m_numberTransports) {
          statesTransport[e].setValue(valueTransport);
        }
        elementTransport = elementTransport->NextSiblingElement("transport");
      }

      //C)Reading optional physical entity
      //**********************************
      int physicalEntity(element->IntAttribute("physicalEntity")); //Default value: 0 (-1 for immersed boundaries)
      // Check that immersed boundaries are not used with an unstructured mesh
      if (m_run->m_mesh->getType() == TypeM::UNS && physicalEntity == -1) {
        throw ErrorXMLElement("Immersed boundaries not available for unstructured meshes", fileName.str(), __FILE__, __LINE__);
      }

      //C)Reading domain type
      //*********************
      std::string typeDomaine(element->Attribute("type"));
      Tools::uppercase(typeDomaine);
      std::vector<std::string> NamesParametresDomaine;
      if      (typeDomaine == "ENTIREDOMAIN"){ 
        domains.push_back(new GDEntireDomain(nameDomaine, statesPhases, stateMixture, statesTransport, physicalEntity)); 
        if(physicalEntity==-1) { solidDomains.push_back(new GDEntireDomain(nameDomaine, statesPhases, stateMixture, statesTransport, physicalEntity)); }
      }
      else if (typeDomaine == "ENTIREDOMAINWITHPARTICULARITIES"){
        domains.push_back(new GDEntireDomainWithParticularities(nameDomaine, statesPhases, stateMixture, statesTransport, physicalEntity)); 
        if(physicalEntity==-1) {solidDomains.push_back(new GDEntireDomainWithParticularities(nameDomaine, statesPhases, stateMixture, statesTransport, physicalEntity)); }
      }
      else if (typeDomaine == "HALFSPACE"){
        domains.push_back(new GDHalfSpace(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); 
        if(physicalEntity==-1) { solidDomains.push_back(new GDHalfSpace(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      }
      else if (typeDomaine == "DISC"){ 
        domains.push_back(new GDDisc(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); 
        if(physicalEntity==-1) { solidDomains.push_back(new GDDisc(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      }
      else if (typeDomaine == "ELLIPSE"){
        domains.push_back(new GDEllipse(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); 
        if(physicalEntity==-1) { solidDomains.push_back(new GDEllipse(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      }
      else if (typeDomaine == "RECTANGLE"){
        domains.push_back(new GDRectangle(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str()));
        if(physicalEntity==-1) { solidDomains.push_back(new GDRectangle(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      }
      else if (typeDomaine == "CUBOID"){
        domains.push_back(new GDCuboid(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str()));
        if(physicalEntity==-1) { solidDomains.push_back(new GDCuboid(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      }
      else if (typeDomaine == "SPHERE"){
        domains.push_back(new GDSphere(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str()));
        if(physicalEntity==-1) { solidDomains.push_back(new GDSphere(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      }
      else if (typeDomaine == "ELLIPSOID"){
        domains.push_back(new GDEllipsoid(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str()));
        if(physicalEntity==-1) { solidDomains.push_back(new GDEllipsoid(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      }
      else if (typeDomaine == "CYLINDER"){
        domains.push_back(new GDCylinder(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str()));
        if(physicalEntity==-1) { solidDomains.push_back(new GDCylinder(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }        
      }
      else{ throw ErrorXMLDomaineUnknown(typeDomaine, fileName.str(), __FILE__, __LINE__); } //Cas ou le domain n a pas ete implemente
      //Domaine suivant
      element = element->NextSiblingElement("domain");

      for (unsigned int k = 0; k < statesPhases.size(); k++) { delete statesPhases[k]; }
      delete stateMixture;
    }

    //------------------------------ CONDITIONS AUX LIMITES -------------------------------

    //4) Get root of XML document
    //---------------------------
    XMLElement* elementLimites(xmlNode->FirstChildElement("boundaryConditions"));
    if (elementLimites == NULL) throw ErrorXMLElement("boundaryConditions", fileName.str(), __FILE__, __LINE__);

    //5) Get defined boundary conditions
    //-----------------------------------
    std::string nameBoundCond, stateLimite, typeBoundCond;
    int numBoundCond;
    element = elementLimites->FirstChildElement("boundCond");
    while (element != NULL)
    {
      //A)Reading name boundCond
      //************************
      nameBoundCond = element->Attribute("name");
      if (nameBoundCond == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(nameBoundCond);

      //B)Reading type boundCond
      //************************
      typeBoundCond = element->Attribute("type");
      if (typeBoundCond == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeBoundCond);

      //C)boundCond number related to mesh
      //**********************************
      error = element->QueryIntAttribute("number", &numBoundCond);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName.str(), __FILE__, __LINE__);

      //D)Reading boundCond specific data
      //*********************************
      if (typeBoundCond == "NONREFLECTING") { boundCond.push_back(new BoundCondNonReflecting(numBoundCond)); }
      else if (typeBoundCond == "SYMMETRY") { 
        if (m_run->m_order == "FIRSTORDER") { boundCond.push_back(new BoundCondSymmetry(numBoundCond)); }
        else {
          if (m_run->m_mesh->getType() != TypeM::UNS) { // Cartesian, AMR
            boundCond.push_back(new BoundCondSymmetryO2Cartesian(numBoundCond));
          }
          else { // Unstructured
            boundCond.push_back(new BoundCondSymmetryO2NS(numBoundCond)); 
          }
        }
      }
      else if (typeBoundCond == "WALL") { 
        if (m_run->m_viscous) {
          if (m_run->m_order == "FIRSTORDER") { boundCond.push_back(new BoundCondWall(numBoundCond, element, fileName.str())); } 
          else {
            if (m_run->m_mesh->getType() != TypeM::UNS) { // Cartesian, AMR
              boundCond.push_back(new BoundCondWallO2Cartesian(numBoundCond, element, fileName.str()));
            }
            else { // Unstructured
              boundCond.push_back(new BoundCondWallO2NS(numBoundCond, element, fileName.str()));
            }
          }
        }
        else {
          if (m_run->m_order == "FIRSTORDER") { boundCond.push_back(new BoundCondSymmetry(numBoundCond)); }
          else {
            if (m_run->m_mesh->getType() != TypeM::UNS) { // Cartesian, AMR
              boundCond.push_back(new BoundCondSymmetryO2Cartesian(numBoundCond));
            }
            else { // Unstructured
              boundCond.push_back(new BoundCondSymmetryO2NS(numBoundCond)); 
            }
          }
        }
      }
      else if (typeBoundCond == "INLETTANK") { boundCond.push_back(new BoundCondInletTank(numBoundCond, element, m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_nameGTR, m_run->m_eos, fileName.str())); }
      else if (typeBoundCond == "INLETINJSTAGSTATE") { boundCond.push_back(new BoundCondInletInjStagState(numBoundCond, element, m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_nameGTR, m_run->m_eos, fileName.str())); }
      else if (typeBoundCond == "INLETINJTEMP") { boundCond.push_back(new BoundCondInletInjTemp(numBoundCond, element, m_run->m_numberPhases, fileName.str())); }
      else if (typeBoundCond == "OUTLETPRESSURE") { boundCond.push_back(new BoundCondOutletPressure(numBoundCond, element, m_run->m_numberTransports, m_run->m_nameGTR, fileName.str())); }
      else if (typeBoundCond == "OUTLETMASSFLOW") { boundCond.push_back(new BoundCondOutletMassflow(numBoundCond, element, fileName.str())); }
      else if (typeBoundCond == "NULLFLUX") { boundCond.push_back(new BoundCondNullflux(numBoundCond)); }
      else { throw ErrorXMLBoundCondUnknown(typeBoundCond, fileName.str(), __FILE__, __LINE__); } //If boundary condition not implemented
      //Next domain
      element = element->NextSiblingElement("boundCond"); 
    }
  }
  catch (ErrorXML &){ throw; }
}

//***********************************************************************

void Input::verifyCompatibilityInput(std::string testCase)
{
  // Checks parameters cross compatibility between input files
  try {
    
    // Moving Reference Frame computation restrictions
    if (m_run->m_MRF != -1) {
      if (m_run->m_mesh->getType() != TypeM::UNS) {
        std::stringstream fileName(testCase + m_nameMesh);
        throw ErrorXMLMessage("MRF is not compatible with this mesh type", fileName.str(), __FILE__, __LINE__);
      }
      if (m_run->m_addPhys.size() > 0 && m_run->m_gradient->getType() != TypeGrad::GG) {
        std::stringstream fileName(testCase + m_nameMain);
        throw ErrorXMLMessage("MRF is not compatible with this gradient method when using additionnal physics", fileName.str(), __FILE__, __LINE__);
      }
    }

    // Recording of Psat restrictions
    if (m_run->m_recordPsat == true && m_run->m_numberPhases != 2) {
      std::stringstream fileName(testCase + m_nameMain);
      throw ErrorXMLMessage("Recording of saturation pressure is possible only for a liquid/vapor fluid", fileName.str(), __FILE__, __LINE__);
    }

    // Check second-order compatibility on unstructured meshes
    if (m_run->m_order == "SECONDORDER") {
      if (m_run->m_mesh->getType() == TypeM::UNS) {
        // Check limiter compatibility
        if (m_run->m_globalLimiter->getType() == 0) { // Limiter is not Minmod or Superbee
          std::stringstream fileName(testCase + m_nameMain);
          throw ErrorXMLMessage("Selected global second-order limiter not compatible with unstructured meshes", fileName.str(), __FILE__, __LINE__);
        }
        // Check gradient compatibility
        if (m_run->m_gradient->getType() == TypeGrad::FD) {
          std::stringstream fileName(testCase + m_nameMain);
          throw ErrorXMLMessage("Selected gradient method not compatible with second-order method on unstructured meshes", fileName.str(), __FILE__, __LINE__);
        }
      }
    }
    
  }
  catch (ErrorXML &) { throw; }
}

//***********************************************************************