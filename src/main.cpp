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

#include "Config.h"
#include "Run.h"
#include "Errors.h"
#include "libTierces/tinyxml2.h"

#include <time.h>

using namespace tinyxml2;

void displayHeader();

//***********************************************************************

int main(int argc, char* argv[])
{
  Run* run(0);

  //Error redirection
  std::ofstream Out("Error.out");
  std::streambuf* strBackup(std::cerr.rdbuf());
  std::cerr.rdbuf(Out.rdbuf());
  int errorCode(0);

  //Parallel initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankCpu);
  MPI_Comm_size(MPI_COMM_WORLD, &Ncpu);

  if(rankCpu == 0) displayHeader();
  MPI_Barrier(MPI_COMM_WORLD);

  //Parsing of the XML file by the tinyxml2 library
  //-----------------------------------------------
  try {
    //Working on specific folder
    if (argc > 1) {
      config.setWorkFolder(std::string(argv[1]));
    }
    std::stringstream fileName;
    fileName << config.getWorkFolder();
    fileName << "ECOGEN.xml";

    tinyxml2::XMLDocument xmlEcogen;
    XMLError error(xmlEcogen.LoadFile(fileName.str().c_str())); //The file is parse here
    if (error != XML_NO_ERROR) throw ErrorXML(fileName.str().c_str(), __FILE__, __LINE__);
    //Find the root of the XML folder
    XMLNode *xmlNode = xmlEcogen.FirstChildElement("ecogen");

    //Loop on the test cases to execute
    //---------------------------------
    int numTestCase(0);
    XMLElement* elementTestCase = xmlNode->FirstChildElement("testCase");
    while (elementTestCase != NULL) {
      try {
        XMLNode* xmlNode2 = elementTestCase->FirstChild();
        if (xmlNode2 == NULL) throw ErrorXMLElement("testCase", fileName.str(), __FILE__, __LINE__);
        XMLText* xmlText = xmlNode2->ToText();
        if (xmlText == NULL) throw ErrorXMLElement("testCase", fileName.str(), __FILE__, __LINE__);
        std::stringstream testLocation;
        testLocation << config.getWorkFolder();
        testLocation << xmlText->Value();

        //1) Creation of the test case
        numTestCase++;
        run = new Run(testLocation.str(), numTestCase);
        MPI_Barrier(MPI_COMM_WORLD);
        if (rankCpu == 0) {
          std::cout << "************************************************************" << std::endl;
          std::cout << "           EXECUTION OF THE TEST CASE NUMBER: " << numTestCase << std::endl;
          std::cout << "************************************************************" << std::endl;
          std::cout << "T" << numTestCase << " | Test case: " << testLocation.str() << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //2) Execution of the test case
        run->initialize();
        run->solver();
        run->finalize();
        //3) Removal of the test case
        delete run;
      }
    
      //Input exceptions (XML + mesh file)
      //----------------------------------
      catch (ErrorInput &e) {
        errorCode = e.getErrorCode();
        if (rankCpu == 0) std::cout << e.infoError() << std::endl;
        if (run) {
          run->finalize();
          delete run;
        }
      }
      //General or runtime exceptions
      //-----------------------------
      catch (ErrorECOGEN &e) {
        errorCode = e.getErrorCode();
       // if(rankCpu==0) std::cerr << "T" << numTestCase << " | " << e.infoError() << std::endl;
        if (run) {
          run->finalize();
          delete run;
        }
      }
      elementTestCase = elementTestCase->NextSiblingElement("testCase");
    }//End of the loop on test cases
  }
  catch (ErrorXML& e) {
    errorCode = e.getErrorCode();
    if (rankCpu == 0) std::cout << e.infoError() << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  std::cerr.rdbuf(strBackup);
  Out.close();
  return errorCode;
}

//***********************************************************************

void displayHeader()
{
  std::cout << "************************************************************" << std::endl;
  std::cout << std::endl;
  std::cout << "    ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. " << std::endl;
  std::cout << "    | .-'   .' .')   / .-. )  .' .'     | .-'    |  \\| | " << std::endl;
  std::cout << "    | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | " << std::endl;
  std::cout << "    | .-'   \\  \\     | | | |  \\  \\ ( _) | .-'    | |\\  | " << std::endl;
  std::cout << "    |  `--.  \\  `-.  \\ `-' /   \\  `-) ) |  `--.  | | |)| " << std::endl;
  std::cout << "    /( __.'   \\____\\  )---'    )\\____/  /( __.'  /(  (_) " << std::endl;
  std::cout << "   (__)              (_)      (__)     (__)     (__)     " << std::endl;
  std::cout << "   Official webSite: https://code-mphi.github.io/ECOGEN/" << std::endl;
  std::cout << std::endl;
}

//***********************************************************************