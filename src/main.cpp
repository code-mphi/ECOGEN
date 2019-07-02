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

//! \file      main.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      December 7 2017

/*! \mainpage ECOGEN API documentation
 *
 * This official API documentation for ECOGEN.
 *
 * To begin with ECOGEN, visit the official webSite : https://code-mphi.github.io/ECOGEN/
 */

#include "Run.h"
#include "Errors.h"
#include "libTierces/tinyxml2.h"

using namespace std;
using namespace tinyxml2;

void displayHeader();

//***********************************************************************

int main(int argc, char* argv[])
{
  Run* run(0);

  //Parallel initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankCpu);
  MPI_Comm_size(MPI_COMM_WORLD, &Ncpu);

  if(rankCpu == 0) displayHeader();
  MPI_Barrier(MPI_COMM_WORLD);

  //Parsing of the XML file by the tinyxml2 library
  //-----------------------------------------------
  stringstream fileName("ECOGEN.xml");
  XMLDocument xmlEcogen;
  XMLError error(xmlEcogen.LoadFile(fileName.str().c_str())); //The file is parse here
  //if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);
  //Find the root of the XML folder
  XMLNode *xmlNode = xmlEcogen.FirstChildElement("ecogen");
  //if (xmlNode == NULL) throw ErrorXMLRacine("ecogen", fileName.str(), __FILE__, __LINE__);

  //Loop on the test cases to execute
  //---------------------------------
  int numTestCase(0);
  XMLElement *elementTestCase = xmlNode->FirstChildElement("testCase");
  while (elementTestCase != NULL) {
    try {
      XMLNode* xmlNode2 = elementTestCase->FirstChild();
      if (xmlNode2 == NULL) throw ErrorXMLElement("testCase", fileName.str(), __FILE__, __LINE__);
      XMLText* xmlText = xmlNode2->ToText();
      if (xmlText == NULL) throw ErrorXMLElement("testCase", fileName.str(), __FILE__, __LINE__);

      //1) Creation of the test case
      numTestCase++;
      run = new Run(xmlText->Value(), numTestCase);
      MPI_Barrier(MPI_COMM_WORLD);
      if (rankCpu == 0) {
        cout << "************************************************************" << endl;
        cout << "          EXECUTION OF THE TEST CASE NUMBER : " << numTestCase << endl;
        cout << "************************************************************" << endl;
        cout << "T" << numTestCase << " | Test case : " << xmlText->Value() << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      //2) Execution of the test case
      //if (rankCpu == 0) { cout << "wait launch test" << endl;  system("pause"); }
      //MPI_Barrier(MPI_COMM_WORLD);
      run->initialize(argc, argv);
      run->solver();
      //if (rankCpu == 0) { cout << "wait end test" << endl; system("pause"); }
      //MPI_Barrier(MPI_COMM_WORLD);
      run->finalize();
      //3) Removal of the test case
      delete run;
      //if (rankCpu == 0) { cout << "test finished" << endl; system("pause"); }
      //MPI_Barrier(MPI_COMM_WORLD);
    }
    //Gestion of the exceptions
    //-------------------------
    catch (ErrorXML &) {
      if (!run) {
        run->finalize();
        delete run;
      }
    }
    catch (ErrorECOGEN &e) {
      if(rankCpu==0) cerr << e.infoError() << endl;
      //cerr << "T" << numTestCase << " | " << e.infoError() << endl;
      if (!run) {
        run->finalize();
        delete run;
      }
    }
    elementTestCase = elementTestCase->NextSiblingElement("testCase");
  }//End of the loop on test cases
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

//***********************************************************************

void displayHeader()
{
  cout << "************************************************************" << endl;
  cout << "    ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. " << endl;
  cout << "    | .-'   .' .')   / .-. )  .' .'     | .-'    |  \\| | " << endl;
  cout << "    | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | " << endl;
  cout << "    | .-'   \\  \\     | | | |  \\  \\ ( _) | .-'    | |\\  | " << endl;
  cout << "    |  `--.  \\  `-.  \\ `-' /   \\  `-) ) |  `--.  | | |)| " << endl;
  cout << "    /( __.'   \\____\\  )---'    )\\____/  /( __.'  /(  (_) " << endl;
  cout << "   (__)              (_)      (__)     (__)     (__)     " << endl;
  cout << "************************************************************" << endl;
}

//***********************************************************************