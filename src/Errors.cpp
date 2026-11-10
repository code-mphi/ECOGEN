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

#include "Errors.h"
#include "Run.h"
#include <mpi.h>

std::vector<Errors> errors;
std::vector<Errors> warnings;

//***********************************************************************

Errors::Errors() : m_message("not specified"), m_state(0), m_line(0), m_value(0.) {}

//***********************************************************************

Errors::Errors(const std::string& message, const char* sourceFile, int lineNumber) :
  m_message(message), m_file(sourceFile), m_line(lineNumber), m_value(0.)
{
  m_state = 1;
}

//***********************************************************************

Errors::~Errors() {}

//***********************************************************************

void Errors::errorMessage(const std::string& message)
{
  std::stringstream numCPU;
  numCPU << rankCpu;

  std::cout << std::endl << "-------------------------------------------------------" << std::endl;
  std::cout << "ERROR on CPU " + numCPU.str() + ": " << message.c_str() << std::endl;
  std::cout << "=> Fix error and restart run" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
}

//***********************************************************************

void Errors::errorMessage(const std::string& message, const char* sourceFile, int lineNumber)
{
  std::stringstream numCPU;
  numCPU << rankCpu;

  std::cout << std::endl << "-------------------------------------------------------" << std::endl;
  std::cout << "ERROR on CPU " + numCPU.str() << std::endl;
  std::cout << " - file: " << sourceFile << " line : " << lineNumber << std::endl;
  std::cout << " - message: " << message.c_str() << std::endl;
  std::cout << "=> Fix error and restart run" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
}

//***********************************************************************

void Errors::errorMessage(const std::string& message, double value)
{
  std::stringstream numCPU;
  numCPU << rankCpu;

  std::cout << std::endl << "-------------------------------------------------------" << std::endl;
  std::cout << "ERROR on CPU " + numCPU.str() + " : " << message.c_str() << " " << value << std::endl;
  std::cout << "=> Fix error and restart run" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
}

//***********************************************************************

void Errors::setError(const std::string& message, const char* sourceFile, int lineNumber)
{
  m_message = message;
  m_state   = 1;
  m_file    = sourceFile;
  m_line    = lineNumber;
}

//***********************************************************************

void Errors::setError(const std::string& message, const double value)
{
  m_message = message;
  m_state   = 1;
  m_value   = value;
}

//***********************************************************************

int Errors::getState() { return m_state; }

//***********************************************************************
std::string Errors::errorToString(const int& num, const int& ErrorType, const int iteration) const
{

  std::string errormsg;
  std::string msg = "\n-------------------------------------------------------\n";

  if (ErrorType == ERROR) {
    msg      += "        ERROR " + std::to_string(num) + " REQUIRING PROGRAM SHUTDOWN\n";
    errormsg  = "=> Fix error and restart run\n";
  }
  else if (ErrorType == WARNING) {
    msg      += "        WARNING " + std::to_string(num) + " REQUIRING ATTENTION\n";
    errormsg  = "";
  }

  if (iteration != -1) {
    msg += " - iteration: " + std::to_string(iteration) + "\n";
  }

  msg += " - CPU number: " + std::to_string(rankCpu) + "\n";
  msg += " - file: " + m_file + " line: " + std::to_string(m_line) + "\n";
  msg += " - message: " + m_message + "\n";
  msg += errormsg;
  msg += "-------------------------------------------------------\n";

  return msg;
}

//***********************************************************************

void Errors::displayError(const int& num) const
{
  std::string msg{errorToString(num, ERROR)};
  std::cout << msg;
}

//***********************************************************************

constexpr int Errors::defaultInt;
constexpr int Errors::defaultIntNeg;
constexpr double Errors::defaultDouble;
const std::string Errors::defaultString = "NA";

//***********************************************************************
