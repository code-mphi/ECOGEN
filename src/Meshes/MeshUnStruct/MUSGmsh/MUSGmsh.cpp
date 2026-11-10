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

#include "MUSGmsh.h"

//***********************************************************************

MUSGmsh::MUSGmsh(const std::string& meshFile, const std::string& meshExtension) : MeshUnStruct(meshFile, meshExtension) {}

//***********************************************************************

MUSGmsh::~MUSGmsh() {}

//***********************************************************************

std::string MUSGmsh::readVersion(const std::string& meshFile)
{
  try {
    std::string pathMeshFile(meshFile);
    std::string currentLine;
    std::ifstream mesh(pathMeshFile.c_str(), std::ios::in);
    if (!mesh) {
      throw ErrorMeshNS("mesh file not found : " + pathMeshFile, __FILE__, __LINE__);
    }

    // Exclusion of unwanted sections
    std::string fileVersion;
    while (currentLine != "$MeshFormat") {
      mesh >> currentLine; //To exclude the newline caracter
      if (mesh.eof()) {
        throw ErrorMeshNS("mesh file format is not compatible, see mesh file: " + meshFile, __FILE__, __LINE__);
      }
    }

    // Get file version
    mesh >> fileVersion; //To exclude the newline caracter

    mesh.close();
    return fileVersion;
  }
  catch (ErrorMeshNS&) {
    throw;
  }
}

//***********************************************************************

void MUSGmsh::initGeometryMonoCPU_common(TypeMeshContainer<Cell*>& cells)
{

  m_numberInnerFaces -= m_numberBoundFaces; // We remove the boundaries
  m_numberInnerFaces /= 2;                  // The internal faces are all counted twice => we restore the truth
  m_numberFacesTotal  = m_numberInnerFaces + m_numberBoundFaces;

  double volTot(0.);
  for (int i = 0; i < m_numberCellsCalcul; i++) {
    volTot += cells[i]->getElement()->getVolume();
  }
  std::cout << "Total volume : " << volTot << std::endl;

  // 1) Building connectivity table
  // -------------------------------
  // Sizing faces array
  m_faces = new FaceNS*[m_numberFacesTotal];
  int** facesBuff;
  int* sumNodesBuff; // We build a temporary array of faces to speed up the search process
  facesBuff    = new int*[m_numberFacesTotal + 1];
  sumNodesBuff = new int[m_numberFacesTotal + 1];
  // Determination of the maximum number of nodes for the faces
  int sizeFace(1); // Will be initialize at maximale size
  if (m_numberElements3D != 0) {
    if (m_numberQuadrangles != 0) {
      sizeFace = 4;
    }
    else if (m_numberTriangles != 0) {
      sizeFace = 3;
    }
    else {
      Errors::errorMessage("Issue in initGeometryMonoCPU for initialization of facesBuff array");
    }
  }
  else if (m_numberElements2D != 0) {
    sizeFace = 2;
  }
  for (int i = 0; i < m_numberFacesTotal + 1; i++) {
    facesBuff[i] = new int[sizeFace];
  }

  // Inner faces
  int indexMaxFaces{0};
  double tTemp{MPI_Wtime()};

  std::cout << "  1/Building faces ..." << std::endl;
  int printFrequency(std::max((m_numberElements - m_numberBoundFaces) / 10, 1));
  for (int i = m_numberBoundFaces; i < m_numberElements; i++) {
    if ((i - m_numberBoundFaces) % printFrequency == 0) {
      std::cout << "    " << (100 * static_cast<uint64_t>(i - m_numberBoundFaces) / (m_numberElements - m_numberBoundFaces)) << "% ... " << std::endl;
    }
    m_elements[i]->construitFaces(m_nodes, m_faces, indexMaxFaces, facesBuff, sumNodesBuff);
  }
  for (int i = 0; i < m_numberFacesTotal + 1; i++) {
    delete[] facesBuff[i];
  }
  delete[] facesBuff;
  delete[] sumNodesBuff;
  tTemp = MPI_Wtime() - tTemp;
  std::cout << "    OK in " << tTemp << " seconds" << std::endl;

  // Boundaries
  std::cout << "  2/Boundary elements attribution to boundary faces ..." << std::endl;
  tTemp          = MPI_Wtime();
  printFrequency = std::max(m_numberBoundFaces / 10, 1);
  for (int i = 0; i < m_numberBoundFaces; i++) {
    if (i % printFrequency == 0) {
      std::cout << "    " << (100 * static_cast<uint64_t>(i) / m_numberBoundFaces) << "% ... " << std::endl;
    }
    // Assigning the boundary
    m_elements[i]->attributFaceLimite(m_faces, indexMaxFaces);
  }
  tTemp = MPI_Wtime() - tTemp;
  std::cout << "    OK in " << tTemp << " seconds" << std::endl;
}
