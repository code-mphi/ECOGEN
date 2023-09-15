#include "CellO2Cartesian.h"
#include "CellInterfaceO2Cartesian.h"

//***********************************************************************

CellO2Cartesian::CellO2Cartesian() : CellO2()
{
}

//***********************************************************************

CellO2Cartesian::CellO2Cartesian(int lvl): CellO2(lvl)
{}

//***********************************************************************

CellO2Cartesian::~CellO2Cartesian()
{
}

//***********************************************************************

void CellO2Cartesian::computeLocalSlopes(CellInterface& cellInterfaceRef, Limiter& globalLimiter, Limiter& interfaceLimiter,
  Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter,
  double& alphaCellAfterOppositeSide, double& alphaCell, double& alphaCellOtherInterfaceSide, double& epsInterface)
{
	//Mise a zero des slopes locales
	//------------------------------
  double coeff(0.), posCellInterfaceRef(0.);
	double sumCoeff(0.), sumCoeff2(0.);
	for (int k = 0; k < numberPhases; k++) {
		slopesPhasesLocal1[k]->setToZero();
		slopesPhasesLocal2[k]->setToZero();
	}
	slopesMixtureLocal1->setToZero();
	slopesMixtureLocal2->setToZero();
	for (int k = 0; k < numberTransports; k++) {
		slopesTransportLocal1[k] = 0.;
		slopesTransportLocal2[k] = 0.;
	}

	//Loop on the cell interfaces for the determination of the slopes on each side of the cell
	//----------------------------------------------------------------------------------------
  int phase0(0);
	for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
    //Calcul de la slope a gauche et a droite de la cell (AMR) en se basant sur celle de reference (cell interface a gauche ou a droite, unknown)
    if (m_cellInterfaces[b]->getSlopesPhase(0) != 0) { //Cell interface de type CellInterfaceO2 ou BoundCondWallO2
      if (m_cellInterfaces[b] == &cellInterfaceRef) {
				for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), 1.); }
				slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), 1.);
				for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += m_cellInterfaces[b]->getSlopesTransport(k)->getValue(); }
        sumCoeff += 1.;
      }
      else {
				if (!m_cellInterfaces[b]->getSplit()) {
				//Produit scalar des normals avec celle de reference
					coeff = std::fabs(m_cellInterfaces[b]->getFace()->getNormal().scalar(cellInterfaceRef.getFace()->getNormal()));
					if (coeff > 1.e-6) {
						//Face majoritement selon X
						if (std::fabs(cellInterfaceRef.getFace()->getNormal().getX()) > 0.5) {
							posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getX();
							//Cote cellInterfaceRef
							if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getX()) <= std::fabs(posCellInterfaceRef - m_element->getPosition().getX())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sumCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sumCoeff2 += coeff;
                if (m_cellInterfaces[b]->getCellLeft() == this) {
                  if (m_cellInterfaces[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellRight()->getPhase(phase0)->getAlpha(); } //Cell interface of type CellInterface/O2 (inner)
                  else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellLeft()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellLeft()->getPhase(phase0)->getAlpha(); }
							}
						}
						//Face majoritement selon Y
						else if (std::fabs(cellInterfaceRef.getFace()->getNormal().getY()) > 0.5) {
							posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getY();
							//Cote cellInterfaceRef
							if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getY()) <= std::fabs(posCellInterfaceRef - m_element->getPosition().getY())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sumCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sumCoeff2 += coeff;
                if (m_cellInterfaces[b]->getCellLeft() == this) {
                  if (m_cellInterfaces[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellRight()->getPhase(phase0)->getAlpha(); } //Cell interface of type CellInterface/O2 (inner)
                  else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellLeft()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellLeft()->getPhase(phase0)->getAlpha(); }
							}
						}
						//Face majoritement selon Z
						else {
							posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getZ();
							//Cote cellInterfaceRef
							if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getZ()) <= std::fabs(posCellInterfaceRef - m_element->getPosition().getZ())) {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sumCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
								slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
								for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
								sumCoeff2 += coeff;
                if (m_cellInterfaces[b]->getCellLeft() == this) {
                  if (m_cellInterfaces[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellRight()->getPhase(phase0)->getAlpha(); } //Cell interface of type CellInterface/O2 (inner)
                  else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellLeft()->getPhase(phase0)->getAlpha(); }
                }
                else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellLeft()->getPhase(phase0)->getAlpha(); }
							}
						}
					}
				}
      }
    }
	} //End loop on cell interfaces

	//Normalisation des slopes
	//------------------------
	if (sumCoeff > 1.e-8) {
		for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->divide(sumCoeff);	}
		slopesMixtureLocal1->divide(sumCoeff);
		for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] /= sumCoeff; }
	}
	if (sumCoeff2 > 1.e-8) {
		for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->divide(sumCoeff2);	}
		slopesMixtureLocal2->divide(sumCoeff2);
		for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] /= sumCoeff2; }
    alphaCellAfterOppositeSide /= sumCoeff2;
	}

	//Limitations des slopes
	//----------------------
  //Detection of the interface location
  if ((alphaCell >= epsInterface) && (alphaCell <= 1. - epsInterface) && ((alphaCellOtherInterfaceSide - alphaCell)*(alphaCell - alphaCellAfterOppositeSide) >= 1.e-8)) {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *slopesPhasesLocal2[k], interfaceLimiter, interfaceVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *slopesMixtureLocal2, interfaceLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = interfaceVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], slopesTransportLocal2[k]);
    }
  }
  else {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *slopesPhasesLocal2[k], globalLimiter, globalVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *slopesMixtureLocal2, globalLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = globalVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], slopesTransportLocal2[k]);
    }
  }
}

//***********************************************************************

void CellO2Cartesian::computeLocalSlopesLimite(CellInterface& cellInterfaceRef, Limiter& globalLimiter, Limiter& interfaceLimiter,
  Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, double& epsInterface)
{
  //Solution pour multiD Cartesian (peut etre une ebauche pour le NS, a voir...)

  //Set to zero the local slopes
  //----------------------------
  double coeff(0.), posCellInterfaceRef(0.);
  double sumCoeff2(0.);
  for (int k = 0; k < numberPhases; k++) {
    slopesPhasesLocal1[k]->setToZero();
    slopesPhasesLocal2[k]->setToZero();
  }
  slopesMixtureLocal1->setToZero();
  slopesMixtureLocal2->setToZero();
  for (int k = 0; k < numberTransports; k++) {
    slopesTransportLocal1[k] = 0.;
    slopesTransportLocal2[k] = 0.;
  }

  //Get the slope on the CL side
  //----------------------------
  for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*cellInterfaceRef.getSlopesPhase(k), 1.); }
  slopesMixtureLocal1->multiplyAndAdd(*cellInterfaceRef.getSlopesMixture(), 1.);
  for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += cellInterfaceRef.getSlopesTransport(k)->getValue(); }

  //Loop on the cell interfaces for the determination of the slopes on the opposite side of the BC
  //----------------------------------------------------------------------------------------------
  for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
    //Calcul de la slope a gauche et a droite de la cell (AMR) en se basant sur celle de reference (cell interface a gauche ou a droite, unknown)
    if (m_cellInterfaces[b]->getSlopesPhase(0) != 0) { //Cell interface de type CellInterface/O2
      if (m_cellInterfaces[b] != &cellInterfaceRef) {
        if (!m_cellInterfaces[b]->getSplit()) {
          //Produit scalar des normals avec celle de reference
          coeff = std::fabs(m_cellInterfaces[b]->getFace()->getNormal().scalar(cellInterfaceRef.getFace()->getNormal()));
          if (coeff > 1.e-6) {
            //Face majoritement selon X
            if (std::fabs(cellInterfaceRef.getFace()->getNormal().getX()) > 0.5) {
              posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getX();
              //Autre cote
              if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getX()) >= std::fabs(posCellInterfaceRef - m_element->getPosition().getX())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
                sumCoeff2 += coeff;
              }
            }
            //Face majoritement selon Y
            else if (std::fabs(cellInterfaceRef.getFace()->getNormal().getY()) > 0.5) {
              posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getY();
              //Autre cote
              if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getY()) >= std::fabs(posCellInterfaceRef - m_element->getPosition().getY())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
                sumCoeff2 += coeff;
              }
            }
            //Face majoritement selon Z
            else {
              posCellInterfaceRef = cellInterfaceRef.getFace()->getPos().getZ();
              //Autre cote
              if (std::fabs(posCellInterfaceRef - m_cellInterfaces[b]->getFace()->getPos().getZ()) >= std::fabs(posCellInterfaceRef - m_element->getPosition().getZ())) {
                for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), coeff); }
                slopesMixtureLocal2->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), coeff);
                for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] += coeff*(m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
                sumCoeff2 += coeff;
              }
            }
          }
        }
      }
    }
  } //End loop on cell interfaces

  //Normalisation de la slope
  //-------------------------
  if (sumCoeff2 > 1.e-8) {
    for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal2[k]->divide(sumCoeff2); }
    slopesMixtureLocal2->divide(sumCoeff2);
    for (int k = 0; k < numberTransports; k++) { slopesTransportLocal2[k] /= sumCoeff2; }
  }

  //Limitations des slopes
  //----------------------
  //Detection of the interface location
  int phase0(0);
  if ((m_vecPhases[phase0]->getAlpha() >= epsInterface) && (m_vecPhases[phase0]->getAlpha() <= 1. - epsInterface)) {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *slopesPhasesLocal2[k], interfaceLimiter, interfaceVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *slopesMixtureLocal2, interfaceLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = interfaceVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], slopesTransportLocal2[k]);
    }
  }
  else {
    for (int k = 0; k < numberPhases; k++) {
      slopesPhasesLocal1[k]->limitSlopes(*slopesPhasesLocal1[k], *slopesPhasesLocal2[k], globalLimiter, globalVolumeFractionLimiter);
    }
    slopesMixtureLocal1->limitSlopes(*slopesMixtureLocal1, *slopesMixtureLocal2, globalLimiter);
    for (int k = 0; k < numberTransports; k++) {
      slopesTransportLocal1[k] = globalVolumeFractionLimiter.limiteSlope(slopesTransportLocal1[k], slopesTransportLocal2[k]);
    }
  }
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

void CellO2Cartesian::createChildCell(const int& lvl)
{
  m_childrenCells.push_back(new CellO2Cartesian(lvl + 1));
}

//****************************************************************************
//********************** Methode Ordre 2 Parallele ***************************
//****************************************************************************

void CellO2Cartesian::fillBufferSlopes(double* buffer, int& counter, const int& lvl, const int& neighbour) const
{
	if (m_lvl == lvl) {
    std::vector<CellInterface*> cellInterfacesWithNeighboringGhostCell;

    for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
      if (m_cellInterfaces[b]->whoAmI() == 0) { //Cell interface of type CellInterface/O2 (inner)
        if (m_cellInterfaces[b]->getLvl() == m_lvl) {
          if (this == m_cellInterfaces[b]->getCellLeft()) {
            if (m_cellInterfaces[b]->getCellRight()->isCellGhost()) {
              if (m_cellInterfaces[b]->getCellRight()->getRankOfNeighborCPU() == neighbour) {
                cellInterfacesWithNeighboringGhostCell.push_back(m_cellInterfaces[b]);
              }
            }
          }
          else {
            if (m_cellInterfaces[b]->getCellLeft()->isCellGhost()) {
              if (m_cellInterfaces[b]->getCellLeft()->getRankOfNeighborCPU() == neighbour) {
                cellInterfacesWithNeighboringGhostCell.push_back(m_cellInterfaces[b]);
              }
            }
          }
        }
      }
    }

    double epsilon(1.e-08), sumCoeff(0.), scalarDiff(0.), alphaCellAfterOppositeSide(0.);
    Coord coordBuffer(0.);
    int phase0(0);
    for (unsigned int b2 = 0; b2 < cellInterfacesWithNeighboringGhostCell.size(); b2++) {
  		//Reset local slope to send
  		//-------------------------
  		sumCoeff = 0.;
  		for (int k = 0; k < numberPhases; k++) {
  			slopesPhasesLocal1[k]->setToZero();
  		}
  		slopesMixtureLocal1->setToZero();
  		for (int k = 0; k < numberTransports; k++) {
  			slopesTransportLocal1[k] = 0.;
  		}

      //Loop over cell interfaces to determine the slope to send
      //--------------------------------------------------------
      alphaCellAfterOppositeSide = 0.;
      int slopeIndex=-1;
      for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
        if (m_cellInterfaces[b] != cellInterfacesWithNeighboringGhostCell[b2]) {
          if (m_cellInterfaces[b]->getSlopesPhase(0) != 0) { //Cell interface de type CellInterfaceO2 ou BoundCondWallO2
            if (!m_cellInterfaces[b]->getSplit()) {
              coordBuffer = m_cellInterfaces[b]->getFace()->getNormal().abs() - cellInterfacesWithNeighboringGhostCell[b2]->getFace()->getNormal();
              if (coordBuffer.norm() < epsilon) { //Face in the same direction than the reference face
                scalarDiff = m_cellInterfaces[b]->getFace()->getPos().scalar(m_cellInterfaces[b]->getFace()->getNormal()) -
                              cellInterfacesWithNeighboringGhostCell[b2]->getFace()->getPos().scalar(cellInterfacesWithNeighboringGhostCell[b2]->getFace()->getNormal());
                if (std::fabs(scalarDiff) > epsilon) { //Face on the opposite side of the cell in regards to the reference face
                  for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesPhase(k), 1.); }
                  slopesMixtureLocal1->multiplyAndAdd(*m_cellInterfaces[b]->getSlopesMixture(), 1.);
                  for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] += (m_cellInterfaces[b]->getSlopesTransport(k)->getValue()); }
                  sumCoeff += 1.;
                  if (m_cellInterfaces[b]->getCellLeft() == this) {
                    if (m_cellInterfaces[b]->whoAmI() == 0) { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellRight()->getPhase(phase0)->getAlpha(); } //Cell interface of type CellInterface/O2 (inner)
                    else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellLeft()->getPhase(phase0)->getAlpha(); }
                  }
                  else { alphaCellAfterOppositeSide += m_cellInterfaces[b]->getCellLeft()->getPhase(phase0)->getAlpha(); }
                  coordBuffer=scalarDiff*cellInterfacesWithNeighboringGhostCell[b2]->getFace()->getNormal();
                  if      (coordBuffer.getX() < -epsilon) slopeIndex=0;
                  else if (coordBuffer.getX() >  epsilon) slopeIndex=1;
                  else if (coordBuffer.getY() < -epsilon) slopeIndex=2;
                  else if (coordBuffer.getY() >  epsilon) slopeIndex=3;
                  else if (coordBuffer.getZ() < -epsilon) slopeIndex=4;
                  else if (coordBuffer.getZ() >  epsilon) slopeIndex=5;
                }
              }
            }
          }
        }
      }

  		//Normalization of the slope
  		//--------------------------
  		if (sumCoeff > 1.e-8) {
  			for (int k = 0; k < numberPhases; k++) { slopesPhasesLocal1[k]->divide(sumCoeff); }
  			slopesMixtureLocal1->divide(sumCoeff);
  			for (int k = 0; k < numberTransports; k++) { slopesTransportLocal1[k] /= sumCoeff; }
        alphaCellAfterOppositeSide /= sumCoeff;
  		}

  		//Fill buffer to send
  		//-------------------
  		for (int k = 0; k < numberPhases; k++) {
  			slopesPhasesLocal1[k]->fillBufferSlopes(buffer, counter);
  		}
  		slopesMixtureLocal1->fillBufferSlopes(buffer, counter);
  		for (int k = 0; k < numberTransports; k++) {
  			buffer[++counter] = slopesTransportLocal1[k];
  		}
      buffer[++counter] = alphaCellAfterOppositeSide;
      buffer[++counter] = static_cast<double>(slopeIndex);
    }
	}

	else {
    for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
      if (m_childrenCells[i]->hasNeighboringGhostCellOfCPUneighbour(neighbour)) {
        m_childrenCells[i]->fillBufferSlopes(buffer, counter, lvl, neighbour);
      }
    }
	}
}

//***********************************************************************