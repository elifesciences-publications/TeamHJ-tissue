//
// Filename     : directionReaction.cc
// Description  : Classes describing some reaction updates related to directions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : May 2008
// Revision     : $Id:$
//

#include"directionReaction.h"
#include"tissue.h"
#include"baseReaction.h"
#include"myMath.h"
#include<cmath>


ContinousMTDirection::ContinousMTDirection(std::vector<double> &paraValue,
					   std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "ContinousMTDirection::ContinousMTDirection() " 
	      << "Uses one parameter: k_rate" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 2 || indValue[0].size() != 1 ||  indValue[1].size() != 1) {
    std::cerr << "ContinousMTDirection::ContinousMTDirection() " << std::endl
	      << "First level gives target direction index (input)." << std::endl
	      << "Second level gives real direction index." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("ContinousMTDirection");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  
  setParameterId(tmp);
}

void ContinousMTDirection::derivs(Tissue &T,
				  DataMatrix &cellData,
				  DataMatrix &wallData,
				  DataMatrix &vertexData,
				  DataMatrix &cellDerivs,
				  DataMatrix &wallDerivs,
				  DataMatrix &vertexDerivs)
{ size_t dimension=vertexData[0].size();
  if (dimension!=2) {
    std::cerr << "ContinuosMTDirection::derivs() Only implemented for two dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }
  size_t target = variableIndex(0, 0);
  size_t real = variableIndex(1, 0);
  double k_rate = parameter(0);
  

  for (size_t n = 0; n < T.numCell(); ++n) {

    Cell cell = T.cell(n);
    size_t index = cell.index();
    
    double x = cellData[index][real + 0];
    double y = cellData[index][real + 1];
    
    double sigma = std::atan2(y, x);
    
    while (sigma > 0.5 * myMath::pi() || sigma <= -0.5 * myMath::pi()) {
      if (sigma > 0.5 * myMath::pi()) {
        sigma -= myMath::pi();
      }
      if (sigma <= -0.5 * myMath::pi()) {
        sigma += myMath::pi();
      }
    }
    
    double dx = cellData[index][target + 0];
    double dy = cellData[index][target + 1];
    
    double dsigma = std::atan2(dy, dx);
    
    while (dsigma > 0.5 * myMath::pi() || dsigma <= -0.5 * myMath::pi()) {
      if (dsigma > 0.5 * myMath::pi()) {
        dsigma -= myMath::pi();
      }
      if (dsigma <= -0.5 * myMath::pi()) {
        dsigma += myMath::pi();
      }
    }
    
    double angle = dsigma - sigma;
    
    while (angle > myMath::pi() || angle <= -myMath::pi()) {
      if (angle > myMath::pi()) {
        angle -= 2.0 * myMath::pi();
      }
      if (angle <= -myMath::pi()) {
        angle += 2.0 * myMath::pi();
      }
    }
    
    double speed = k_rate * std::abs(angle) / (0.25 * myMath::pi() + std::abs(angle));
    
    speed *= (angle >= 0) ? +1 : -1;
    
    cellDerivs[index][real + 0] += -y * speed;
    cellDerivs[index][real + 1] += x * speed;
  }
}


ContinousMTDirection3d::ContinousMTDirection3d(std::vector<double> &paraValue,
					   std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "ContinousMTDirection3d::ContinousMTDirection3d() " 
	      << "Uses one parameter: k_rate" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 2 || indValue[0].size() != 1 ||  indValue[1].size() != 1) {
    std::cerr << "ContinousMTDirection3d::ContinousMTDirection() " << std::endl
	      << "First level gives target direction index (input)." << std::endl
	      << "Second level gives real direction index." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("ContinousMTDirection3d");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  
  setParameterId(tmp);
}

void ContinousMTDirection3d::derivs(Tissue &T,
				  DataMatrix &cellData,
				  DataMatrix &wallData,
				  DataMatrix &vertexData,
				  DataMatrix &cellDerivs,
				  DataMatrix &wallDerivs,
				  DataMatrix &vertexDerivs)
{
  size_t target = variableIndex(0, 0);
  size_t real = variableIndex(1, 0);
  double k_rate = parameter(0);
  

  for (size_t n = 0; n < T.numCell(); ++n) {

    Cell cell = T.cell(n);
    size_t index = cell.index();
    
    double tmp=std::sqrt( cellData[index][real + 0]*cellData[index][real + 0]+
                          cellData[index][real + 1]*cellData[index][real + 1]+
                          cellData[index][real + 2]*cellData[index][real + 2]);
    if (tmp==0) {
      cellData[index][real + 0]=1;
      tmp=1;
    }
    cellData[index][real + 0]/=tmp;
    cellData[index][real + 1]/=tmp;
    cellData[index][real + 2]/=tmp;

    double x = cellData[index][real + 0];
    double y = cellData[index][real + 1];
    double z = cellData[index][real + 2];

    tmp=std::sqrt( cellData[index][target + 0]*cellData[index][target + 0]+
                   cellData[index][target + 1]*cellData[index][target + 1]+
                   cellData[index][target + 2]*cellData[index][target + 2]);
    if (tmp==0) {
      cellData[index][target + 0]=1;
      tmp=1;
    }
    cellData[index][target + 0]/=tmp;
    cellData[index][target + 1]/=tmp;
    cellData[index][target + 2]/=tmp;
    
    double tx = cellData[index][target + 0];
    double ty = cellData[index][target + 1];
    double tz = cellData[index][target + 2];
    
    double inner=tx*x+ty*y+tz*z;
    if (inner<0){ 
      cellData[index][real + 0] *=-1;
      cellData[index][real + 1] *=-1;
      cellData[index][real + 2] *=-1;
      x *=-1;
      y *=-1;
      z *=-1;     
    }

    
    double dx =tx-x;
    double dy =ty-y;
    double dz =tz-z;
    
    
    cellDerivs[index][real + 0] +=k_rate * dx;
    cellDerivs[index][real + 1] +=k_rate * dy;
    cellDerivs[index][real + 2] +=k_rate * dz;
  }

}

UpdateMTDirection::UpdateMTDirection(std::vector<double> &paraValue,
				     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "UpdateMTDirection::UpdateMTDirection() " 
	      << "Uses one parameter: k_rate" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 2 || indValue[0].size() != 1 ||  indValue[1].size() != 1) {
    std::cerr << "UpdateMTDirection::UpdateMTDirection() " << std::endl
	      << "First level gives target direction index (input)." << std::endl
	      << "Second level gives real direction index." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("UpdateMTDirection");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  
  setParameterId(tmp);
}

void UpdateMTDirection::initiate(Tissue &T,
				 DataMatrix &cellData,
				 DataMatrix &wallData,
				 DataMatrix &vertexData) 
{
  size_t numCell=cellData.size();
  size_t dimension=vertexData[0].size();
  size_t inIndex=variableIndex(0,0);
  size_t outIndex=variableIndex(1,0);
  for (size_t i=0; i<numCell; ++i)
    for (size_t d=0; d<dimension; ++d)
      cellData[i][outIndex+d] = cellData[i][inIndex+d];
}

void UpdateMTDirection::derivs(Tissue &T,
			       DataMatrix &cellData,
			       DataMatrix &wallData,
			       DataMatrix &vertexData,
			       DataMatrix &cellDerivs,
			       DataMatrix &wallDerivs,
			       DataMatrix &vertexDerivs ) {}

void UpdateMTDirection::update(Tissue &T,
			       DataMatrix &cellData,
			       DataMatrix &wallData,
			       DataMatrix &vertexData,
			       double h) 
{
  size_t numCell=cellData.size();
  size_t dimension=vertexData[0].size();
  size_t inIndex=variableIndex(0,0);
  size_t outIndex=variableIndex(1,0);
  if (parameter(0)==0.0)
    return;
  for (size_t i=0; i<numCell; ++i) {
    for (size_t d=0; d<dimension; ++d)
      cellData[i][outIndex+d] += parameter(0)*h*(cellData[i][inIndex+d]-cellData[i][outIndex+d]);
    // Normalize
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += cellData[i][outIndex+d]*cellData[i][outIndex+d];
    norm = 1.0/std::sqrt(norm);
    for (size_t d=0; d<dimension; ++d)
      cellData[i][outIndex+d] *= norm;
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

UpdateMTDirectionConcenHill::UpdateMTDirectionConcenHill(std::vector<double> &paraValue,
				     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 3) {
    std::cerr << "UpdateMTDirectionConcenHill::UpdateMTDirectionConcenHill() " 
	      << "Uses three parameters: k_rate, k_Hill and n_Hill " << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 3 || indValue[0].size() != 1 ||  indValue[1].size() != 1 || indValue[2].size() != 1) {
    std::cerr << "UpdateMTDirectionConcenHill::UpdateMTDirectionConcenHill() " << std::endl
	      << "First level gives target direction index (input)." << std::endl
	      << "Second level gives real direction index." << std::endl
              << "Third level gives concentration(anisotropy) index." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("UpdateMTDirectionConcenHill");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  tmp[1] = "k_Hill";
  tmp[2] = "n_Hill";
  
  setParameterId(tmp);
}

void UpdateMTDirectionConcenHill::initiate(Tissue &T,
				 DataMatrix &cellData,
				 DataMatrix &wallData,
				 DataMatrix &vertexData) 
{
  // size_t numCell=cellData.size();
  // size_t dimension=vertexData[0].size();
  // size_t inIndex=variableIndex(0,0);
  // size_t outIndex=variableIndex(1,0);
  // for (size_t i=0; i<numCell; ++i)
  //   for (size_t d=0; d<dimension; ++d)
  //     cellData[i][outIndex+d] = cellData[i][inIndex+d];
}

void UpdateMTDirectionConcenHill::derivs(Tissue &T,
			       DataMatrix &cellData,
			       DataMatrix &wallData,
			       DataMatrix &vertexData,
			       DataMatrix &cellDerivs,
			       DataMatrix &wallDerivs,
			       DataMatrix &vertexDerivs ) {

  size_t numCell=cellData.size();
  size_t dimension=vertexData[0].size();
  size_t inIndex=variableIndex(0,0);
  size_t outIndex=variableIndex(1,0);
  size_t concIndex=variableIndex(2,0);
  if (parameter(0)==0.0)
    return;
  for (size_t i=0; i<numCell; ++i) {
    double kh=parameter(1);
    double nh=parameter(2);    
    double ConcFactor=(std::pow(cellData[i][concIndex],nh)/(std::pow(((1-cellData[i][concIndex])*kh),nh)+std::pow(cellData[i][concIndex],nh)));
    // std::cerr<<i<<" "<< ConcFactor<<std::endl;
    if (cellData[i][inIndex]*cellData[i][outIndex]+
        cellData[i][inIndex+1]*cellData[i][outIndex+1]+
        cellData[i][inIndex+2]*cellData[i][outIndex+2]<0){
     for (size_t d=0; d<dimension; ++d)
      cellData[i][outIndex+d] *=-1;
    }
    for (size_t d=0; d<dimension; ++d)
      cellData[i][outIndex+d] += parameter(0)*ConcFactor*(cellData[i][inIndex+d]-cellData[i][outIndex+d]);
    // Normalize
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += cellData[i][outIndex+d]*cellData[i][outIndex+d];
    norm = 1.0/std::sqrt(norm);
    for (size_t d=0; d<dimension; ++d)
      cellData[i][outIndex+d] *= norm;
  }

}

void UpdateMTDirectionConcenHill::update(Tissue &T,
			       DataMatrix &cellData,
			       DataMatrix &wallData,
			       DataMatrix &vertexData,
			       double h) 
{
  // size_t numCell=cellData.size();
  // size_t dimension=vertexData[0].size();
  // size_t inIndex=variableIndex(0,0);
  // size_t outIndex=variableIndex(1,0);
  // size_t concIndex=variableIndex(2,0);
  // if (parameter(0)==0.0)
  //   return;
  // for (size_t i=0; i<numCell; ++i) {
  //   double kh=parameter(1);
  //   double nh=parameter(2);    
  //   double ConcFactor=(std::pow(cellData[i][concIndex],nh)/(std::pow(((1-cellData[i][concIndex])*kh),nh)+std::pow(cellData[i][concIndex],nh)));
  //   // std::cerr<<i<<" "<< ConcFactor<<std::endl;
  //   if (cellData[i][inIndex]*cellData[i][outIndex]+
  //       cellData[i][inIndex+1]*cellData[i][outIndex+1]+
  //       cellData[i][inIndex+2]*cellData[i][outIndex+2]<0){
  //    for (size_t d=0; d<dimension; ++d)
  //     cellData[i][outIndex+d] *=-1;
  //   }
  //   for (size_t d=0; d<dimension; ++d)
  //     cellData[i][outIndex+d] += parameter(0)*ConcFactor*(cellData[i][inIndex+d]-cellData[i][outIndex+d]);
  //   // Normalize
  //   double norm=0.0;
  //   for (size_t d=0; d<dimension; ++d)
  //     norm += cellData[i][outIndex+d]*cellData[i][outIndex+d];
  //   norm = 1.0/std::sqrt(norm);
  //   for (size_t d=0; d<dimension; ++d)
  //     cellData[i][outIndex+d] *= norm;
  // }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

RotatingDirection::RotatingDirection(std::vector<double> &paraValue,
				     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "RotatingDirection::RotatingDirection() " 
	      << "Uses one parameter: k_rate" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "RotatingDirection::RotatingDirection() \n"
	      << "First level gives cell MT index.\n";
    std::exit(EXIT_FAILURE);
  }
  
  setId("RotatingDirection");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  
  setParameterId(tmp);
}

void RotatingDirection::derivs(Tissue &T,
			       DataMatrix &cellData,
			       DataMatrix &wallData,
			       DataMatrix &vertexData,
			       DataMatrix &cellDerivs,
			       DataMatrix &wallDerivs,
			       DataMatrix &vertexDerivs)
{
  const size_t xIndex = variableIndex(0, 0) + 0;
  const size_t yIndex = variableIndex(0, 0) + 1;
  
  const double k = parameter(0);
  
  for (size_t cellIndex = 0; cellIndex < T.numCell(); ++cellIndex)
    {
      const double x = cellData[cellIndex][xIndex];
      const double y = cellData[cellIndex][yIndex];
      
      cellDerivs[cellIndex][xIndex] -= k * y;
      cellDerivs[cellIndex][yIndex] += k* x;
    }
}
