//
// Filename     : euler.cc
// Description  : Different Euler solvers
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2012
// Revision     : $Id:$
//
#include <cmath>
#include "euler.h"

Euler::Euler(Tissue *T,std::ifstream &IN)
  : BaseSolver(T,IN)
{
  readParameterFile(IN);
}

void Euler::readParameterFile(std::ifstream &IN)
{
  IN >> startTime_;
  t_= startTime_;
  IN >> endTime_;
  
  IN >> printFlag_;
  IN >> numPrint_;
  
  IN >> h_;
}

void Euler::simulate(size_t verbose) 
{
  //
  // Check that h and endTime-startTime are > 0
  //
  if( !(h_>0. && (endTime_-startTime_)>0.) ) {
    std::cerr << "Euler::simulate() Wrong time borders or time step for "
	      << "simulation. No simulation performed." << std::endl;
    return;
  }
  std::cerr << "Simulating using explicit Euler." << std::endl;
  
  //
  // Check that sizes of permanent data is ok
  //
  if( cellData_.size() && cellData_.size() != cellDerivs_.size() ) {
    cellDerivs_.resize( cellData_.size(),cellData_[0]);
  }
  if( wallData_.size() && wallData_.size() != wallDerivs_.size() ) {
    wallDerivs_.resize( wallData_.size(),wallData_[0]);
  }
  if( vertexData_.size() && vertexData_.size() != vertexDerivs_.size() ) {
    vertexDerivs_.resize( vertexData_.size(),vertexData_[0]);
  }
  
  // Initiate reactions and direction for those where it is applicable
  T_->initiateReactions(cellData_, wallData_, vertexData_, cellDerivs_, 
			wallDerivs_, vertexDerivs_);
  if (cellData_.size()!=cellDerivs_.size())
    cellDerivs_.resize(cellData_.size(),cellDerivs_[0]);
  if (wallData_.size()!=wallDerivs_.size())
    wallDerivs_.resize(wallData_.size(),wallDerivs_[0]);
  if (vertexData_.size()!=vertexDerivs_.size())
    vertexDerivs_.resize(vertexData_.size(),vertexDerivs_[0]);
  T_->initiateDirection(cellData_, wallData_, vertexData_, cellDerivs_, 
			wallDerivs_, vertexDerivs_);
  
  assert( cellData_.size() == T_->numCell() && 
	  cellData_.size()==cellDerivs_.size() );
  assert( wallData_.size() == T_->numWall() && 
	  wallData_.size()==wallDerivs_.size() );
  assert( vertexData_.size() == T_->numVertex() && 
	  vertexData_.size()==vertexDerivs_.size() );
  
  // Initiate print times
  //
  double tiny = 1e-10;
  double printTime=endTime_+tiny;
  double printDeltaTime=endTime_+2.*tiny;
  if( numPrint_<=0 )//No printing
    printFlag_=0;
  else if( numPrint_==1 ) {//Print last point (default)
  }
  else if( numPrint_==2 ) {//Print first/last point
    printTime=startTime_-tiny;
  } 
  else {//Print first/last points and spread the rest uniformly
    printTime=startTime_-tiny;
    printDeltaTime=(endTime_-startTime_)/double(numPrint_-1);
  }
  // Go
  //
  t_=startTime_;
  numOk_ = numBad_ = 0;
  while( t_<endTime_ ) {
    if (debugFlag()) {
      cellDataCopy_[debugCount()] = cellData_;
    } 
    //Update the derivatives
    T_->derivs(cellData_,wallData_,vertexData_,cellDerivs_,wallDerivs_,
	       vertexDerivs_);
    
    //Print if applicable 
    //if( printFlag_ && t_ >= printTime ) {
    if( t_ >= printTime ) {
      printTime += printDeltaTime;
      print();
    }
    
    //Check if step is larger than max allowed
    //max step end is min of endTime_ and printTime
    //double tMin= endTime_<printTime ? endTime_ : printTime;
    //if( t_+h>tMin ) h=tMin-t_;
    
    //Update
    eulerStep();
    numOk_++;
    
    //
    // Check for discrete and reaction updates
    //
    T_->updateDirection(h_,cellData_,wallData_,vertexData_,cellDerivs_,
			wallDerivs_,vertexDerivs_);
    T_->updateReactions(cellData_,wallData_,vertexData_,h_);
    T_->checkCompartmentChange(cellData_,wallData_,vertexData_,
			       cellDerivs_,wallDerivs_,vertexDerivs_ );
    
    // Check the tissue connectivity in each step
    T_->checkConnectivity(1);
        
    //update time variable
    if( (t_+h_)==t_ ) {
      std::cerr << "Euler::simulate() Step size too small.";
      exit(-1);
    }
    t_ += h_;
  }
  //if( printFlag_ ) {
  if (1) {
    //Update the derivatives
    T_->derivs(cellData_,wallData_,vertexData_,cellDerivs_,wallDerivs_,
	       vertexDerivs_);
    print();
  }
  std::cerr << "Simulation done.\n"; 
  return;
}

void Euler::eulerStep()
{  
  // Take step
  // Is this derivs calculation needed?
  T_->derivs(cellData_,wallData_,vertexData_,cellDerivs_,wallDerivs_,vertexDerivs_);
  for( size_t i=0 ; i<cellData_.size() ; ++i )
    for( size_t j=0 ; j<cellData_[i].size() ; ++j )
      cellData_[i][j] = cellData_[i][j] + h_*cellDerivs_[i][j];
  for( size_t i=0 ; i<wallData_.size() ; ++i )
    for( size_t j=0 ; j<wallData_[i].size() ; ++j )
      wallData_[i][j] = wallData_[i][j] + h_*wallDerivs_[i][j];
  for( size_t i=0 ; i<vertexData_.size() ; ++i )
    for( size_t j=0 ; j<vertexData_[i].size() ; ++j )
      vertexData_[i][j] = vertexData_[i][j] + h_*vertexDerivs_[i][j];
}

