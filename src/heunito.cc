//
// Filename     : heunito.cc
// Description  : Heun numerical solver in the Ito interpretation (Carrillo et al. 2003 PRE)
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : December 2014
// Revision     : $Id:$
//
#include <cmath>
#include "heunito.h"
#include "myRandom.h"

HeunIto::HeunIto(Tissue *T,std::ifstream &IN)
  :BaseSolver(T,IN)
{
  readParameterFile(IN);
}

void HeunIto::readParameterFile(std::ifstream &IN)
{
  //Read in the needed parameters
  IN >> startTime_;
  t_=startTime_;
  IN >> endTime_;
  
  IN >> printFlag_; // output format
  IN >> numPrint_;  // number of time points printed to output
  
  IN >> h_;   // time step
  IN >> vol_; // effective volume
}

void HeunIto::simulate(size_t verbose) 
{
  //
  // Check that h1 and endTime-startTime are > 0
  //
  if( !(h_>0. && (endTime_-startTime_)>0.) ) {
    std::cerr << "HeunIto::simulate() Wrong time borders or time step for "
	      << "simulation. No simulation performed.\n";
    return;
  }
  std::cerr << "Simulating using an HeunIto solver\n";

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
  //
  // Create all vectors that will be needed by the HeunIto algorithm
  //
  size_t Nc=T_->numCell(),Nw=T_->numWall(),Nv=T_->numVertex();
  DataMatrix sdydtCell(Nc),stCell(Nc),y1Cell(Nc),dydt2Cell(Nc),
    sdydtWall(Nw),stWall(Nw),y1Wall(Nw),dydt2Wall(Nw),
    sdydtVertex(Nv),stVertex(Nv),y1Vertex(Nv),dydt2Vertex(Nv);
  //Resize each vector
  //size_t Ncvar=T_->cell(0).numVariable();
  for( size_t i=0 ; i<Nc ; ++i ) {
    sdydtCell[i].resize(cellData_[i].size());
    stCell[i].resize(cellData_[i].size());
    y1Cell[i].resize(cellData_[i].size());
    dydt2Cell[i].resize(cellData_[i].size());
  }
  //size_t Nwvar=T_->wall(0).numVariable()+1;
  for( size_t i=0 ; i<Nw ; ++i ) {
    sdydtWall[i].resize(wallData_[i].size());
    stWall[i].resize(wallData_[i].size());
    y1Wall[i].resize(wallData_[i].size());
    dydt2Wall[i].resize(wallData_[i].size());
  }
  //size_t Nvvar=T_->vertex(0).numPosition();
  for( size_t i=0 ; i<Nv ; ++i ) {
    sdydtVertex[i].resize(vertexData_[i].size());
    stVertex[i].resize(vertexData_[i].size());
    y1Vertex[i].resize(vertexData_[i].size());
    dydt2Vertex[i].resize(vertexData_[i].size());
  }
  
  // Initiate print times
  //
  double tiny = 1e-10;
  double printTime=endTime_+tiny;
  double printDeltaTime=endTime_+2.*tiny;
	int doPrint=1;
  if( numPrint_<=0 )//No printing
    doPrint=0;
  else if( numPrint_==1 ) {//Print last point (default)
  }
  else if( numPrint_==2 ) {//Print first/last point
    printTime=startTime_-tiny;
  } 
  else {//Print first/last points and spread the rest uniformly
    printTime=startTime_-tiny;
    printDeltaTime=(endTime_-startTime_)/double(numPrint_-1);
  }
  //
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
    if( doPrint && t_ >= printTime ) {
      printTime += printDeltaTime;
      print();
    }

    //Check if step is larger than max allowed
    //max step end is min of endTime_ and printTime
    //double tMin= endTime_<printTime ? endTime_ : printTime;
    //if( t_+h>tMin ) h=tMin-t_;

    //Update
    heunito(sdydtCell,sdydtWall,sdydtVertex,stCell,stWall,stVertex,
	    y1Cell,y1Wall,y1Vertex,dydt2Cell,dydt2Wall,dydt2Vertex);
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
    
    // Resize temporary containers as well
    if(cellData_.size() != sdydtCell.size() ) {
      sdydtCell.resize( cellData_.size(), sdydtCell[0] );
      stCell.resize( cellDerivs_.size(), stCell[0] );
      y1Cell.resize( cellDerivs_.size(), y1Cell[0] );
      dydt2Cell.resize( cellDerivs_.size(), dydt2Cell[0] );
      sdydtWall.resize( wallData_.size(), sdydtWall[0] );
      stWall.resize( wallDerivs_.size(), stWall[0] );
      y1Wall.resize( wallDerivs_.size(), y1Wall[0] );
      dydt2Wall.resize( wallDerivs_.size(), dydt2Wall[0] );
      sdydtVertex.resize( vertexData_.size(), sdydtVertex[0] );
      stVertex.resize( vertexDerivs_.size(), stVertex[0] );
      y1Vertex.resize( vertexDerivs_.size(), y1Vertex[0] );
      dydt2Vertex.resize( vertexDerivs_.size(), dydt2Vertex[0] );
    }
    
    //update time variable
    if( (t_+h_)==t_ ) {
      std::cerr << "HeunIto::simulate() Step size too small.";
      exit(-1);
    }
    t_ += h_;
  }
  if( doPrint ) {
    //Update the derivatives
    T_->derivs(cellData_,wallData_,vertexData_,cellDerivs_,wallDerivs_,
	       vertexDerivs_);
    print();
  }
  std::cerr << "Simulation done.\n"; 
  return;
}

void HeunIto::heunito(DataMatrix &sdydtCell,
		      DataMatrix &sdydtWall,
		      DataMatrix &sdydtVertex,
		      DataMatrix &stCell,
		      DataMatrix &stWall,
		      DataMatrix &stVertex,
		      DataMatrix &y1Cell,
		      DataMatrix &y1Wall,
		      DataMatrix &y1Vertex,
		      DataMatrix &dydt2Cell,
		      DataMatrix &dydt2Wall,
		      DataMatrix &dydt2Vertex)
{ 
  double hh=0.5*h_;
  
  T_->derivsWithAbs(cellData_,wallData_,vertexData_,cellDerivs_,wallDerivs_,vertexDerivs_,
		    sdydtCell,sdydtWall,sdydtVertex);// first step
  DataMatrix randCell( sdydtCell.size() ),randWall( sdydtWall.size() ),randVertex( sdydtVertex.size() );
  for(size_t i=0 ; i<sdydtCell.size() ; ++i ) {
    randCell[i].resize( sdydtCell[i].size() );
    for( size_t j=0 ; j<sdydtCell[i].size() ; ++j ) {      
      randCell[i][j] = myRandom::Grand();
      stCell[i][j] = sqrt(sdydtCell[i][j]*h_/vol_);
      y1Cell[i][j] = cellData_[i][j]+h_*cellDerivs_[i][j]+stCell[i][j]*randCell[i][j];
    }
  }
  for(size_t i=0 ; i<sdydtWall.size() ; ++i ) {
    randWall[i].resize( sdydtWall[i].size() );
    for( size_t j=0 ; j<sdydtWall[i].size() ; ++j ) {      
      randWall[i][j] = myRandom::Grand();
      stWall[i][j] = sqrt(sdydtWall[i][j]*h_/vol_);
      y1Wall[i][j] = wallData_[i][j]+h_*wallDerivs_[i][j]+stWall[i][j]*randWall[i][j];
    }
  }
  for(size_t i=0 ; i<sdydtVertex.size() ; ++i ) {
    randVertex[i].resize( sdydtVertex[i].size() );
    for( size_t j=0 ; j<sdydtVertex[i].size() ; ++j ) {      
      randVertex[i][j] = myRandom::Grand();
      stVertex[i][j] = sqrt(sdydtVertex[i][j]*h_/vol_);
      y1Vertex[i][j] = vertexData_[i][j]+h_*vertexDerivs_[i][j]+stVertex[i][j]*randVertex[i][j];
    }
  }
  
  T_->derivs(y1Cell,y1Wall,y1Vertex,dydt2Cell,dydt2Wall,dydt2Vertex);//second step
  for(size_t i=0 ; i<sdydtCell.size() ; ++i ) {
    for( size_t j=0 ; j<sdydtCell[i].size() ; ++j ) {      
      cellData_[i][j] = y1Cell[i][j]+hh*(cellDerivs_[i][j]+dydt2Cell[i][j])+stCell[i][j]*randCell[i][j];      
      if(cellData_[i][j]<0.0 )// Setting and absortive barrier at 0
	cellData_[i][j]=0.0; 
    }
  }
  for(size_t i=0 ; i<sdydtWall.size() ; ++i ) {
    for( size_t j=0 ; j<sdydtWall[i].size() ; ++j ) {      
      wallData_[i][j] = y1Wall[i][j]+hh*(wallDerivs_[i][j]+dydt2Wall[i][j])+stWall[i][j]*randWall[i][j];      
      if(wallData_[i][j]<0.0 )// Setting and absortive barrier at 0
	wallData_[i][j]=0.0; 
    }
  }
  for(size_t i=0 ; i<sdydtVertex.size() ; ++i ) {
    for( size_t j=0 ; j<sdydtVertex[i].size() ; ++j ) {      
      vertexData_[i][j] = y1Vertex[i][j]+hh*(vertexDerivs_[i][j]+dydt2Vertex[i][j])+
	stVertex[i][j]*randVertex[i][j];      
      if(vertexData_[i][j]<0.0 )// Setting and absortive barrier at 0
	vertexData_[i][j]=0.0; 
    }
  }
}

