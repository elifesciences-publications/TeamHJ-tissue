//
// Filename     : rungeKutta.cc
// Description  : Different Runge Kutta solvers
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
//              : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include <cmath>
#include "rungeKutta.h"

RK5Adaptive::RK5Adaptive(Tissue *T,std::ifstream &IN)
  : BaseSolver(T,IN)
{
  readParameterFile(IN);
}

void RK5Adaptive::readParameterFile(std::ifstream &IN)
{
  IN >> startTime_;
  t_= startTime_;
  IN >> endTime_;
  
  IN >> printFlag_;
  IN >> numPrint_;
  
  IN >> h1_;
  IN >> eps_;
}

void RK5Adaptive::simulate(size_t verbose)
{ 
  // double tiny = 1e-30; // Using NR definition
  double tiny = 1e-9*eps_; // Caveat! Using new definition
  double  h, hNext, hDid;
  
  //
  //Check that h1 and endTime - startTime are > 0
  //
  if (h1_ > 0.0 && (endTime_ - startTime_) > 0.0)
    h = h1_;
  else {//either h or (endTime-startTime) <=0
    std::cerr << "RK5Adaptive::simulate() - "
	      << "Wrong time borders or time step for simulation. "
	      << "No simulation performed.\n";
    exit(-1);
  }
  
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
  // Create all vectors that will be needed here and by rkqs and rkck!
  //
  size_t Nc=T_->numCell(),Nw=T_->numWall(),Nv=T_->numVertex();
  //Used here
  std::vector< std::vector<double> > yScalC(Nc),yScalW(Nw),yScalV(Nv);
  //Used by rkqs
  std::vector< std::vector<double> > yTempC(Nc),yTempW(Nw),yTempV(Nv),
    yErrC(Nc),yErrW(Nw),yErrV(Nv);
  //used by rkck
  std::vector< std::vector<double> > ak2C(Nc),ak2W(Nw),ak2V(Nv),
    ak3C(Nc),ak3W(Nw),ak3V(Nv),
    ak4C(Nc),ak4W(Nw),ak4V(Nv),
    ak5C(Nc),ak5W(Nw),ak5V(Nv),
    ak6C(Nc),ak6W(Nw),ak6V(Nv),
    yTempRkckC(Nc),yTempRkckW(Nw),yTempRkckV(Nv);
  //Resize each vector
  size_t Ncvar=T_->cell(0).numVariable();
  for (size_t i=0; i<Nc; ++i) {
    yScalC[i].resize(cellData_[i].size());
    yErrC[i].resize(cellData_[i].size());
    yTempC[i].resize(cellData_[i].size());
    ak2C[i].resize(cellData_[i].size());
    ak3C[i].resize(cellData_[i].size());
    ak4C[i].resize(cellData_[i].size());
    ak5C[i].resize(cellData_[i].size());
    ak6C[i].resize(cellData_[i].size());
    yTempRkckC[i].resize(cellData_[i].size());
  }
  size_t Nwvar=T_->wall(0).numVariable()+1;
  for (size_t i=0; i<Nw; ++i) {
    yScalW[i].resize(wallData_[i].size());
    yErrW[i].resize(wallData_[i].size());
    yTempW[i].resize(wallData_[i].size());
    ak2W[i].resize(wallData_[i].size());
    ak3W[i].resize(wallData_[i].size());
    ak4W[i].resize(wallData_[i].size());
    ak5W[i].resize(wallData_[i].size());
    ak6W[i].resize(wallData_[i].size());
    yTempRkckW[i].resize(wallData_[i].size());
  }
  size_t Nvvar=T_->vertex(0).numPosition();
  for (size_t i=0; i<Nv; ++i) {
    yScalV[i].resize(vertexData_[i].size());
    yErrV[i].resize(vertexData_[i].size());
    yTempV[i].resize(vertexData_[i].size());
    ak2V[i].resize(vertexData_[i].size());
    ak3V[i].resize(vertexData_[i].size());
    ak4V[i].resize(vertexData_[i].size());
    ak5V[i].resize(vertexData_[i].size());
    ak6V[i].resize(vertexData_[i].size());
    yTempRkckV[i].resize(vertexData_[i].size());
  }
  
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
  double printTime = endTime_ + tiny;
  double printDeltaTime = endTime_ + 2.0 * tiny;
  if (numPrint_ <= 0) //No printing
    printFlag_ = 0;
  else if (numPrint_ == 1) { // Print last point (default)
  }
  else if (numPrint_ == 2) { //Print first/last point
    printTime = startTime_ - tiny;
  } 
  else { //Print first/last points and spread the rest uniformly
    printTime = startTime_ - tiny;
    printDeltaTime = (endTime_ - startTime_) / ((double) (numPrint_ - 1));
  }
  
  // Go
  //////////////////////////////////////////////////////////////////////
  t_ = startTime_;
  numOk_ = numBad_ = 0;
  for (unsigned int nstp = 0;; nstp++) {
    if (debugFlag()) {
      cellDataCopy_[debugCount()] = cellData_;
    } 
    // Update the derivatives
    T_->derivs(cellData_,wallData_,vertexData_,cellDerivs_,wallDerivs_,
	       vertexDerivs_);
    
    // Calculate 'scaling' for error measure
    Nc = yScalC.size();
    Nw = yScalW.size();
    Nv = yScalV.size();    
    for (size_t i=0; i<Nc; ++i) {
      Ncvar = yScalC[i].size();
      for (size_t j = 0; j<Ncvar; ++j)
        yScalC[i][j] = std::fabs(cellData_[i][j]) + 
	  std::fabs(cellDerivs_[i][j] * h) + tiny;
    }
    for (size_t i=0; i<Nw; ++i) {
      Nwvar = yScalW[i].size();
      for (size_t j = 0; j<Nwvar; ++j)
        yScalW[i][j] = std::fabs(wallData_[i][j]) + 
	  std::fabs(wallDerivs_[i][j] * h) + tiny;
    }
    for (size_t i=0; i<Nv; ++i) {
      Nvvar = yScalV[i].size();
      for (size_t j = 0; j<Nvvar; ++j)
        yScalV[i][j] = std::fabs(vertexData_[i][j]) + 
	  std::fabs(vertexDerivs_[i][j] * h) + tiny;
    }
    // Print if applicable 
    //if (printFlag_ && t_ >= printTime) {
    if (t_ >= printTime) {
      printTime += printDeltaTime;
      print();
    }
    
    // Check if step is larger than max allowed
    // max step end is min of endTime_ and printTime
    double tMin = endTime_< printTime ? endTime_ : printTime;
    if (t_+h > tMin) h = tMin - t_;
    
    // Update
    rkqs(h,hDid,hNext,yScalC,yScalW,yScalV,yTempC,yTempW,yTempV, 
	 yErrC,yErrW,yErrV,ak2C,ak2W,ak2V,ak3C,ak3W,ak3V,
	 ak4C,ak4W,ak4V,ak5C,ak5W,ak5V,ak6C,ak6W,ak6V, 
	 yTempRkckC,yTempRkckW,yTempRkckV);
    if (hDid == h) ++numOk_; else ++numBad_;
    
    //
    // Check for discrete and reaction updates
    //
    T_->updateDirection(h,cellData_,wallData_,vertexData_,cellDerivs_,
			wallDerivs_,vertexDerivs_);
    T_->updateReactions(cellData_,wallData_,vertexData_,h);
    T_->checkCompartmentChange(cellData_,wallData_,vertexData_,
			       cellDerivs_,wallDerivs_,vertexDerivs_ );
    
    // Check the tissue connectivity in each step
    T_->checkConnectivity(1);
    
    // Rescale all temporary vectors as well
    if (cellData_.size() != yScalC.size() ||
	wallData_.size() != yScalW.size() ||
	vertexData_.size() != yScalV.size())
      {
	yScalC.resize(cellData_.size(), yScalC[0]);
	yScalW.resize(wallData_.size(), yScalW[0]);
	yScalV.resize(vertexData_.size(), yScalV[0]);
	yTempC.resize(cellData_.size(), yTempC[0]);
	yTempW.resize(wallData_.size(), yTempW[0]);
	yTempV.resize(vertexData_.size(), yTempV[0]);
	yErrC.resize(cellData_.size(), yErrC[0]);
	yErrW.resize(wallData_.size(), yErrW[0]);
	yErrV.resize(vertexData_.size(), yErrV[0]);
	ak2C.resize(cellData_.size(), ak2C[0]);
	ak2W.resize(wallData_.size(), ak2W[0]);
	ak2V.resize(vertexData_.size(), ak2V[0]);
	ak3C.resize(cellData_.size(), ak3C[0]);
	ak3W.resize(wallData_.size(), ak3W[0]);
	ak3V.resize(vertexData_.size(), ak3V[0]);
	ak4C.resize(cellData_.size(), ak4C[0]);
	ak4W.resize(wallData_.size(), ak4W[0]);
	ak4V.resize(vertexData_.size(), ak4V[0]);
	ak5C.resize(cellData_.size(), ak5C[0]);
	ak5W.resize(wallData_.size(), ak5W[0]);
	ak5V.resize(vertexData_.size(), ak5V[0]);
	ak6C.resize(cellData_.size(), ak6C[0]);
	ak6W.resize(wallData_.size(), ak6W[0]);
	ak6V.resize(vertexData_.size(), ak6V[0]);
	yTempRkckC.resize(cellData_.size(), yTempRkckC[0]);
	yTempRkckW.resize(wallData_.size(), yTempRkckW[0]);
	yTempRkckV.resize(vertexData_.size(), yTempRkckV[0]);
      }
    
    // If the end t is passed return (print if applicable)
    if (t_ >= endTime_) {
      //if (printFlag_) {
      if (1) {
	// Update the derivatives
	T_->derivs(cellData_,wallData_,vertexData_,cellDerivs_,wallDerivs_,
		   vertexDerivs_);
	print();
      }
      std::cerr << "Simulation done.\n"; 
      return;
    }
    //Warn for small step sizes...
    //  if (fabs(hNext) <= hMin) {
    //        std::cerr << "Warning: Step size small (" << hNext
    //  		<< ") in rk5Adaptive::simulate at time " 
    //  		<< t << "\n"; /*exit(-1);*/ }
    h = hNext;
    //Do not take larger steps than h1
    if (h > h1_)
      h = h1_;
  }
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
void RK5Adaptive::rkqs(double hTry, double &hDid, double &hNext,
		       std::vector< std::vector<double> > &yScalC,
		       std::vector< std::vector<double> > &yScalW,
		       std::vector< std::vector<double> > &yScalV,
		       std::vector< std::vector<double> > &yTempC,
		       std::vector< std::vector<double> > &yTempW,
		       std::vector< std::vector<double> > &yTempV,
		       std::vector< std::vector<double> > &yErrC,
		       std::vector< std::vector<double> > &yErrW,
		       std::vector< std::vector<double> > &yErrV,
		       std::vector< std::vector<double> > &ak2C,
		       std::vector< std::vector<double> > &ak2W,
		       std::vector< std::vector<double> > &ak2V,
		       std::vector< std::vector<double> > &ak3C,
		       std::vector< std::vector<double> > &ak3W,
		       std::vector< std::vector<double> > &ak3V,
		       std::vector< std::vector<double> > &ak4C,
		       std::vector< std::vector<double> > &ak4W,
		       std::vector< std::vector<double> > &ak4V,
		       std::vector< std::vector<double> > &ak5C,
		       std::vector< std::vector<double> > &ak5W,
		       std::vector< std::vector<double> > &ak5V,
		       std::vector< std::vector<double> > &ak6C,
		       std::vector< std::vector<double> > &ak6W,
		       std::vector< std::vector<double> > &ak6V,
		       std::vector< std::vector<double> > &yTempRkckC,
		       std::vector< std::vector<double> > &yTempRkckW,
		       std::vector< std::vector<double> > &yTempRkckV)
{
  double errMax, h, hTemp, tNew, aux;
  h = hTry;
  for (;;) {
    rkck(h,yTempC,yTempW,yTempV,yErrC,yErrW,yErrV,ak2C,ak2W,ak2V,
	 ak3C,ak3W,ak3V,ak4C,ak4W,ak4V,ak5C,ak5W,ak5V,
	 ak6C,ak6W,ak6V,yTempRkckC,yTempRkckW,yTempRkckV);
    errMax = 0.0;
    size_t N=cellData_.size();
    for (size_t i=0; i<N; ++i) {
      size_t Nv=cellData_[i].size();
      for (size_t j=0; j<Nv; ++j) {
        aux = std::fabs(yErrC[i][j] / yScalC[i][j]);
        if (aux > errMax)
	  errMax = aux;
      }
    }
    N=wallData_.size();
    for (size_t i=0; i<N; ++i) {
      size_t Nv=wallData_[i].size();
      for (size_t j=0; j<Nv; ++j) {
        aux = std::fabs(yErrW[i][j] / yScalW[i][j]);
        if (aux > errMax)
	  errMax = aux;
      }
    }
    N=vertexData_.size();
    for (size_t i=0; i<N; ++i) {
      Nv=vertexData_[i].size();
      for (size_t j=0; j<Nv; ++j) {
        aux = std::fabs(yErrV[i][j] / yScalV[i][j]);
        if (aux > errMax)
	  errMax = aux;
      }
    }
    errMax /= eps_;
    if (errMax <= 1.0) break;
    hTemp = SAFETY * h * pow(errMax, PSHRNK);
    if (h >= 0.0)
      h = hTemp > 0.1 * h ? hTemp : 0.1 * h;
    else
      h = hTemp > 0.1 * h ? 0.1 * h : hTemp;
    tNew = t_ + h;
    if (tNew == t_) { 
      std::cerr << "Warning stepsize underflow in solverRk5Adaptive::rkqs\n"; 
      exit(-1);
    }
  }
  if (errMax > ERRCON) hNext = SAFETY * h * pow(errMax, PGROW);
  else hNext = 5.0 * h;
  t_ += (hDid = h);
  
  size_t Nc=cellData_.size();
  for (size_t i=0; i<Nc; ++i) {
    size_t Ncvar=cellData_[i].size();
    for (size_t j=0; j<Ncvar; ++j)
      cellData_[i][j] = yTempC[i][j];
  }
  size_t Nw=wallData_.size();
  for (size_t i=0; i<Nw; ++i) {
    size_t Nwvar=wallData_[i].size();
    for (size_t j=0; j<Nwvar; ++j)
      wallData_[i][j] = yTempW[i][j];
  }
  size_t Nv=vertexData_.size();
  for (size_t i=0; i<Nv; ++i) {
    size_t Nvvar=vertexData_[i].size();
    for (size_t j=0; j<Nvvar; ++j)
      vertexData_[i][j] = yTempV[i][j];
  }
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

void RK5Adaptive::
rkck(double h,
     std::vector< std::vector<double> > &yOutC,
     std::vector< std::vector<double> > &yOutW,
     std::vector< std::vector<double> > &yOutV,
     std::vector< std::vector<double> > &yErrC,
     std::vector< std::vector<double> > &yErrW,
     std::vector< std::vector<double> > &yErrV,
     std::vector< std::vector<double> > &ak2C,
     std::vector< std::vector<double> > &ak2W,
     std::vector< std::vector<double> > &ak2V,
     std::vector< std::vector<double> > &ak3C,
     std::vector< std::vector<double> > &ak3W,
     std::vector< std::vector<double> > &ak3V,
     std::vector< std::vector<double> > &ak4C,
     std::vector< std::vector<double> > &ak4W,
     std::vector< std::vector<double> > &ak4V,
     std::vector< std::vector<double> > &ak5C,
     std::vector< std::vector<double> > &ak5W,
     std::vector< std::vector<double> > &ak5V,
     std::vector< std::vector<double> > &ak6C,
     std::vector< std::vector<double> > &ak6W,
     std::vector< std::vector<double> > &ak6V,
     std::vector< std::vector<double> > &yTempRkckC, 
     std::vector< std::vector<double> > &yTempRkckW, 
     std::vector< std::vector<double> > &yTempRkckV ) 
{  
  //static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875;
  static double b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  
  size_t Nc=cellData_.size(),Nw=wallData_.size(),Nv=vertexData_.size(),
    Nwvar=wallData_[0].size(),
    Nvvar=vertexData_[0].size();
	
  for (size_t i=0; i< Nc; ++i) {
		size_t Ncvar = cellData_[i];
    for (size_t j=0; j<Ncvar; ++j)
      yTempRkckC[i][j]=cellData_[i][j]+b21*h*cellDerivs_[i][j];
	}
  for (size_t i=0; i< Nw; ++i) {
    for (size_t j=0; j<Nwvar; ++j)
      yTempRkckW[i][j]=wallData_[i][j]+b21*h*wallDerivs_[i][j];
	}
  for (size_t i=0; i< Nv; ++i) {
    for (size_t j=0; j<Nvvar; ++j)
      yTempRkckV[i][j]=vertexData_[i][j]+b21*h*vertexDerivs_[i][j];
  }
  
  T_->derivs(yTempRkckC,yTempRkckW,yTempRkckV,ak2C,ak2W,ak2V); // t + a2h
  for (size_t i=0; i< Nc; ++i) {
		size_t Ncvar = cellData_[i];
    for (size_t j=0; j<Ncvar; ++j)
      yTempRkckC[i][j]=cellData_[i][j]+h*(b31*cellDerivs_[i][j]+b32*ak2C[i][j]);
	}
  for (size_t i=0; i< Nw; ++i) {
    for (size_t j=0; j<Nwvar; ++j)
      yTempRkckW[i][j]=wallData_[i][j]+h*(b31*wallDerivs_[i][j]+b32*ak2W[i][j]);
	}
  for (size_t i=0; i< Nv; ++i) {
    for (size_t j=0; j<Nvvar; ++j)
      yTempRkckV[i][j]=vertexData_[i][j]+h*(b31*vertexDerivs_[i][j]+b32*ak2V[i][j]);
  }

  T_->derivs(yTempRkckC,yTempRkckW,yTempRkckV,ak3C,ak2W,ak3V); // t + a3h
  for (size_t i=0; i< Nc; ++i) {
		size_t Ncvar = cellData_[i];
    for (size_t j=0; j<Ncvar; ++j)
      yTempRkckC[i][j]=cellData_[i][j]+h*(b41*cellDerivs_[i][j]+b42*ak2C[i][j]+b43*ak3C[i][j]);
	}
  for (size_t i=0; i< Nw; ++i) {
    for (size_t j=0; j<Nwvar; ++j)
      yTempRkckW[i][j]=wallData_[i][j]+h*(b41*wallDerivs_[i][j]+b42*ak2W[i][j]+b43*ak3W[i][j]);
	}
  for (size_t i=0; i< Nv; ++i) {
    for (size_t j=0; j<Nvvar; ++j)
      yTempRkckV[i][j]=vertexData_[i][j]+h*(b41*vertexDerivs_[i][j]+b42*ak2V[i][j]+b43*ak3V[i][j]);
  }

  T_->derivs(yTempRkckC,yTempRkckW,yTempRkckV,ak4C,ak4W,ak4V); // t + a4 * h
  for (size_t i=0; i< Nc; ++i) {
		size_t Ncvar = cellData_[i];
		for (size_t j=0; j<Ncvar; ++j)
      yTempRkckC[i][j]=cellData_[i][j]+h*
        (b51*cellDerivs_[i][j]+b52*ak2C[i][j]+b53*ak3C[i][j]+b54*ak4C[i][j]);
	}
  for (size_t i=0; i< Nw; ++i) {
    for (size_t j=0; j<Nwvar; ++j)
      yTempRkckW[i][j]=wallData_[i][j]+h*
        (b51*wallDerivs_[i][j]+b52*ak2W[i][j]+b53*ak3W[i][j]+b54*ak4W[i][j]);
	}
  for (size_t i=0; i< Nv; ++i) {
    for (size_t j=0; j<Nvvar; ++j)
      yTempRkckV[i][j]=vertexData_[i][j]+h*
        (b51*vertexDerivs_[i][j]+b52*ak2V[i][j]+b53*ak3V[i][j]+b54*ak4V[i][j]);
  }

  T_->derivs(yTempRkckC,yTempRkckW,yTempRkckV,ak5C,ak5W,ak5V); // t + a5 * h
  for (size_t i=0; i< Nc; ++i) {
		size_t Ncvar = cellData_[i];
    for (size_t j=0; j<Ncvar; ++j)
      yTempRkckC[i][j]=cellData_[i][j]+h*
        (b61*cellDerivs_[i][j]+b62*ak2C[i][j]+b63*ak3C[i][j]+b64*ak4C[i][j]+
         b65*ak5C[i][j]);
	}
  for (size_t i=0; i< Nw; ++i) {
    for (size_t j=0; j<Nwvar; ++j)
      yTempRkckW[i][j]=wallData_[i][j]+h*
        (b61*wallDerivs_[i][j]+b62*ak2W[i][j]+b63*ak3W[i][j]+b64*ak4W[i][j]+
         b65*ak5W[i][j]);
	}
  for (size_t i=0; i< Nv; ++i) {
    for (size_t j=0; j<Nvvar; ++j)
      yTempRkckV[i][j]=vertexData_[i][j]+h*
        (b61*vertexDerivs_[i][j]+b62*ak2V[i][j]+b63*ak3V[i][j]+b64*ak4V[i][j]+
         b65*ak5V[i][j]);
  }

  T_->derivs(yTempRkckC,yTempRkckW,yTempRkckV,ak6C,ak6W,ak6V); // t + a6 * h
  for (size_t i=0; i< Nc; ++i) {
		size_t Ncvar = cellData_[i];
    for (size_t j=0; j<Ncvar; ++j)
      yOutC[i][j]=cellData_[i][j]+h*
        (c1*cellDerivs_[i][j]+c3*ak3C[i][j]+c4*ak4C[i][j]+c6*ak6C[i][j]);
	}
  for (size_t i=0; i< Nw; ++i) {
    for (size_t j=0; j<Nwvar; ++j)
      yOutW[i][j]=wallData_[i][j]+h*
        (c1*wallDerivs_[i][j]+c3*ak3W[i][j]+c4*ak4W[i][j]+c6*ak6W[i][j]);
	}
  for (size_t i=0; i< Nv; ++i) {
    for (size_t j=0; j<Nvvar; ++j)
      yOutV[i][j]=vertexData_[i][j]+h*
        (c1*vertexDerivs_[i][j]+c3*ak3V[i][j]+c4*ak4V[i][j]+c6*ak6V[i][j]);
  }

  for (size_t i=0; i< Nc; ++i) {
		size_t Ncvar = cellData_[i];
    for (size_t j=0; j<Ncvar; ++j)
      yErrC[i][j]=h*
        (dc1*cellDerivs_[i][j]+dc3*ak3C[i][j]+dc4*ak4C[i][j]+
         dc5*ak5C[i][j]+dc6*ak6C[i][j]);
	}
  for (size_t i=0; i< Nw; ++i) {
    for (size_t j=0; j<Nwvar; ++j)
      yErrW[i][j]=h*
        (dc1*wallDerivs_[i][j]+dc3*ak3W[i][j]+dc4*ak4W[i][j]+
         dc5*ak5W[i][j]+dc6*ak6W[i][j]);
	}
  for (size_t i=0; i< Nv; ++i) {
    for (size_t j=0; j<Nvvar; ++j)
      yErrV[i][j]=h*
        (dc1*vertexDerivs_[i][j]+dc3*ak3V[i][j]+dc4*ak4V[i][j]+
         dc5*ak5V[i][j]+dc6*ak6V[i][j]);
	}
}

//
// CLASS SOLVERRK4
//
RK4::RK4(Tissue *T,std::ifstream &IN)
  :BaseSolver(T,IN)
{
  readParameterFile(IN);
}

void RK4::readParameterFile(std::ifstream &IN)
{
  //Read in the needed parameters
  IN >> startTime_;
  t_=startTime_;
  IN >> endTime_;
  
  IN >> printFlag_;
  IN >> numPrint_;
  
  IN >> h_;
}

void RK4::simulate(size_t verbose) 
{
  //
  // Check that h1 and endTime-startTime are > 0
  //
  if( !(h_>0. && (endTime_-startTime_)>0.) ) {
    std::cerr << "Rk4::simulate() Wrong time borders or time step for "
	      << "simulation. No simulation performed.\n";
    return;
  }
  std::cerr << "Simulating using fourth-order Runge-Kutta\n";
  
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
  T_->initiateDirection(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
			vertexDerivs_);
  
  assert( cellData_.size() == T_->numCell() && 
	  cellData_.size()==cellDerivs_.size() );
  assert( wallData_.size() == T_->numWall() && 
	  wallData_.size()==wallDerivs_.size() );
  assert( vertexData_.size() == T_->numVertex() && 
	  vertexData_.size()==vertexDerivs_.size() );
  
  //
  // Create all vectors that will be needed by rk4()!
  //
  size_t Nc=T_->numCell(),Nw=T_->numWall(),Nv=T_->numVertex();
  std::vector< std::vector<double> > ytCell(Nc),dytCell(Nc),dymCell(Nv),
    ytWall(Nw),dytWall(Nw),dymWall(Nw),
    ytVertex(Nw),dytVertex(Nw),dymVertex(Nw);
  //Resize each vector
  size_t Ncvar=T_->cell(0).numVariable();
  for( size_t i=0 ; i<Nc ; ++i ) {
    ytCell[i].resize(Ncvar);
    dytCell[i].resize(Ncvar);
    dymCell[i].resize(Ncvar);
  }
  size_t Nwvar=T_->wall(0).numVariable()+1;
  for( size_t i=0 ; i<Nw ; ++i ) {
    ytWall[i].resize(Nwvar);
    dytWall[i].resize(Nwvar);
    dymWall[i].resize(Nwvar);
  }
  size_t Nvvar=T_->vertex(0).numPosition();
  for( size_t i=0 ; i<Nv ; ++i ) {
    ytVertex[i].resize(Nvvar);
    dytVertex[i].resize(Nvvar);
    dymVertex[i].resize(Nvvar);
  }
  
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////
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
    rk4(ytCell,ytWall,ytVertex,dytCell,dytWall,dytVertex,
	dymCell,dymWall,dymVertex);
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
    if(cellData_.size() != ytCell.size() ) {
      ytCell.resize( cellData_.size(), ytCell[0] );
      ytWall.resize( wallData_.size(), ytWall[0] );
      ytVertex.resize( vertexData_.size(), ytVertex[0] );
      dytCell.resize( cellDerivs_.size(), dytCell[0] );
      dymCell.resize( cellDerivs_.size(), dymCell[0] );
      dytWall.resize( wallDerivs_.size(), dytWall[0] );
      dymWall.resize( wallDerivs_.size(), dymWall[0] );
      dytVertex.resize( vertexDerivs_.size(), dytVertex[0] );
      dymVertex.resize( vertexDerivs_.size(), dymVertex[0] );
    }
    
    //update time variable
    if( (t_+h_)==t_ ) {
      std::cerr << "Rk4::simulate() Step size too small.";
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

void RK4::rk4(std::vector< std::vector<double> > &ytCell,
	      std::vector< std::vector<double> > &ytWall,
	      std::vector< std::vector<double> > &ytVertex,
	      std::vector< std::vector<double> > &dytCell,
	      std::vector< std::vector<double> > &dytWall,
	      std::vector< std::vector<double> > &dytVertex,
	      std::vector< std::vector<double> > &dymCell,
	      std::vector< std::vector<double> > &dymWall,
	      std::vector< std::vector<double> > &dymVertex )
{  
  double hh=0.5*h_;
  double h6=h_/6.0;
  // Take first half step
  // Is this first derivs calculation needed?
  T_->derivs(cellData_,wallData_,vertexData_,cellDerivs_,wallDerivs_,vertexDerivs_);
  for( size_t i=0 ; i<ytCell.size() ; ++i )
    for( size_t j=0 ; j<ytCell[i].size() ; ++j )
      ytCell[i][j] = cellData_[i][j] + hh*cellDerivs_[i][j];
  for( size_t i=0 ; i<ytWall.size() ; ++i )
    for( size_t j=0 ; j<ytWall[i].size() ; ++j )
      ytWall[i][j] = wallData_[i][j] + hh*wallDerivs_[i][j];
  for( size_t i=0 ; i<ytVertex.size() ; ++i )
    for( size_t j=0 ; j<ytVertex[i].size() ; ++j )
      ytVertex[i][j] = vertexData_[i][j] + hh*vertexDerivs_[i][j];
  
  // Take second half step
  T_->derivs(ytCell,ytWall,ytVertex,dytCell,dytWall,dytVertex);    
  for( size_t i=0 ; i<ytCell.size() ; ++i )
    for( size_t j=0 ; j<ytCell[i].size() ; ++j )
      ytCell[i][j] = cellData_[i][j] + hh*dytCell[i][j];
  for( size_t i=0 ; i<ytWall.size() ; ++i )
    for( size_t j=0 ; j<ytWall[i].size() ; ++j )
      ytWall[i][j] = wallData_[i][j] + hh*dytWall[i][j];
  for( size_t i=0 ; i<ytVertex.size() ; ++i )
    for( size_t j=0 ; j<ytVertex[i].size() ; ++j )
      ytVertex[i][j] = vertexData_[i][j] + hh*dytVertex[i][j];
  
  // Take temporary 'full' step
  T_->derivs(ytCell,ytWall,ytVertex,dymCell,dymWall,dymVertex);
  for( size_t i=0 ; i<cellData_.size() ; ++i )
    for( size_t j=0 ; j<cellData_[i].size() ; ++j ) {
      ytCell[i][j] = cellData_[i][j] + h_*dymCell[i][j];
      dymCell[i][j] += dytCell[i][j];
    }
  for( size_t i=0 ; i<wallData_.size() ; ++i )
    for( size_t j=0 ; j<wallData_[i].size() ; ++j ) {
      ytWall[i][j] = wallData_[i][j] + h_*dymWall[i][j];
      dymWall[i][j] += dytWall[i][j];
    }
  for( size_t i=0 ; i<vertexData_.size() ; ++i )
    for( size_t j=0 ; j<vertexData_[i].size() ; ++j ) {
      ytVertex[i][j] = vertexData_[i][j] + h_*dymVertex[i][j];
      dymVertex[i][j] = dytVertex[i][j];
    }
  // Take full step
  T_->derivs(ytCell,ytWall,ytVertex,dytCell,dytWall,dytVertex);
  for( size_t i=0 ; i<cellData_.size() ; ++i )
    for( size_t j=0 ; j<cellData_[i].size() ; ++j )
      cellData_[i][j] = cellData_[i][j] +
	h6*(cellDerivs_[i][j]+dytCell[i][j]+2.0*dymCell[i][j]);
  for( size_t i=0 ; i<wallData_.size() ; ++i )
    for( size_t j=0 ; j<wallData_[i].size() ; ++j )
      wallData_[i][j] = wallData_[i][j] +
	h6*(wallDerivs_[i][j]+dytWall[i][j]+2.0*dymWall[i][j]);
  for( size_t i=0 ; i<vertexData_.size() ; ++i )
    for( size_t j=0 ; j<vertexData_[i].size() ; ++j )
      vertexData_[i][j] = vertexData_[i][j] + 
	h6*(vertexDerivs_[i][j]+dytVertex[i][j]+2.0*dymVertex[i][j]);
}

//!Finds the maximal |dydt|/|y| for the system
double RK4::maxDerivative() {  
  std::cerr << "RK4::maxDerivative()" << std::endl;
  exit(-1);
  //   double max=0.0,val;
  //   for( size_t i=0 ; i<N() ; i++ ) {
  //     for( size_t j=0 ; j<M() ; j++ ) {
  //       if( y_[i][j]!=0 )
  // 				val = fabs( dydt_[i][j]/y_[i][j] );
  //       else
  // 				val = fabs( dydt_[i][j] );
  
  //       if( val>max ) 
  // 				max = val;
  //     }
  //   }
  //   return max;
}

