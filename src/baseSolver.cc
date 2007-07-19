//
// Filename     : baseSolver.cc
// Description  : Base class for solvers 
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
//                Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include <cmath>

#include "baseSolver.h"
#include "rungeKutta.h"
#include "myConfig.h"
#include "myFiles.h"
#include "myTimes.h"

BaseSolver::BaseSolver()
{
	//C_=0;
}

BaseSolver::BaseSolver(Tissue *T,std::ifstream &IN)
{
	//C_=0;
	setTissue(T);
	getInit();

  //check debugging status
	std::string debugCheck = myConfig::getValue("debug_output", 0);
	if(!debugCheck.empty()) {
		std::cerr << "Performing simulation in debug-mode\n";
		debugFlag_ = true;
		if (debugFlag_) {
			int numCopy = 10;
			cellDataCopy_.resize(numCopy);
		}
	}
	else 
		debugFlag_=false;
}

BaseSolver::~BaseSolver()
{
	
}

size_t BaseSolver::debugCount() const
{
	static size_t count = 0;
	if (debugFlag()) 
		return count++ % cellDataCopy_.size(); 
	std::cerr << "Warning  BaseSolver::debugCount() should never be" 
						<< " called when debugFlag == " << debugFlag() << "\n"; 
	exit(-1);
}

BaseSolver *BaseSolver::getSolver(Tissue *T, const std::string &file)
{
	std::istream *IN = myFiles::openFile(file);
	if (!IN) {
		std::cerr << "BaseSolver::BaseSolver() - "
							<< "Cannot open file " << file << std::endl;
		exit(EXIT_FAILURE);
	}
	std::string idValue;
	*IN >> idValue;

	BaseSolver *solver;
	if (idValue == "RK5Adaptive")
		solver = new RK5Adaptive(T,(std::ifstream &) *IN);
	else if (idValue == "RK4")
		solver = new RK4(T,(std::ifstream &) *IN);
	else {
		std::cerr << "BaseSolver::BaseSolver() - "
							<< "Unknown solver: " << idValue << std::endl;
		delete IN;
		exit(EXIT_FAILURE);
	}
	delete IN;
	return solver;
}

void BaseSolver::getInit()
{
	//
	// Resize data vectors
	//
	cellData_.resize(T_->numCell());
	cellDerivs_.resize(T_->numCell());
	wallData_.resize(T_->numWall());
	wallDerivs_.resize(T_->numWall());
	vertexData_.resize(T_->numVertex());
	vertexDerivs_.resize(T_->numVertex());
  
	//
	// Copy variable values from tissue
	//
  for( size_t i=0 ; i<T_->numCell() ; ++i ) {
    cellData_[i].resize(T_->cell(i).numVariable());
    cellDerivs_[i].resize(T_->cell(i).numVariable());
    for( size_t j=0 ; j<cellData_[i].size() ; ++j )
      cellData_[i][j]=T_->cell(i).variable(j);
  }
  
  for( size_t i=0 ; i<T_->numWall() ; ++i ) {
    wallData_[i].resize(T_->wall(i).numVariable()+1);
    wallDerivs_[i].resize(T_->wall(i).numVariable()+1);
    wallData_[i][0]=T_->wall(i).length();
    for( size_t j=1 ; j<wallData_[i].size() ; ++j )
      wallData_[i][j]=T_->wall(i).variable(j-1);
  }
  
  for( size_t i=0 ; i<T_->numVertex() ; ++i ) {
    vertexData_[i].resize(T_->vertex(i).numPosition());
    vertexDerivs_[i].resize(T_->vertex(i).numPosition());
    for( size_t j=0 ; j<vertexData_[i].size() ; ++j )
      vertexData_[i][j] = T_->vertex(i).position(j);
  }
}

void BaseSolver::readParameterFile(std::ifstream &IN)
{
  std::cerr << "BaseSolver::readParameterFile(std::ifstream &IN)\n";
  exit(-1);
}

void BaseSolver::simulate(void)
{
  std::cerr << "BaseSolver::simulate(void)\n";
  exit(-1);
}

void BaseSolver::print(std::ostream &os) 
{  
  static int tCount=0;
  static int NOld=0,okOld=0,badOld=0;
  static double tOld=0.0;
	double time=myTimes::getDiffTime();
  std::cerr << tCount << " " << t_ << " " << cellData_.size() << " " 
						<< wallData_.size() << " " << vertexData_.size() << " "
						<< numOk_ << " " << numBad_ << "  " << t_-tOld << " " 
						<< static_cast<int>(cellData_.size())-static_cast<int>(NOld) 
						<< " " << numOk_-okOld << " " << numBad_-badOld << " "
						<< time << std::endl;; 
  tOld = t_;
  NOld = cellData_.size();
  okOld = numOk_;
  badOld = numBad_;
	
  if( printFlag_==1 ) {//Print vertex and cell variables
		if( tCount==0 )
			os << numPrint_ << "\n";
	  size_t Nv = vertexData_.size(); 
		if( !Nv ) {
			os << "0 0" << std::endl << "0 0" << std::endl;
			return;
		}
		//Print the vertex positions
		size_t dimension = vertexData_[0].size();
		os << Nv << " " << dimension << std::endl;
		for( size_t i=0 ; i<Nv ; ++i ) {
			for( size_t d=0 ; d<dimension ; ++d )
				os << vertexData_[i][d] << " ";
			os << std::endl;
		}
		//os << std::endl;
		//Print the cells, first connected vertecis and then variables
		size_t Nc = cellData_.size();
		int numPrintVar=cellData_[0].size()+3;
		os << Nc << " " << numPrintVar << std::endl;
		for( size_t i=0 ; i<Nc ; ++i ) {
			size_t Ncv = T_->cell(i).numVertex(); 
			os << Ncv << " ";
			for( size_t k=0 ; k<Ncv ; ++k )
				os << T_->cell(i).vertex(k)->index() << " ";
			
			for( size_t k=0 ; k<cellData_[i].size() ; ++k )
				os << cellData_[i][k] << " ";
			os << i << " " << T_->cell(i).calculateVolume(vertexData_) << " " 
				 << T_->cell(i).numWall() << std::endl;
		}		
		os << std::endl;
  }
  else if( printFlag_==2 ) {//Print vertex and wall variables
		if( tCount==0 )
			os << numPrint_ << "\n";
		size_t Nv = vertexData_.size(); 
		if( !Nv ) {
			os << "0 0" << std::endl << "0 0" << std::endl;
			return;
		}
		//Print the vertex positions
		size_t dimension = vertexData_[0].size();
		os << Nv << " " << dimension << std::endl;
		for( size_t i=0 ; i<Nv ; ++i ) {
			for( size_t d=0 ; d<dimension ; ++d )
				os << vertexData_[i][d] << " ";
			os << std::endl;
		}
		//os << std::endl;
		// Print the walls, first connected vertecis and then variables
		size_t Nw = wallData_.size();
		//
		// Print wall variables
		//
		int numPrintVar=wallData_[0].size()+3;
		os << Nw << " " << numPrintVar << std::endl;
		for( size_t i=0 ; i<Nw ; ++i ) {
			os << "2 ";
			os << T_->wall(i).vertex1()->index() << " " 
				 << T_->wall(i).vertex2()->index() << " ";
			for( size_t k=0 ; k<wallData_[i].size() ; ++k )
				os << wallData_[i][k] << " ";
			os << i << " " << T_->wall(i).lengthFromVertexPosition(vertexData_)
				 << " " << T_->wall(i).lengthFromVertexPosition(vertexData_)-wallData_[i][0]
				 << std::endl;
		}		
		os << std::endl;
  }
  else if( printFlag_==3 ) {//Print cell variables for gnuplot
		//Print the cells, first connected vertecis and then variables
		size_t Nc = cellData_.size();
		//os << Nc << " " << numPrintVar << std::endl;
		for( size_t i=0 ; i<Nc ; ++i ) {
			for( size_t k=0 ; k<cellData_[i].size() ; ++k )
				os << cellData_[i][k] << " ";
			os << i << " " << T_->cell(i).calculateVolume(vertexData_) << " " 
				 << T_->cell(i).numWall() << std::endl;
		}		
  }
  else
    std::cerr << "BaseSolver::print() Wrong printFlag value\n";
  tCount++;
}

void BaseSolver::printInit(std::ostream &os) const
{
	assert( T_->numCell()==cellData_.size() && 
					T_->numWall()==wallData_.size() &&
					T_->numVertex()==vertexData_.size() );

  os << T_->numCell() << " " << T_->numWall() << " " << T_->numVertex() << std::endl;
	
  //Print the connectivity from walls
  for( size_t i=0 ; i<T_->numWall() ; ++i ) {
    os << i << " ";
		if( T_->wall(i).cell1()->index()<T_->numCell() )
			os << T_->wall(i).cell1()->index() << " " ;
		else
			os << "-1 ";
		if( T_->wall(i).cell2()->index()<T_->numCell() )
			os << T_->wall(i).cell2()->index() << " ";
		else
			os << "-1 ";
		os << T_->wall(i).vertex1()->index() 
			 << " " << T_->wall(i).vertex2()->index() << std::endl;
	}
  os << std::endl;
  
  //Print the vertex positions
  os << T_->numVertex() << " " << T_->vertex(0).numPosition() << std::endl;
  for( size_t i=0 ; i<T_->numVertex() ; ++i ) {
		assert( T_->vertex(i).numPosition()==vertexData_[i].size() );
    for( size_t j=0 ; j<T_->vertex(i).numPosition() ; ++j )
      os << vertexData_[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;
  
  //Print wall data
  os << T_->numWall() << " 1 " << T_->wall(0).numVariable() << std::endl;
  for( size_t i=0 ; i<T_->numWall() ; ++i ) {
		assert( wallData_[i].size() );
    for( size_t j=0 ; j<wallData_[i].size() ; ++j )
      os << wallData_[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;

  //Print cell data
  os << T_->numCell() << " " << T_->cell(0).numVariable() << std::endl;
  if( T_->cell(0).numVariable() ) {
    for( size_t i=0 ; i<T_->numCell() ; ++i ) {
			assert( cellData_[i].size() );
			for( size_t j=0 ; j<cellData_[i].size() ; ++j )
				os << cellData_[i][j] << " ";
      os << std::endl;
    }
    os << std::endl;
  }  
}

void BaseSolver::printDebug(std::ostream &os) const
{
// 	os << yCopy_.size() << "\n";
// 	size_t startElement = debugCount();
// 	for (size_t c=startElement; c<startElement+yCopy_.size(); ++c) {
// 		size_t n = c%yCopy_.size();
// 		os << yCopy_[n].size() << " " << M() << " " << c-startElement << "\n";
// 		for (size_t i=0; i<yCopy_[n].size(); i++) {
// 			for (size_t j=0; j<M(); j++) {
// 				os << yCopy_[n][i][j] << " ";
// 			}
// 			os << "\n";
// 		}
// 		os << "\n\n";
// 	}
}

