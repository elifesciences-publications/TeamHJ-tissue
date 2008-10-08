//
// Filename     : baseSolver.cc
// Description  : Base class for solvers 
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
//                Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include <cmath>
#include <set>
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

void BaseSolver::setTissueVariables()
{
	//
	// Check size of data vectors
	//
	if (cellData_.size() != T_->numCell() ||
			wallData_.size() != T_->numWall() ||
			vertexData_.size() != T_->numVertex()) {
		std::cerr << "BaseSolver::setTissueVariables wrong size of data:" << std::endl
							<< "Data\tTissue\tInternal" << std::endl
							<< "Cell\t" << T_->numCell() << "\t" << cellData_.size() << std::endl
							<< "Wall\t" << T_->numWall() << "\t" << wallData_.size() << std::endl
							<< "Vertex\t" << T_->numVertex() << "\t" << vertexData_.size() << std::endl;
		exit(-1);
  }
	//
	// Copy variable values to tissue
	//
  for (size_t i=0; i<T_->numCell(); ++i) {
    for (size_t j=0; j<cellData_[i].size(); ++j)
      T_->cell(i).setVariable(j,cellData_[i][j]);
  }
  
  for (size_t i=0; i<T_->numWall(); ++i) {
    T_->wall(i).setLength(wallData_[i][0]);
		T_->wall(i).setVariable(wallData_[i]);
  }
  
  for (size_t i=0; i<T_->numVertex(); ++i) {
		T_->vertex(i).setPosition(vertexData_[i]);
  }
}

void BaseSolver::readParameterFile(std::ifstream &IN)
{
  std::cerr << "BaseSolver::readParameterFile(std::ifstream &IN)\n";
  exit(-1);
}

void BaseSolver::simulate(size_t verbose)
{
	if (verbose)
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
	///
	/// Print vertex, cell, and wall variables
  ///
	if( printFlag_==0 ) {
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
		// Print wall variables, first the two connected vertices and then the variables
		numPrintVar=wallData_[0].size()+4;
		size_t Nw = wallData_.size();
		os << Nw << " " << numPrintVar << std::endl;
		for( size_t i=0 ; i<Nw ; ++i ) {
			//os << "2 ";
			os << T_->wall(i).vertex1()->index() << " " 
				 << T_->wall(i).vertex2()->index() << " ";
			for( size_t k=0 ; k<wallData_[i].size() ; ++k )
				os << wallData_[i][k] << " ";
			double distance = T_->wall(i).lengthFromVertexPosition(vertexData_);
			os << i << " " << distance
				 << " " << distance-wallData_[i][0] << " " << (distance-wallData_[i][0])/wallData_[i][0]
				 << std::endl;
		}		
		os << std::endl;
  }
	///
	/// Print vertex and cell variables
	///
  else if( printFlag_==1 ) {
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
	///
	/// Print vertex and wall variables
	///
  else if( printFlag_==2 ) {
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
		int numPrintVar=wallData_[0].size()+4;
		os << Nw << " " << numPrintVar << std::endl;
		for( size_t i=0 ; i<Nw ; ++i ) {
			os << "2 ";
			os << T_->wall(i).vertex1()->index() << " " 
				 << T_->wall(i).vertex2()->index() << " ";
			for( size_t k=0 ; k<wallData_[i].size() ; ++k )
				os << wallData_[i][k] << " ";
			double distance = T_->wall(i).lengthFromVertexPosition(vertexData_);
			os << i << " " << distance
				 << " " << distance-wallData_[i][0] << " " << (distance-wallData_[i][0])/wallData_[i][0]
				 << std::endl;
		}		
		os << std::endl;
  }
	///
	/// Print cell variables for gnuplot
	///
  else if( printFlag_==3 ) {
		//Print the cells, first connected vertecis and then variables
		size_t Nc = cellData_.size();
		//os << Nc << " " << numPrintVar << std::endl;
		for( size_t i=0 ; i<Nc ; ++i ) {
			os << "0 " << i << " " << t_ << " ";
			for( size_t k=0 ; k<cellData_[i].size() ; ++k )
				os << cellData_[i][k] << " ";
			os << i << " " << T_->cell(i).calculateVolume(vertexData_) << " " 
				 << T_->cell(i).numWall() << std::endl;
		}		
		size_t Nw = wallData_.size();
		for( size_t i=0 ; i<Nw ; ++i ) {
			os << "1 " << i << " " << t_ << " ";
			for( size_t k=0 ; k<wallData_[i].size() ; ++k )
				os << wallData_[i][k] << " ";
			os << i << " " << T_->wall(i).lengthFromVertexPosition(vertexData_)
				 << " " << T_->wall(i).lengthFromVertexPosition(vertexData_)-wallData_[i][0]
				 << std::endl;
		}				
		os << std::endl;
  }

  else if (printFlag_ == 97) {
	  for (size_t i = 0; i < T_->numCell(); ++i) {
		  Cell &cell = T_->cell(i);
		  
		  if (cell.isNeighbor(T_->background())) {
			  continue;
		  }
		  
		  double length = 0.0;
		  
		  for (size_t j = 0; j < cell.numWall(); ++j) {
			  Wall *wall = cell.wall(j);
			  
			  length += wall->lengthFromVertexPosition(vertexData_);
		  }
		  
		  double area = cell.calculateVolume(vertexData_, 0);
		  
		  std::cout << std::pow(length, 2.0) / area << "\n";
	  }
  }

  else if (printFlag_ == 98) {
	  std::list<int> neighbors;

	  for (size_t i = 0; i < T_->numCell(); ++i) {
		  Cell &cell = T_->cell(i);

		  if (cell.isNeighbor(T_->background())) {
			  continue;
		  } else {
			  std::set<Cell *> candidates;

			  for (size_t j = 0; j < cell.numWall(); ++j)
			  {
				  Wall *wall = cell.wall(j);
				  
				  if (wall->cell1() != &cell) {
					  candidates.insert(wall->cell1());
				  }
				  
				  if (wall->cell2() != &cell) {
					  candidates.insert(wall->cell2());
					  
				  }
			  }
			  
			  if (!candidates.empty()) {
				  neighbors.push_back(candidates.size());
			  }
		  }
	  }
	  
	  if (neighbors.empty()) {
		  std::cout << "> 0\n";
		  return;
	  } 
	  
	  std::vector<int> histogram(*(std::max_element(neighbors.begin(), neighbors.end())) + 1, 0);
	  
	  for (std::list<int>::const_iterator i = neighbors.begin(), e = neighbors.end(); i != e; ++i) {
		  ++histogram[*i];
	  }
	  
	  double sum = 0.0;

	  for (std::vector<int>::const_iterator i = histogram.begin(), e = histogram.end(); i != e; ++i) {
		  sum += *i;
	  }
	  
	  std::cout << "> " <<  histogram.size() << "\n";
	  for (size_t i = 0; i < histogram.size(); ++i) {
		  std::cout << i << " " << (double) histogram[i] / sum << "\n";
	  }
  }

  else if (printFlag_ == 99) {
	  size_t dimensions = vertexData_[0].size();
	  
	  for (size_t i = 0; i < T_->numVertex(); ++i) {
		  Vertex *vertex = T_->vertexP(i);
		  
		  std::vector<double> vertexPosition(dimensions);
		  
		  std::cout << t_ << " ";
		  
		  for (size_t j = 0; j < dimensions; ++j) {
			  std::cout << (vertexPosition[j] = vertexData_[vertex->index()][j]) << " ";
		  }

// 		  std::cout << "\n";

// 		  std::cout << t_ << " ";

		  std::vector<double> stressDirection = vertex->getStressDirection();

		  double A = 0.0;
		  
		  for (size_t j = 0; j < stressDirection.size(); ++j) {
			  A += (stressDirection[j] * stressDirection[j]);
		  }

		  A = std::sqrt(A);

		  for (size_t j = 0; j < stressDirection.size(); ++j) {
			  std::cout << (stressDirection[j] / A) << " ";
		  }

		  std::cout << "\n";
		  
	  }
	  std::cout << "\n";
  }

	///
	/// For printing pin1 also in membranes
	///
	else if (printFlag_==4) {
		if( tCount==0 )
			os << numPrint_ << "\n";
		size_t Nv = vertexData_.size(); 
		if( !Nv ) {
			os << "0 0" << std::endl << "0 0" << std::endl;
			return;
		}
		//Print the vertex positions
		size_t dimension = vertexData_[0].size();
		os << T_->numVertex() << " " << dimension << std::endl;
		for( size_t i=0 ; i<Nv ; ++i ) {
			for( size_t d=0 ; d<dimension ; ++d )
				os << vertexData_[i][d] << " ";
			os << std::endl;
		}
		os << std::endl;
		//Print the cells, first connected vertecis and then variables
		size_t Nc = cellData_.size();
		// For membrane PIN1 printing (version3)
		//////////////////////////////////////////////////////////////////////
		os << Nc << std::endl;
		size_t pinI=8,xI=9;
		//size_t auxinI=4,mI=7;
		std::vector<double> parameter(8);
		parameter[0]=0.1;
		parameter[3]=0.01;
		for( size_t i=0 ; i<Nc ; ++i ) {
			size_t Ncv = T_->cell(i).numVertex(); 
			os << Ncv << " ";
			for( size_t k=0 ; k<Ncv ; ++k )
				os << T_->cell(i).vertex(k)->index() << " ";
			
			//pin polarization
			//////////////////////////////////////////////////////////////////////
			size_t numWalls=T_->cell(i).numWall();
			
			//Polarization coefficient normalization constant
			double sum=0.0;
			std::vector<double> Pij(numWalls);
			for( size_t n=0 ; n<numWalls ; ++n ) {
				if( T_->cell(i).wall(n)->cell1() != T_->background() &&
						T_->cell(i).wall(n)->cell2() != T_->background() ) { 
					size_t neighI;
					if( T_->cell(i).wall(n)->cell1()->index()==i )
						neighI = T_->cell(i).wall(n)->cell2()->index();
					else
						neighI = T_->cell(i).wall(n)->cell1()->index();
					//double powX = std::pow(cellData_[ neighI ][ xI ],parameter[5));
					//double Cij = powX/(std::pow(parameter[4),parameter[5))+powX);
					sum += Pij[n] = cellData_[ neighI ][ xI ];
					//sum += Pij[n] = (1.0-parameter[2]) + 
					//parameter[2]*cellData_[ neighI ][xI];
				}
				else 
					sum += Pij[n] = parameter[0];
			}
			//sum /= numWalls;//For adjusting for different num neigh
			
			if( sum >= 0.0 )
				os << parameter[3]*cellData_[i][pinI] / sum << " ";
			else
				os << "0.0 ";
			
			for( size_t n=0 ; n<numWalls ; ++n ) {
				double pol=0.0;
				if( sum != 0.0 )
					pol = cellData_[i][pinI] * Pij[n] / sum;
				os << pol << " ";
			}
			os << std::endl;
		}
    //////////////////////////////////////////////////////////////////////
    //End Pij printing version3 
	}
	///
	/// For printing pin1 also in membranes
	///
	else if (printFlag_==5) {
		if( tCount==0 )
			os << numPrint_ << "\n";
		size_t Nv = vertexData_.size(); 
		if( !Nv ) {
			os << "0 0" << std::endl << "0 0" << std::endl;
			return;
		}
		//Print the vertex positions
		size_t dimension = vertexData_[0].size();
		os << T_->numVertex() << " " << dimension << std::endl;
		for( size_t i=0 ; i<Nv ; ++i ) {
			for( size_t d=0 ; d<dimension ; ++d )
				os << vertexData_[i][d] << " ";
			os << std::endl;
		}
		os << std::endl;
		//Print the cells, first connected vertecis and then variables
		size_t Nc = cellData_.size();
		// For membrane PIN1 printing wall version
		//////////////////////////////////////////////////////////////////////
		os << Nc << std::endl;
		size_t pinI=8,xI=1;
		//size_t auxinI=4,mI=7;
		std::vector<double> parameter(8);
		parameter[0]=0.1;
		parameter[3]=0.01;
		for( size_t i=0 ; i<Nc ; ++i ) {
			size_t Ncv = T_->cell(i).numVertex(); 
			os << Ncv << " ";
			for( size_t k=0 ; k<Ncv ; ++k )
				os << T_->cell(i).vertex(k)->index() << " ";
			
			//pin polarization
			//////////////////////////////////////////////////////////////////////
			size_t numWalls=T_->cell(i).numWall();
			
			//Polarization coefficient normalization constant
			double sum=0.0;
			std::vector<double> Pij(numWalls);
			for( size_t n=0 ; n<numWalls ; ++n ) {
				size_t neighI = T_->cell(i).wall(n)->index();
				sum += Pij[n] = wallData_[ neighI ][ xI ];
				//sum += Pij[n] = (1.0-parameter[2]) + 
				//parameter[2]*cellData_[ neighI ][xI];
			}
			if( sum >= 0.0 )
				os << parameter[3]*cellData_[i][pinI] / sum << " ";
			else
				os << "0.0 ";
			
			for( size_t n=0 ; n<numWalls ; ++n ) {
				double pol=0.0;
				if( sum != 0.0 )
					pol = cellData_[i][pinI] * Pij[n] / sum;
				os << pol << " ";
			}
			os << std::endl;
		}
    //////////////////////////////////////////////////////////////////////
    //End Pij printing wall version
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

	// Increase resolution 
	unsigned int oldPrecision = os.precision(); 
	os.precision(20);
	std::cerr << "Tissue::printInit(): old precision: " << oldPrecision << " new " 
						<< os.precision() << std::endl;	

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
  os << T_->numWall() << " 1 " << wallData_[0].size()-1 << std::endl;
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
	os.precision(oldPrecision);
}

void BaseSolver::printInitFem(std::ostream &os) const
{
	assert( T_->numCell()==cellData_.size() && 
					T_->numWall()==wallData_.size() &&
					T_->numVertex()==vertexData_.size() );
	
	// Increase resolution 
	unsigned int oldPrecision = os.precision(); 
	os.precision(20);
	//std::cerr << "Tissue::printInit(): old precision: " << oldPrecision << " new " 
	//				<< os.precision() << std::endl;	
	
  //Print the vertex positions
	size_t numV=vertexData_.size();
	os << numV << " nodes" << std::endl;
  for (size_t i=0; i<numV; ++i) {
		assert( T_->vertex(i).numPosition()==vertexData_[i].size() );
		os << i << " : ";
    for (size_t j=0; j<T_->vertex(i).numPosition(); ++j )
      os << vertexData_[i][j] << " ";
    os << std::endl;
  }
  //os << std::endl;
  
  //Print cell connection data
	size_t numC=T_->numCell();
  os << numC << " faces" << std::endl;
	for (size_t i=0; i<numC; ++i) {
		os << i << " : ";
		size_t numV=T_->cell(i).numVertex();
		os << numV << ", ";
		for (size_t k=0; k<numV; ++k)
			os << T_->cell(i).vertex(k)->index() << " ";
		os << std::endl;
	}
	os.precision(oldPrecision);
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

