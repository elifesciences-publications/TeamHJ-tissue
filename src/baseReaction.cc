/**
 * Filename     : baseReaction.cc
 * Description  : A base class describing variable updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : October 2003
 * Revision     : $Id: baseReaction.cc,v 1.25 2006/03/18 00:05:14 henrik Exp $
 */
#include<vector>

#include"baseReaction.h"
#include"growth.h"
#include"mechanical.h"
#include"network.h"
//#include"massAction.h"

BaseReaction::~BaseReaction(){}

//!Factory creator, all creation should be mapped onto this one 
/*! Given the idValue a reaction of the defined type is returned
(using new Class).*/
BaseReaction *
BaseReaction::createReaction(std::vector<double> &paraValue,
			       std::vector< std::vector<size_t> > &indValue, 
			       std::string idValue ) {
  
  //Growth related updates
  //growth.h,growth.cc
  if(idValue=="WallGrowthExponentialTruncated")
    return new WallGrowthExponentialTruncated(paraValue,indValue);
  else if(idValue=="WallGrowthExponentialStressTruncated")
    return new WallGrowthExponentialStressTruncated(paraValue,indValue);
  else if(idValue=="WallGrowthConstantStress")
    return new WallGrowthConstantStress(paraValue,indValue);
  else if(idValue=="WallGrowthConstantStressEpidermalAsymmetric")
    return new WallGrowthConstantStressEpidermalAsymmetric(paraValue,indValue);
  else if(idValue=="MoveVertexRadially")
    return new MoveVertexRadially(paraValue,indValue);
  else if (idValue == "WallLengthGrowExperimental")
	  return new WallLengthGrowExperimental(paraValue, indValue);

  //Mechanical interactions between vertices
  //mechanical.h,mechanical.cc
  else if(idValue=="VertexFromWallSpringAsymmetric")
    return new VertexFromWallSpringAsymmetric(paraValue,indValue);
  else if(idValue=="VertexFromWallSpringMT")
    return new VertexFromWallSpringMT(paraValue,indValue);
  else if(idValue=="VertexFromEpidermalWallSpringAsymmetric")
    return new VertexFromEpidermalWallSpringAsymmetric(paraValue,indValue);
  else if(idValue=="VertexFromEpidermalCellWallSpringAsymmetric")
    return new VertexFromEpidermalCellWallSpringAsymmetric(paraValue,indValue);
  else if(idValue=="VertexFromCellPowerdiagram")
    return new VertexFromCellPowerdiagram(paraValue,indValue);
  else if(idValue=="VertexFromCellPressure")
    return new VertexFromCellPressure(paraValue,indValue);
  else if(idValue=="VertexFromCellPressureVolumeNormalized")
    return new VertexFromCellPressureVolumeNormalized(paraValue,indValue);
  else if(idValue=="VertexFromCellPressureThresholdFromMaxPos")
    return new VertexFromCellPressureThresholdFromMaxPos(paraValue,indValue);
  else if(idValue=="VertexFromCellInternalPressure")
    return new VertexFromCellInternalPressure(paraValue,indValue);
  else if(idValue=="VertexNoUpdateFromPosition")
    return new VertexNoUpdateFromPosition(paraValue,indValue); 
  else if(idValue=="VertexForceOrigoFromIndex")
    return new VertexForceOrigoFromIndex(paraValue,indValue); 
  else if(idValue=="CellForceOrigoFromIndex")
    return new CellForceOrigoFromIndex(paraValue,indValue); 
  else if(idValue=="CylinderForce")
    return new CylinderForce(paraValue,indValue); 
  else if(idValue=="SphereCylinderForce")
    return new SphereCylinderForce(paraValue,indValue); 
  else if(idValue=="SphereCylinderForceFromRadius")
    return new SphereCylinderForceFromRadius(paraValue,indValue); 
  else if(idValue=="InfiniteWallForce")
    return new InfiniteWallForce(paraValue,indValue); 
  else if(idValue=="EpidermalVertexForce")
    return new EpidermalVertexForce(paraValue,indValue); 
  else if (idValue == "VertexFromPressureExperimental")
	  return new VertexFromPressureExperimental(paraValue, indValue);
  else if (idValue == "VertexFromWallSpringExperimental")
	  return new VertexFromWallSpringExperimental(paraValue, indValue);
  else if (idValue == "CellVolumeExperimental")
	  return new CellVolumeExperimental(paraValue, indValue);
  else if (idValue == "EpidermalRadialForce")
	  return new EpidermalRadialForce(paraValue, indValue);

	//network.h,network.cc
  else if(idValue=="AuxinModelSimple1")
    return new AuxinModelSimple1(paraValue,indValue); 
  else if(idValue=="AuxinModelSimple2")
    return new AuxinModelSimple2(paraValue,indValue); 
  else if(idValue=="AuxinModelSimple3")
    return new AuxinModelSimple3(paraValue,indValue); 
  else if(idValue=="AuxinModel4")
    return new AuxinModel4(paraValue,indValue); 
  //Default, if nothing found
  else {
    std::cerr << "\nBaseReaction::createReaction() WARNING: Reactiontype " 
							<< idValue << " not known, no reaction created.\n\7";
    exit(-1);
  }
}

//!This creator reads from an open file and then calls for the main creator
BaseReaction* 
BaseReaction::createReaction(std::istream &IN ) {
  
  std::string idVal;
  size_t pNum,levelNum;
  IN >> idVal;
  IN >> pNum;
  IN >> levelNum;
  std::vector<size_t> varIndexNum( levelNum );
  for( size_t i=0 ; i<levelNum ; i++ )
    IN >> varIndexNum[i];
  
  std::vector<double> pVal( pNum );
  for( size_t i=0 ; i<pNum ; i++ )
    IN >> pVal[i];
  
  std::vector< std::vector<size_t> > varIndexVal( levelNum );
  for( size_t i=0 ; i<levelNum ; i++ )
    varIndexVal[i].resize( varIndexNum[i] );
  
  for( size_t i=0 ; i<levelNum ; i++ )
    for( size_t j=0 ; j<varIndexNum[i] ; j++ )
      IN >> varIndexVal[i][j];
  
  return createReaction(pVal,varIndexVal,idVal);
}

void BaseReaction::
derivs(Tissue &T,std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &walldata,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellderivs, 
       std::vector< std::vector<double> > &wallderivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  std::cerr << "BaseReaction::derivs() should not be used. "
	    << "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseReaction::initiate(Tissue &T) {}

void BaseReaction::update(Tissue &T,double h) {}

void BaseReaction::print( std::ofstream &os ) {
  std::cerr << "BaseReaction::print(ofstream) should not be used. "
						<< "Should always be mapped onto one of the real types.\n";
  exit(0);
}
