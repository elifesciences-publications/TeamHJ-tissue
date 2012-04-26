//
// Filename     : baseReaction.cc
// Description  : A base class describing variable updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2003
// Revision     : $Id: baseReaction.cc,v 1.25 2006/03/18 00:05:14 henrik Exp $
//
#include <vector>

#include "baseReaction.h"
#include "adhocReaction.h"
#include "cellTime.h"
#include "creation.h"
#include "degradation.h"
#include "directionReaction.h"
#include "growth.h"
#include "mechanical.h"
#include "mechanicalSpring.h"
#include "mechanicalTRBS.h"
#include "network.h"

//#include"massAction.h"

BaseReaction::~BaseReaction(){}

BaseReaction *
BaseReaction::createReaction(std::vector<double> &paraValue,
			       std::vector< std::vector<size_t> > &indValue, 
			       std::string idValue ) {
  
  //Growth related updates
  //growth.h,growth.cc
  if(idValue == "WallGrowthExponentialTruncated")
    return new WallGrowthExponentialTruncated(paraValue, indValue);
  else if(idValue == "WallGrowthExponentialStressTruncated")
    return new WallGrowthExponentialStressTruncated(paraValue, indValue);
  else if(idValue == "WallGrowthStress")
    return new WallGrowthStress(paraValue, indValue);
  else if(idValue == "WallGrowthStresscenterTriangulation")
    return new WallGrowthStresscenterTriangulation(paraValue, indValue);
  else if(idValue == "WallGrowthStressSpatial")
    return new WallGrowthStressSpatial(paraValue, indValue);
  else if(idValue == "WallGrowthStressSpatialSingle")
    return new WallGrowthStressSpatialSingle(paraValue, indValue);
  else if(idValue == "WallGrowthStressConcentrationHill")
    return new WallGrowthStressConcentrationHill(paraValue, indValue);
  else if(idValue == "WallGrowthConstantStressEpidermalAsymmetric")
    return new WallGrowthConstantStressEpidermalAsymmetric(paraValue, indValue);
  else if(idValue == "MoveVertexRadially")
    return new MoveVertexRadially(paraValue, indValue);
  else if(idValue == "MoveVertexRadiallycenterTriangulation")
    return new MoveVertexRadiallycenterTriangulation(paraValue, indValue);
  else if(idValue == "MoveVertexSphereCylinder")
    return new MoveVertexSphereCylinder(paraValue, indValue);
  else if (idValue == "WallLengthGrowExperimental")
    return new WallLengthGrowExperimental(paraValue, indValue);
  else if (idValue == "WaterVolumeFromTurgor")
    return new WaterVolumeFromTurgor(paraValue, indValue);
  else if (idValue == "DilutionFromVertexDerivs")
    return new DilutionFromVertexDerivs(paraValue, indValue);
  else if (idValue == "WallGrowthConstantStress" || 
	   idValue == "WallGrowthConstantStressConcentrationHill") {
    std::cerr << "BaseReaction::createReaction() WallGrowthConstantStress* has been "
	      << "replaced by WallGrowthStress*." << std::endl;
    exit(-1);
  }
  
  //Mechanical interactions between vertices
  //mechanicalSpring.h,mechanicalSpring.cc
  else if(idValue=="VertexFromWallSpring")
    return new VertexFromWallSpring(paraValue,indValue);
  else if(idValue=="VertexFromDoubleWallSpring")
    return new VertexFromDoubleWallSpring(paraValue,indValue);
  else if(idValue=="VertexFromWallSpringSpatial")
    return new VertexFromWallSpringSpatial(paraValue,indValue);
  else if(idValue=="VertexFromWallSpringConcentrationHill")
    return new VertexFromWallSpringConcentrationHill(paraValue,indValue);
  else if(idValue=="VertexFromWallSpringMT")
    return new VertexFromWallSpringMT(paraValue,indValue);
  else if(idValue=="VertexFromWallSpringMTSpatial")
    return new VertexFromWallSpringMTSpatial(paraValue,indValue);
  else if(idValue=="VertexFromWallSpringMTHistory")
    return new VertexFromWallSpringMTHistory(paraValue,indValue);
  else if(idValue=="VertexFromEpidermalWallSpring")
    return new VertexFromEpidermalWallSpring(paraValue,indValue);
  else if(idValue=="VertexFromEpidermalCellWallSpring")
    return new VertexFromEpidermalCellWallSpring(paraValue,indValue);
  else if (idValue == "VertexFromWallSpringExperimental")
    return new VertexFromWallSpringExperimental(paraValue, indValue);
  else if(idValue=="VertexFromWallSpringMTConcentrationHill")
    return new VertexFromWallSpringMTConcentrationHill(paraValue,indValue);
  else if(idValue=="VertexFromDoubleWallSpringMTConcentrationHill")
    return new VertexFromDoubleWallSpringMTConcentrationHill(paraValue,indValue);
  else if (idValue=="VertexFromWallSpringAsymmetric" ||
	   idValue=="VertexFromEpidermalWallSpringAsymmetric" ||
	   idValue=="VertexFromEpidermalCellWallSpringAsymmetric") {
    std::cerr << "BaseReaction::BaseReaction() All *SpringAsymmetric have been renamed "
	      << "*Spring." << std::endl;
    exit(-1);
  }
  //Mechanical interactions between vertices
  //mechanical.h,mechanical.cc
  else if(idValue=="VertexFromCellPowerdiagram")
    return new VertexFromCellPowerdiagram(paraValue,indValue);
  else if(idValue=="VertexFromCellPressure")
    return new VertexFromCellPressure(paraValue,indValue);
  else if(idValue=="VertexFromCellPressurecenterTriangulation")
    return new VertexFromCellPressurecenterTriangulation(paraValue,indValue);
  else if(idValue=="VertexFromCellPressureVolumeNormalized")
    return new VertexFromCellPressureVolumeNormalized(paraValue,indValue);
  else if(idValue=="VertexFromCellPressureThresholdFromMaxPos")
    return new VertexFromCellPressureThresholdFromMaxPos(paraValue,indValue);
  else if(idValue=="VertexFromCellInternalPressure")
    return new VertexFromCellInternalPressure(paraValue,indValue);
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
  else if (idValue == "CellVolumeExperimental")
    return new CellVolumeExperimental(paraValue, indValue);
  else if (idValue == "EpidermalRadialForce")
    return new EpidermalRadialForce(paraValue, indValue);
  else if (idValue == "PerpendicularWallPressure")
    return new PerpendicularWallPressure(paraValue, indValue);
  else if (idValue == "VertexFromCellPlane")
    return new VertexFromCellPlane(paraValue, indValue);
  else if (idValue == "VertexFromCellPlaneSpatial")
    return new VertexFromCellPlaneSpatial(paraValue, indValue);
  else if (idValue == "VertexFromCellPlaneConcentrationHill")
    return new VertexFromCellPlaneConcentrationHill(paraValue, indValue);
  else if (idValue == "VertexFromCellPlaneNormalized")
    return new VertexFromCellPlaneNormalized(paraValue, indValue);
  else if (idValue == "VertexFromCellPlaneNormalizedSpatial")
    return new VertexFromCellPlaneNormalizedSpatial(paraValue, indValue);
  else if (idValue == "VertexFromCellPlaneSphereCylinder")
    return new VertexFromCellPlaneSphereCylinder(paraValue, indValue);
  else if (idValue == "VertexFromCellPlaneSphereCylinderConcentrationHill")
    return new VertexFromCellPlaneSphereCylinderConcentrationHill(paraValue, indValue);
  else if (idValue == "VertexFromCellPlaneTriangular")
    return new VertexFromCellPlaneTriangular(paraValue, indValue);
  else if(idValue=="VertexFromForce")
    return new VertexFromForce(paraValue,indValue);
  else if(idValue=="VertexFromBall")
    return new VertexFromBall(paraValue,indValue);
  else if (idValue == "DebugReaction")
    return new DebugReaction(paraValue, indValue);
  
  // mechanicalTRBS.h (.cc)
  // Mechanical updates related to triangular (biquadratic) springs
  else if (idValue == "VertexFromTRBS")
    return new VertexFromTRBS(paraValue, indValue);
  else if (idValue == "VertexFromTRBScenterTriangulation")
    return new VertexFromTRBScenterTriangulation(paraValue, indValue);
  else if (idValue == "VertexFromTRBScenterTriangulationConcentrationHill")
    return new VertexFromTRBScenterTriangulationConcentrationHill(paraValue, indValue);
  else if (idValue == "VertexFromTRBSMT")
    return new VertexFromTRBSMT(paraValue, indValue);
  else if (idValue == "VertexFromTRBScenterTriangulationMT")
    return new VertexFromTRBScenterTriangulationMT(paraValue, indValue);
 else if (idValue == "VertexFromTRBScenterTriangulationConcentrationHillMT")
   return new VertexFromTRBScenterTriangulationConcentrationHillMT(paraValue, indValue);
  

  //creation.h,creation.cc
  else if(idValue=="CreationZero")
    return new CreationZero(paraValue,indValue); 
  else if(idValue=="CreationOne")
    return new CreationOne(paraValue,indValue); 
  //degradation.h,degradation.cc
  else if(idValue=="DegradationOne")
    return new DegradationOne(paraValue,indValue); 
  else if(idValue=="DegradationTwo")
    return new DegradationTwo(paraValue,indValue); 

  //network.h,network.cc
  else if(idValue=="AuxinModelSimple1")
    return new AuxinModelSimple1(paraValue,indValue); 
  else if(idValue=="AuxinModelStress")
    return new AuxinModelStress(paraValue,indValue); 
  else if(idValue=="AuxinModelSimpleStress")
    return new AuxinModelSimpleStress(paraValue,indValue); 
  else if(idValue=="AuxinModelSimple1Wall")
    return new AuxinModelSimple1Wall(paraValue,indValue); 
  else if(idValue=="AuxinModelSimple2")
    return new AuxinModelSimple2(paraValue,indValue); 
  else if(idValue=="AuxinModelSimple3")
    return new AuxinModelSimple3(paraValue,indValue); 
  else if(idValue=="AuxinModel4")
    return new AuxinModel4(paraValue,indValue); 
  else if(idValue=="AuxinModel5")
    return new AuxinModel5(paraValue,indValue); 
  else if(idValue=="AuxinModel6")
    return new AuxinModel6(paraValue,indValue); 
  else if(idValue=="AuxinModel7")
    return new AuxinModel7(paraValue,indValue); 
  else if(idValue=="AuxinTransportCellCellNoGeometry")
    return new AuxinTransportCellCellNoGeometry(paraValue,indValue); 
  else if(idValue=="AuxinWallModel")
    return new AuxinWallModel(paraValue,indValue); 
  else if(idValue=="AuxinROPModel")
    return new AuxinROPModel(paraValue,indValue); 
  else if(idValue=="AuxinROPModel2")
    return new AuxinROPModel2(paraValue,indValue); 
  else if(idValue=="AuxinROPModel3")
    return new AuxinROPModel3(paraValue,indValue); 

  //directionReaction.h, directionUpdate.cc
  else if (idValue == "ContinousMTDirection")
    return new ContinousMTDirection(paraValue, indValue);
  else if (idValue == "UpdateMTDirection")
    return new UpdateMTDirection(paraValue, indValue);
  else if (idValue == "RotatingDirection")
    return new RotatingDirection(paraValue, indValue);
  
  //adhocReaction.h,adhocReaction.cc
  else if(idValue=="VertexNoUpdateFromPosition")
    return new VertexNoUpdateFromPosition(paraValue,indValue); 
  else if(idValue=="VertexNoUpdateBoundary")
    return new VertexNoUpdateBoundary(paraValue,indValue); 
  else if(idValue=="VertexTranslateToMax")
    return new VertexTranslateToMax(paraValue,indValue); 
  else if(idValue=="CenterCOM")
    return new CenterCOM(paraValue,indValue); 
  else if(idValue=="CenterCOMcenterTriangulation")
    return new CenterCOMcenterTriangulation(paraValue,indValue); 
  else if(idValue=="CalculatePCAPlane")
    return new CalculatePCAPlane(paraValue,indValue); 
  else if(idValue=="InitiateWallLength")
    return new InitiateWallLength(paraValue,indValue); 
  else if(idValue=="InitiateWallMesh")
    return new InitiateWallMesh(paraValue,indValue); 
  else if(idValue=="StrainTest")
    return new StrainTest(paraValue,indValue); 
  else if(idValue=="CalculateVertexStressDirection")
    return new CalculateVertexStressDirection(paraValue,indValue); 

  // cellTime.h
  else if (idValue=="CellTimeDerivative")
    return new CellTimeDerivative(paraValue, indValue);
	
  // Default, if nothing found
  else {
    std::cerr << "\nBaseReaction::createReaction() WARNING: Reactiontype " 
	      << idValue << " not known, no reaction created.\n\7";
    exit(-1);
  }
}

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
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &walldata,
       DataMatrix &vertexData,
       DataMatrix &cellderivs, 
       DataMatrix &wallderivs,
       DataMatrix &vertexDerivs ) 
{
  std::cerr << "BaseReaction::derivs() should not be used. "
	    << "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseReaction::initiate(Tissue &T,
			    DataMatrix &cellData,
			    DataMatrix &walldata,
			    DataMatrix &vertexData,
			    DataMatrix &cellderivs, 
			    DataMatrix &wallderivs,
			    DataMatrix &vertexDerivs )
{
}

void BaseReaction::update(Tissue &T,
			  DataMatrix &cellData,
			  DataMatrix &walldata,
			  DataMatrix &vertexData,
			  double h) 
{
}

void BaseReaction::print( std::ofstream &os ) {
  std::cerr << "BaseReaction::print(ofstream) should not be used. "
	    << "Should always be mapped onto one of the real types.\n";
  exit(0);
}
