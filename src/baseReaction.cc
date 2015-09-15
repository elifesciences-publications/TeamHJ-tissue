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
#include "bending.h"
#include "cellTime.h"
#include "centerTriangulation.h"
#include "creation.h"
#include "degradation.h"
#include "directionReaction.h"
#include "grn.h"
#include "growth.h"
#include "mechanical.h"
#include "mechanicalSpring.h"
#include "mechanicalTRBS.h"
#include "network.h"
#include "transport.h"
#include "sisterVertex.h"
#include "membraneCycling.h"

#include"massAction.h"

BaseReaction::~BaseReaction(){}

BaseReaction *
BaseReaction::createReaction(std::vector<double> &paraValue,
			       std::vector< std::vector<size_t> > &indValue, 
			       std::string idValue ) {
  
  //Growth related updates
  //growth.h,growth.cc
  if(idValue == "WallGrowthExponentialTruncated" ) {
    std::cerr << "Reaction WallGrowthExponentialTruncated has been replaced by WallGrowth::Constant." 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  else if(idValue == "WallGrowth::Constant" )
    return new WallGrowth::Constant(paraValue, indValue);
  else if(idValue == "WallGrowthExponentialStressTruncated") {
    std::cerr << "Reaction WallGrowthExponentialStressTruncated "
	      << "has been replaced by WallGrowth::Stress (setting the stretch_flag to 1 and provide L_th)." 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  else if(idValue == "WallGrowthStress" || idValue == "WallGrowth::Stress" )
    return new WallGrowth::Stress(paraValue, indValue);
  else if(idValue == "WallGrowthStrain" || idValue == "WallGrowth::Strain" )
    return new WallGrowth::Strain(paraValue, indValue);
  else if (idValue == "WallGrowth::CenterTriangulation::Constant" ||
	   idValue == "CenterTriangulation::WallGrowth::Constant")
    return new WallGrowth::CenterTriangulation::Constant(paraValue, indValue);
  else if (idValue == "WallGrowthStresscenterTriangulation" ||
	   idValue == "WallGrowth::CenterTriangulation::Stress" ||
	   idValue == "CenterTriangulation::WallGrowth::Stress")
    return new WallGrowth::CenterTriangulation::Stress(paraValue, indValue);
  else if (idValue == "CenterTriangulation::WallGrowth::StrainTRBS")
    return new WallGrowth::CenterTriangulation::StrainTRBS(paraValue, indValue);    
  else if(idValue == "WallGrowthStressSpatial" || idValue == "WallGrowth::StressSpatial")
    return new WallGrowth::StressSpatial(paraValue, indValue);
  else if(idValue == "WallGrowthStressSpatialSingle" || idValue == "WallGrowth::StressSpatialSingle")
    return new WallGrowth::StressSpatialSingle(paraValue, indValue);
  else if(idValue == "WallGrowthStressConcentrationHill" || idValue == "WallGrowth::StressConcentrationHill")
    return new WallGrowth::StressConcentrationHill(paraValue, indValue);
  else if(idValue == "WallGrowthConstantStressEpidermalAsymmetric" || 
	  idValue == "WallGrowth::ConstantStressEpidermalAsymmetric")
    return new WallGrowth::ConstantStressEpidermalAsymmetric(paraValue, indValue);
  else if (idValue == "WallLengthGrowExperimental") {
    std::cerr << "Reaction WallLengthGrowExperimental "
	      << "has been replaced by WallGrowth::Force. Better is to use the WallGrowth::Stress "
	      << "(setting the stretch_flag to 0 and not provide L_th)." 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  else if (idValue == "WallGrowthConstantStress" || 
	   idValue == "WallGrowthConstantStressConcentrationHill") {
    std::cerr << "BaseReaction::createReaction() WallGrowthConstantStress* has been "
	      << "replaced by WallGrowth::Stress*." << std::endl;
    exit(EXIT_FAILURE);
  }
  else if (idValue == "WallGrowth::Force")
    return new WallGrowth::Force(paraValue, indValue);
  else if(idValue == "MoveVertexRadially")
    return new MoveVertexRadially(paraValue, indValue);
  else if(idValue == "MoveEpidermalVertexRadially")
    return new MoveEpidermalVertexRadially(paraValue, indValue);
  else if(idValue == "MoveVerteX")
    return new MoveVerteX(paraValue, indValue);
  else if(idValue == "MoveVertexY")
      return new MoveVertexY(paraValue, indValue);
  else if(idValue == "MoveVertexRadiallycenterTriangulation")
    return new MoveVertexRadiallycenterTriangulation(paraValue, indValue);
  else if(idValue == "MoveVertexSphereCylinder")
    return new MoveVertexSphereCylinder(paraValue, indValue);
  else if (idValue == "WaterVolumeFromTurgor")
    return new WaterVolumeFromTurgor(paraValue, indValue);
  else if (idValue == "DilutionFromVertexDerivs")
    return new DilutionFromVertexDerivs(paraValue, indValue);
  
  //Mechanical interactions between vertices
  //mechanicalSpring.h,mechanicalSpring.cc
  else if(idValue=="VertexFromWallSpring")
    return new VertexFromWallSpring(paraValue,indValue);
  else if(idValue=="VertexFromWallSpringMTnew")
    return new VertexFromWallSpringMTnew(paraValue,indValue);
  else if(idValue=="VertexFromWallBoundarySpring")
    return new VertexFromWallBoundarySpring(paraValue,indValue);
  else if(idValue=="CenterTriangulation::EdgeSpring")
    return new CenterTriangulation::EdgeSpring(paraValue,indValue);
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
  else if(idValue=="VertexFromExternalSpring")
    return new VertexFromExternalSpring(paraValue,indValue);
  else if(idValue=="VertexFromExternalSpringFromPerpVertex")
    return new VertexFromExternalSpringFromPerpVertex(paraValue,indValue);
 else if(idValue=="VertexFromExternalSpringFromPerpVertexDynamic")
    return new VertexFromExternalSpringFromPerpVertexDynamic(paraValue,indValue);
 else if(idValue=="cellcellRepulsion")
    return new cellcellRepulsion(paraValue,indValue);
 else if(idValue=="vertexFromSubstrate")
   return new vertexFromSubstrate(paraValue,indValue);
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
 else if(idValue=="CenterTriangulation::VertexFromCellPressure" ||
         idValue=="VertexFromCellPressurecenterTriangulation")
   return new CenterTriangulation::VertexFromCellPressure(paraValue,indValue);
 else if(idValue=="CenterTriangulation::VertexFromCellPressureLinear" ||
         idValue=="VertexFromCellPressurecenterTriangulationLinear")
   return new CenterTriangulation::VertexFromCellPressureLinear(paraValue,indValue);
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
  else if (idValue == "VertexFromCellPlaneLinear")
    return new VertexFromCellPlaneLinear(paraValue, indValue);
  else if (idValue == "VertexFromCellPlaneLinearCenterTriangulation")
    return new VertexFromCellPlaneLinearCenterTriangulation(paraValue, indValue);
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
  else if(idValue=="VertexFromForceLinear")
    return new VertexFromForceLinear(paraValue,indValue);
  else if(idValue=="VertexFromBall")
    return new VertexFromBall(paraValue,indValue);
  else if(idValue=="VertexFromParabolid")
    return new VertexFromParabolid(paraValue,indValue);
  else if(idValue=="VertexFromExternalWall")
    return new VertexFromExternalWall(paraValue,indValue);
  else if(idValue=="TemplateVolumeChange")
    return new TemplateVolumeChange(paraValue,indValue);
  else if(idValue=="CalculateAngleVectors")
    return new CalculateAngleVectors(paraValue,indValue);
  else if(idValue=="CalculateAngleVectorXYplane")
    return new CalculateAngleVectorXYplane(paraValue,indValue);
  else if(idValue=="AngleVector")
    return new AngleVector(paraValue,indValue);
  else if(idValue=="VertexFromHypocotylGrowth")
    return new VertexFromHypocotylGrowth(paraValue,indValue);
  else if(idValue=="maxVelocity")
    return new maxVelocity(paraValue,indValue);
  else if (idValue == "DebugReaction")
    return new DebugReaction(paraValue, indValue);
  
  // centerTriangulation.h (.cc)
  // Reactions related to a center triangulation of cells
  else if (idValue == "CenterTriangulation::Initiate")
    return new CenterTriangulation::Initiate(paraValue, indValue);
  
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
  else if (idValue == "FiberModel")
    return new FiberModel(paraValue, indValue);
  else if (idValue == "FiberDeposition")
    return new FiberDeposition(paraValue, indValue);

  // bending.h (.cc)
  else if (idValue == "Bending::NeighborCenter")
    return new Bending::NeighborCenter(paraValue, indValue);
  else if (idValue == "Bending::Angle")
    return new Bending::Angle(paraValue, indValue);
  else if (idValue == "Bending::AngleInitiate")
    return new Bending::AngleInitiate(paraValue, indValue);
  else if (idValue == "Bending::AngleRelax")
    return new Bending::AngleRelax(paraValue, indValue);

  //creation.h,creation.cc
  else if(idValue=="CreationZero")
    return new CreationZero(paraValue,indValue); 
  else if(idValue=="CreationOne")
    return new CreationOne(paraValue,indValue); 
  else if(idValue=="CreationTwo")
    return new CreationTwo(paraValue,indValue); 
  else if(idValue=="CreationSpatialSphere")
    return new CreationSpatialSphere(paraValue,indValue); 
  else if(idValue=="CreationSpatialRing")
    return new CreationSpatialRing(paraValue,indValue); 
  else if(idValue=="CreationSpatialCoordinate")
    return new CreationSpatialCoordinate(paraValue,indValue); 
  else if(idValue=="CreationFromList")
    return new CreationFromList(paraValue,indValue); 
  else if(idValue=="CreationOneGeometric")
    return new CreationOneGeometric(paraValue,indValue); 
  
  
  //degradation.h,degradation.cc
  else if(idValue=="DegradationOne")
    return new DegradationOne(paraValue,indValue); 
  else if(idValue=="DegradationTwo")
    return new DegradationTwo(paraValue,indValue); 
  else if(idValue=="DegradationTwoGeometric")
    return new DegradationTwoGeometric(paraValue,indValue); 
  else if(idValue=="DegradationHill")
    return new DegradationHill(paraValue,indValue); 
  //grn.h,grn.cc
  else if(idValue=="Hill")
    return new Hill(paraValue,indValue); 
  else if(idValue=="HillGeneralOne")
    return new HillGeneralOne(paraValue,indValue); 
  else if(idValue=="HillGeneralTwo")
    return new HillGeneralTwo(paraValue,indValue); 
  else if(idValue=="HillGeneralThree")
    return new HillGeneralThree(paraValue,indValue); 
  else if(idValue=="Grn")
    return new Grn(paraValue,indValue); 
  else if(idValue=="Gsrn2")
    return new Gsrn2(paraValue,indValue); 



  //transport.h,transport.cc
  else if(idValue=="MembraneDiffusionSimple")
    return new MembraneDiffusionSimple(paraValue,indValue); 
  else if(idValue=="MembraneDiffusionSimple2")
    return new MembraneDiffusionSimple2(paraValue,indValue); 
  else if(idValue=="DiffusionSimple")
    return new DiffusionSimple(paraValue,indValue);
  else if(idValue=="ActiveTransportCellEfflux")
    return new ActiveTransportCellEfflux(paraValue,indValue);
  else if(idValue=="ActiveTransportCellEffluxMM")
    return new ActiveTransportCellEffluxMM(paraValue,indValue);
  else if(idValue=="ActiveTransportWall")
    return new ActiveTransportWall(paraValue,indValue);


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
  else if(idValue=="AuxinModelSimple4")
    return new AuxinModelSimple4(paraValue,indValue); 
  else if(idValue=="AuxinModelSimple5")
    return new AuxinModelSimple5(paraValue,indValue); 
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
  else if(idValue=="AuxinPINBistabilityModel")
    return new AuxinPINBistabilityModel(paraValue,indValue); 
  else if(idValue=="AuxinPINBistabilityModelCell")
    return new AuxinPINBistabilityModelCell(paraValue,indValue); 
  else if(idValue=="AuxinExoBistability")
    return new AuxinExoBistability(paraValue,indValue); 
  else if(idValue=="AuxinPINBistabilityModelCellNew")
    return new AuxinExoBistability(paraValue,indValue); 
  else if(idValue=="SimpleROPModel")
    return new SimpleROPModel(paraValue,indValue); 
  else if(idValue=="SimpleROPModel2")
    return new SimpleROPModel2(paraValue,indValue); 
  else if(idValue=="SimpleROPModel3")
    return new SimpleROPModel3(paraValue,indValue); 
  else if(idValue=="SimpleROPModel4")
    return new SimpleROPModel4(paraValue,indValue);  
  else if(idValue=="SimpleROPModel5")
    return new SimpleROPModel5(paraValue,indValue);  
  else if(idValue=="SimpleROPModel6")
    return new SimpleROPModel6(paraValue,indValue);  
  else if(idValue=="SimpleROPModel7")
    return new SimpleROPModel7(paraValue,indValue); 
  else if(idValue=="UpInternalGradientModel")
    return new UpInternalGradientModel(paraValue,indValue);
    
   else if(idValue=="DownInternalGradientModel")
    return new DownInternalGradientModel(paraValue,indValue);
   else if(idValue=="DownInternalGradientModelGeometric")
    return new DownInternalGradientModelGeometric(paraValue,indValue);
    
  else if(idValue=="DownInternalGradientModelSingleCell")
    return new DownInternalGradientModelSingleCell(paraValue,indValue);
    
  else if(idValue=="UpExternalGradientModel")
    return new UpExternalGradientModel(paraValue,indValue);
    
  else if(idValue=="UpInternalGradientModel")
    return new UpInternalGradientModel(paraValue,indValue);
    
  else if(idValue=="AuxinFluxModel")
    return new AuxinFluxModel(paraValue,indValue);
    
 else if(idValue=="IntracellularPartitioning")
    return new IntracellularPartitioning(paraValue,indValue);
    
 else if(idValue=="IntracellularCoupling")
    return new IntracellularCoupling(paraValue,indValue);
    
 else if(idValue=="IntracellularIndirectCoupling")
    return new IntracellularIndirectCoupling(paraValue,indValue); 

  //directionReaction.h, directionUpdate.cc
  else if (idValue == "ContinousMTDirection")
    return new ContinousMTDirection(paraValue, indValue);
  else if (idValue == "ContinousMTDirection3d")
  return new ContinousMTDirection3d(paraValue, indValue);
  else if (idValue == "UpdateMTDirection")
    return new UpdateMTDirection(paraValue, indValue);
  else if (idValue == "UpdateMTDirectionEquilibrium")
    return new UpdateMTDirectionEquilibrium(paraValue, indValue);
  else if (idValue == "UpdateMTDirectionConcenHill")
    return new UpdateMTDirectionConcenHill(paraValue, indValue);
  else if (idValue == "RotatingDirection")
    return new RotatingDirection(paraValue, indValue);
  
  //sisterVertex.h, sisterVertex.cc
  else if (idValue == "SisterVertex::InitiateFromFile")
    return new SisterVertex::InitiateFromFile(paraValue, indValue);
  else if (idValue == "SisterVertex::InitiateFromDistance")
    return new SisterVertex::InitiateFromDistance(paraValue, indValue);
  else if (idValue == "SisterVertex::Spring")
    return new SisterVertex::Spring(paraValue, indValue);
  else if (idValue == "SisterVertex::CombineDerivatives")
    return new SisterVertex::CombineDerivatives(paraValue, indValue);
  
  //adhocReaction.h,adhocReaction.cc
  else if(idValue=="VertexNoUpdateFromPosition")
    return new VertexNoUpdateFromPosition(paraValue,indValue); 
  else if(idValue=="VertexNoUpdateFromIndex")
    return new VertexNoUpdateFromIndex(paraValue,indValue);
  else if(idValue=="VertexNoUpdateFromList")
    return new VertexNoUpdateFromList(paraValue,indValue);
  else if(idValue=="VertexRandTip")
    return new VertexRandTip(paraValue,indValue); 
  else if(idValue=="VertexNoUpdateBoundary")
    return new VertexNoUpdateBoundary(paraValue,indValue); 
  else if(idValue=="VertexNoUpdateBoundaryPtemplate")
    return new VertexNoUpdateBoundaryPtemplate(paraValue,indValue); 
  else if(idValue=="VertexNoUpdateBoundaryPtemplateStatic")
    return new VertexNoUpdateBoundaryPtemplateStatic(paraValue,indValue); 
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
  else if(idValue=="MoveVerticesRandomlyCapCylinder")
    return new MoveVerticesRandomlyCapCylinder(paraValue,indValue); 
  else if(idValue=="scaleTemplate")
    return new scaleTemplate(paraValue,indValue); 
  else if(idValue=="copyCellVector")
    return new copyCellVector(paraValue,indValue); 
  else if(idValue=="restrictVertexRadially")
    return new restrictVertexRadially(paraValue,indValue); 
  else if(idValue=="VertexFromRotationalForceLinear")
    return new VertexFromRotationalForceLinear(paraValue,indValue); 
  

  // cellTime.h
  else if (idValue=="CellTimeDerivative")
    return new CellTimeDerivative(paraValue, indValue);

 // MembraneCycling.h
  else if (idValue=="MembraneCycling::Constant")
    return new MembraneCycling::Constant(paraValue, indValue);
  else if (idValue=="MembraneCycling::CrossMembraneNonLinear")
    return new MembraneCycling::CrossMembraneNonLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::LocalWallFeedbackNonLinear")
    return new MembraneCycling::LocalWallFeedbackNonLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::LocalWallFeedbackLinear")
    return new MembraneCycling::LocalWallFeedbackLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::CellUpTheGradientNonLinear")
    return new MembraneCycling::CellUpTheGradientNonLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::CellUpTheGradientLinear")
    return new MembraneCycling::CellUpTheGradientLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::PINFeedbackNonLinear")
    return new MembraneCycling::PINFeedbackNonLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::PINFeedbackLinear")
    return new MembraneCycling::PINFeedbackLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::InternalCellNonLinear")
    return new MembraneCycling::InternalCellNonLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::InternalCellLinear")
    return new MembraneCycling::InternalCellLinear(paraValue, indValue);
  else if (idValue=="MembraneCycling::CellFluxExocytosis")
    return new MembraneCycling::CellFluxExocytosis(paraValue, indValue);

  //massAction.h
  else if (idValue=="MassAction::General")
    return new MassAction::General(paraValue, indValue);
  else if (idValue=="MassAction::OneToTwo")
    return new MassAction::OneToTwo(paraValue, indValue);
  else if (idValue=="MassAction::TwoToOne")
    return new MassAction::TwoToOne(paraValue, indValue);
  else if (idValue=="MassAction::GeneralWall")
      return new MassAction::GeneralWall(paraValue, indValue);
  else if (idValue=="MassAction::OneToTwoWall")
      return new MassAction::OneToTwoWall(paraValue, indValue);
  else if (idValue=="MassAction::TwoToOneWall")
      return new MassAction::TwoToOneWall(paraValue, indValue);


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
	    << "Should always be mapped onto one of the real types." << std::endl;
  exit(0);
}  

void BaseReaction::
derivsWithAbs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &walldata,
       DataMatrix &vertexData,
       DataMatrix &cellderivs, 
       DataMatrix &wallderivs,
       DataMatrix &vertexDerivs,
       DataMatrix &sdydtCell, 
       DataMatrix &sdydtWall,
       DataMatrix &sdydtVertex ) 
{
  std::cerr << "BaseReaction::derivsWithAbs() should not be used. "
	    << "Should always be mapped onto one of the real types." << std::endl;
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

void BaseReaction::print( std::ofstream &os ) 
{
  std::cerr << "BaseReaction::print(ofstream) should not be used. "
	    << "Should always be mapped onto one of the real types.\n";
  exit(0);
}

void BaseReaction::printState(Tissue *T,
			      DataMatrix &cellData,
			      DataMatrix &wallData,
			      DataMatrix &vertexData, 
			      std::ostream &os)
{
}
