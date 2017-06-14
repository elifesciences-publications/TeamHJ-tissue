/**
 * Filename     : adhocReaction.cc
 * Description  : Classes describing some ad hoc updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : September 2007
 * Revision     : $Id:$
 */

#include"adhocReaction.h"
#include"tissue.h"
#include"baseReaction.h"
#include<cmath>
#include<cstdlib>
#include<ctime>
#include <fstream>
#include <iostream>
#include <sstream>

VertexNoUpdateFromPosition::
VertexNoUpdateFromPosition(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "Uses two parameters threshold and direction "
	      << "(-1 -> less than)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "direction (second parameter) need to be 1 (greater than) "
	      << "or -1 (less than)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexNoUpdateFromPosition");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "threshold";
  tmp[1] = "direction";
  setParameterId( tmp );
}

void VertexNoUpdateFromPosition::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  size_t dimension = vertexData[0].size();
  assert( posIndex < vertexData[0].size() );
  
  for( size_t i=0 ; i<numVertices ; ++i )
    if( (parameter(1)>0 && vertexData[i][posIndex]>parameter(0)) || 
	(parameter(1)<0 && vertexData[i][posIndex]<parameter(0)) )
      for( size_t d=0 ; d<dimension ; ++d )
	vertexDerivs[i][d] = 0.0;
}

VertexNoUpdateFromIndex::
VertexNoUpdateFromIndex(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "VertexNoUpdateFromIndex::"
	      << "VertexNoUpdateFromIndex() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() < 1 ) {
    std::cerr << "VertexNoUpdateFromIndex::"
	      << "VertexNoUpdateFromIndex() "
	      << "Vertex indices in first level." << std::endl;
    exit(0);
  }
  // Set the variable values
  //
  setId("VertexNoUpdateFromIndex");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
}

void VertexNoUpdateFromIndex::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Check the cancelation for vertices with given indices
  size_t dimension = vertexData[0].size();
  for( size_t i=0 ; i<numVariableIndex(0) ; ++i ) {
    size_t k = variableIndex(0,i);
    assert( k < vertexData.size() );
    for( size_t d=0 ; d<dimension ; ++d )
      vertexDerivs[k][d] = 0.0;
  }
}


VertexNoUpdateFromList::
VertexNoUpdateFromList(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "VertexNoUpdateFromList::"
	      << "VertexNoUpdateFromList() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexNoUpdateFromList::"
	      << "VertexNoUpdateFromList() "
	      << "Uses no variable index." << std::endl;
    exit(0);
  }
  // Set the variable values
  //
  setId("VertexNoUpdateFromList");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
}

void VertexNoUpdateFromList::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Check the cancelation for vertices with given indices
  size_t dimension = vertexData[0].size();
  size_t numVertices = T.numVertex();
  
  for (size_t i=0; i<numVertices; ++i)
    if(std::find(updateVertices.begin(), updateVertices.end(), i) == updateVertices.end())    
      for (size_t d=0; d<dimension; ++d)
        vertexDerivs[i][d] = 0.0;   
}

void VertexNoUpdateFromList:: 
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h){
  updateVertices.clear();
  size_t numVertices = T.numVertex();
  for (size_t i=0; i<numVertices; ++i){
    size_t numWalls=T.vertex(i).numWall();
    if(numWalls==2)
      if((T.vertex(i).wall(0)->cell1()==T.background() ||
          T.vertex(i).wall(0)->cell2()==T.background()) &&
         (T.vertex(i).wall(1)->cell1()==T.background() ||
          T.vertex(i).wall(1)->cell2()==T.background()))
        updateVertices.push_back(i);
  }
    // for(size_t i=0; i<updateVertices.size(); ++i)
    //   std::cerr<< updateVertices[i]<<"  ";
    // std::cerr<<std::endl;  
}




VertexRandTip::
VertexRandTip(std::vector<double> &paraValue, 
              std::vector< std::vector<size_t> > 
              &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexRandTip::"
	      << "VertexRandTip() "
	      << "Uses one parameter for max_randomization angle."
	      << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 2 ) {
    std::cerr << "VertexRandTip::"
	      << "VertexRandTip() "
	      << "Uses two variable indices at the first level for " 
              << "tip_signal_index and mt_vector_starting_index." << std::endl;
    exit(0);
  }
  // Set the variable values
  //
  setId("VertexRandTip");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
}

void VertexRandTip::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  

}

void VertexRandTip:: 
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h){
  size_t rotIndex=variableIndex(0,0);
  size_t mtIndex=variableIndex(0,1);
  size_t numCells = T.numCell();



  for (size_t cellInd=0; cellInd<numCells; ++cellInd){
    size_t numWalls=T.cell(cellInd).numWall();
    size_t baseWall, wallsBackground=0;
    for(size_t wallInd=0; wallInd<numWalls; ++wallInd){
      if(T.cell(cellInd).wall(wallInd)->cell1()==T.background() ||
         T.cell(cellInd).wall(wallInd)->cell2()==T.background())
        wallsBackground++;
      if(T.cell(cellInd).wall(wallInd)->cell1()!=T.background() &&
         T.cell(cellInd).wall(wallInd)->cell2()!=T.background())
        baseWall=wallInd;
    }
    if(cellData[cellInd][rotIndex]!=1 && wallsBackground==3){   
      cellData[cellInd][rotIndex]=1;
      // finding the wall M and its middle point
      std::vector<double> midWall(2);
      size_t v1=T.cell(cellInd).wall(baseWall)->vertex1()->index(),
        v2=T.cell(cellInd).wall(baseWall)->vertex2()->index();
      midWall[0]=0.5*(vertexData[v1][0]+vertexData[v2][0]);
      midWall[1]=0.5*(vertexData[v1][1]+vertexData[v2][1]);
     
      // random rotation angle  
      std::srand (time(NULL));
      double teta=(parameter(0)*3.1415/180)*(1-2*((double) rand() / (RAND_MAX)));

      std::vector<std::vector<double> > rot(2);
      rot[0].resize(2);
      rot[1].resize(2);
      rot[0][0]=std::cos(teta);
      rot[0][1]=-std::sin(teta);
      rot[1][0]=std::sin(teta);
      rot[1][1]=std::cos(teta);
      // vertices of the cell

      // do the rotation for mt_vector
      std::vector<double> vtmp(2); 
      vtmp[0]=cellData[cellInd][mtIndex];
      vtmp[1]=cellData[cellInd][mtIndex+1];
      
      cellData[cellInd][mtIndex]  =rot[0][0]*vtmp[0]+rot[0][1]*vtmp[1];
      cellData[cellInd][mtIndex+1]=rot[1][0]*vtmp[0]+rot[1][1]*vtmp[1];
      
      for (size_t verInd=0; verInd< numWalls; verInd++){
        size_t vi=T.cell(cellInd).vertex(verInd)->index();
        
        
        vtmp[0]=vertexData[vi][0];
        vtmp[1]=vertexData[vi][1];
        // moving the origin to the middle of wall M
        vtmp[0]=vtmp[0]-midWall[0];
        vtmp[1]=vtmp[1]-midWall[1];
        
        // do the rotation for vertices
        vertexData[vi][0]=rot[0][0]*vtmp[0]+rot[0][1]*vtmp[1];
        vertexData[vi][1]=rot[1][0]*vtmp[0]+rot[1][1]*vtmp[1];
       
        
      // moving the origin back
        vertexData[vi][0]= vertexData[vi][0]+midWall[0];
        vertexData[vi][1]= vertexData[vi][1]+midWall[1];
      }
    }
   
  }



    
}


////////////////////////////////
VertexNoUpdateBoundary::
VertexNoUpdateBoundary(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()) {
    std::cerr << "VertexNoUpdateBoundary::"
	      << "VertexNoUpdateBoundary() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=0 && indValue.size()!=1) {
    std::cerr << "VertexNoUpdateBoundary::"
	      << "VertexNoUpdateBoundary() "
	      << "Either no variable indices is used (all directions fixed), " 
	      << "or fixed directions are given in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexNoUpdateBoundary");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}

void VertexNoUpdateBoundary::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t dimension = vertexData[0].size();
  if (numVariableIndexLevel())
    for (size_t dI=0; dI<numVariableIndex(0); ++dI)
      assert(variableIndex(0,dI)<dimension);

  for (size_t i=0; i<numVertices; ++i)
    if (T.vertex(i).isBoundary(T.background())) {
      if (numVariableIndexLevel()) {
	for (size_t dI=0; dI<numVariableIndex(0); ++dI)
	  vertexDerivs[i][variableIndex(0,dI)] = 0.0;
      }
      else {
	for (size_t d=0; d<dimension; ++d)
	  vertexDerivs[i][d] = 0.0;
      }
    }
}



VertexNoUpdateBoundaryPtemplate::  // BB
VertexNoUpdateBoundaryPtemplate(std::vector<double> &paraValue, 
                                std::vector< std::vector<size_t> > 
                                &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()) {
    std::cerr << "VertexNoUpdateBoundaryPtemplate::"
	      << "VertexNoUpdateBoundaryPtemplate() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=1 ) {
    std::cerr << "VertexNoUpdateBoundaryPtemplate::"
	      << "VertexNoUpdateBoundaryPtemplate() "
              << "Start of additional Cell variable indices (center(x,y,z) at first level "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexNoUpdateBoundaryPtemplate");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}



void VertexNoUpdateBoundaryPtemplate::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t comIndex = variableIndex(0,0);  
  size_t numVertices = T.numVertex();
  size_t numWalls = T.numWall();
  size_t dimension = vertexData[0].size();
  
  for (size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex)//check all the vertices
    if (T.vertex(vertexIndex).isBoundary(T.background())){          //if it is at the boundary
      size_t numVertexWalls=T.vertex(vertexIndex).numWall();        //take the number of walls connected to it
      for(size_t wallIndexVertex=0; wallIndexVertex<numVertexWalls; 
          wallIndexVertex++){                                       //for each of those walls
        size_t wallIndex=T.vertex(vertexIndex).wall(wallIndexVertex)->index(); // take global index
        if (T.wall(wallIndex).cell1()==T.background() 
            || T.wall(wallIndex).cell2()==T.background()){          // if the wall is boundary
          size_t cellIndex;


          if(T.wall(wallIndex).cell1()==T.background())             //find the cell which is not background
            cellIndex=T.wall(wallIndex).cell2()->index();
          else
            cellIndex=T.wall(wallIndex).cell1()->index();
          
          // std::cerr<<"com index is   "<<comIndex<<
          //   "cell index is "<<cellIndex<<"  here...............";
          // std::cerr<<cellData[cellIndex][comIndex]<<"  "
          //          <<cellData[cellIndex][comIndex+1]<<"  "
          //          <<cellData[cellIndex][comIndex+2]<<std::endl;
          double COM[3]={cellData[cellIndex][comIndex],
                         cellData[cellIndex][comIndex+1],
                         cellData[cellIndex][comIndex+2]};          // take the position of COM

          
          double CV[3]={vertexData[vertexIndex][0]-COM[0],          
                        vertexData[vertexIndex][1]-COM[1],
                        vertexData[vertexIndex][2]-COM[2]};         // vector vertex-COM
          size_t v1 = T.wall(wallIndex).vertex1()->index();
 	  size_t v2 = T.wall(wallIndex).vertex2()->index();
          double wallVector[3]={vertexData[v1][0]-vertexData[v2][0],
                                vertexData[v1][1]-vertexData[v2][1],
                                vertexData[v1][2]-vertexData[v2][2] };// wall vector
          double temp=std::sqrt(wallVector[0]*wallVector[0]
                                +wallVector[1]*wallVector[1]
                                +wallVector[2]*wallVector[2]);
          
          if(temp<0.00000000001)
            std::cerr<<"VertexNoUpdateBoundaryPtemplate::derivs(), strange wall length at the boundary";
          else
            for(size_t i=0;i<3;i++)                                 // normalize wall vector
              wallVector[i]/=temp;   
          temp=CV[0]*wallVector[0]+CV[1]*wallVector[1]+CV[2]*wallVector[2];
          double edgeNormal[3]={CV[0]-temp*wallVector[0],
                                CV[1]-temp*wallVector[1],
                                CV[2]-temp*wallVector[2]};          // edge normal vector
          temp=std::sqrt(edgeNormal[0]*edgeNormal[0]+
                         edgeNormal[1]*edgeNormal[1]+
                         edgeNormal[2]*edgeNormal[2]);
          
          if(temp<0.00000000001)
            std::cerr<<"VertexNoUpdateBoundaryPtemplate::derivs(), strange edge Normal at the boundary";
          else
            for(size_t i=0;i<3;i++)                             // normalize edge normal vector
              edgeNormal[i]/=temp;   
          
          temp=vertexDerivs[vertexIndex][0]*edgeNormal[0]+      // projection of vertex derivs on edge normal
            vertexDerivs[vertexIndex][1]*edgeNormal[1]+
            vertexDerivs[vertexIndex][2]*edgeNormal[2];
          for (size_t d=0; d<dimension; ++d)
            vertexDerivs[vertexIndex][d] -=temp*edgeNormal[d] ;
          
        }
      }
      // take the list of wall-indices of the walls connected to vertexIndex
      // for(j in the list)
      //     if ( cell1 in background OR cell2 in background)
      //           take the coordinates of the cell-centerCOM which is not background
      //               calculate the vector Vertex-COM
      //               calculate wallvector
      //               calculate A=(vertex-COM)-((wallvector).(Vertex-COM))wallvector/norm(wallvector)
      //               normalize A and make sure it is outward
      //               calculate B=(VderivVector.A)A
      //               VderivVector-=B
      
    }
  
  
  
}

VertexNoUpdateBoundaryPtemplateStatic::  // BB
VertexNoUpdateBoundaryPtemplateStatic(std::vector<double> &paraValue, 
                                std::vector< std::vector<size_t> > 
                                &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()) {
    std::cerr << "VertexNoUpdateBoundaryPtemplateStatic::"
	      << "VertexNoUpdateBoundaryPtemplateStatic() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=1 ) {
    std::cerr << "VertexNoUpdateBoundaryPtemplateStatic::"
	      << "VertexNoUpdateBoundaryPtemplateStatic() "
              << "Start of Cell COM indices at first level "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexNoUpdateBoundaryPtemplateStatic");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}

void VertexNoUpdateBoundaryPtemplateStatic::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
{
  size_t comIndex = variableIndex(0,0);  
  size_t numVertices = T.numVertex();
  size_t dimension = vertexData[0].size();
  
  numBoundaryVertices=0;
  for (size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex)//check all the vertices

    if (T.vertex(vertexIndex).isBoundary(T.background())){          //if it is at the boundary
      boundaryVertices.push_back(vertexIndex);                      // place its index in a vector
      numBoundaryVertices ++;                                          // number of boundary vertices
      size_t numVertexWalls=T.vertex(vertexIndex).numWall();        //take the number of walls connected to it
      //std::cerr<<"boundary vertex "<< vertexIndex  <<" has "<< numVertexWalls << " walls"<< std::endl;
      std::vector<double> cellNormal(3);
      boundaryNormal.push_back(cellNormal);                         // alocate space for the normal vector
      size_t numVertexWallBoundary=0;     // counter for number of boundary walls of the vertex
      for(size_t wallIndexVertex=0; wallIndexVertex<numVertexWalls; 
          wallIndexVertex++){                                       //for each of those walls
        size_t wallIndex=T.vertex(vertexIndex).wall(wallIndexVertex)->index(); // take global index
        
        // std::cerr<<" wall "<< wallIndex  <<" is shared by  "
        //          << T.wall(wallIndex).cell1()->index() << " and "
        //          << T.wall(wallIndex).cell2()->index() << std::endl;
        
        if (T.wall(wallIndex).cell1()==T.background() 
            || T.wall(wallIndex).cell2()==T.background()){          // if the wall is boundary
          numVertexWallBoundary++;
          size_t cellIndex;
          if(T.wall(wallIndex).cell1()==T.background())             //find its cell (which is not background)
            cellIndex=T.wall(wallIndex).cell2()->index();
          else
            cellIndex=T.wall(wallIndex).cell1()->index();
          // std::cerr<<" boundary wall "<< wallIndex  <<" is in the cell "<< cellIndex
          //          << std::endl;

          double COM[3]={cellData[cellIndex][comIndex],
                         cellData[cellIndex][comIndex+1],
                         cellData[cellIndex][comIndex+2]};          // take the position of COM

          double CV[3]={vertexData[vertexIndex][0]-COM[0],          
                        vertexData[vertexIndex][1]-COM[1],
                        vertexData[vertexIndex][2]-COM[2]};         // vector vertex-COM

          size_t v1 = T.wall(wallIndex).vertex1()->index();
 	  size_t v2 = T.wall(wallIndex).vertex2()->index();
          double wallVector[3]={vertexData[v1][0]-vertexData[v2][0],
                                vertexData[v1][1]-vertexData[v2][1],
                                vertexData[v1][2]-vertexData[v2][2] };// wall vector
          // extract cell normal from  (Vertec-COM)x(wallVector)
          cellNormal[0]=CV[1]*wallVector[2]-CV[2]*wallVector[1]; 
          cellNormal[1]=CV[2]*wallVector[0]-CV[0]*wallVector[2];
          cellNormal[2]=CV[0]*wallVector[1]-CV[1]*wallVector[0];

          double temp=std::sqrt(cellNormal[0]*cellNormal[0]
                                +cellNormal[1]*cellNormal[1]
                                +cellNormal[2]*cellNormal[2]);
          
          if(temp<0.00000000001)
            std::cerr<<"VertexNoUpdateBoundaryPtemplateStatic::initiate(), "
                     <<"strange cell normal at the boundary"<<std::endl;
          else
            for(size_t d=0;d<dimension;d++)                                 // normalize cell normal
              cellNormal[d]/=temp;   

          if ( numVertexWallBoundary>1){
            temp=boundaryNormal[numBoundaryVertices-1][0]*cellNormal[0]+
              boundaryNormal[numBoundaryVertices-1][1]*cellNormal[1]+
              boundaryNormal[numBoundaryVertices-1][2]*cellNormal[2];
            if(temp>0)
              for(size_t d=0;d<dimension;d++)  
                boundaryNormal[numBoundaryVertices-1][d]+=cellNormal[d];
            else
              for(size_t d=0;d<dimension;d++)  
                boundaryNormal[numBoundaryVertices-1][d]-=cellNormal[d];
          }
          else
            // add the cell normal to the boundaryNormal
            for(size_t d=0;d<dimension;d++)  
              boundaryNormal[numBoundaryVertices-1][d]+=cellNormal[d];
          
          
        }
        //std::cerr<< numVertexWallBoundary<<std::endl;
        // if(numVertexWallBoundary==0) std::cerr<<"here.................................."<<numBoundaryVertices<<std::endl;;  
        
      }

      if(numVertexWallBoundary!=0){
        double temp=std::sqrt(boundaryNormal[numBoundaryVertices-1][0]
                              *boundaryNormal[numBoundaryVertices-1][0]+
                              boundaryNormal[numBoundaryVertices-1][1]
                              *boundaryNormal[numBoundaryVertices-1][1]+
                              boundaryNormal[numBoundaryVertices-1][2]
                              *boundaryNormal[numBoundaryVertices-1][2]);
        for (size_t d=0; d<dimension; ++d)  // normalize the boundaryNormal
          if (temp!=0) boundaryNormal[numBoundaryVertices-1][d]/=temp;
      }        
      else
        {for (size_t d=0; d<dimension; ++d)  // normalize the boundaryNormal
            boundaryNormal[numBoundaryVertices-1][d]=0;
          std::cerr<< "strange boundary normal"<<std::endl;
        }
      
    }
  
  // std:: cerr<<"number of boundary events  "<<numBoundaryVertices<<std::endl;
  // for (size_t vertex=0; vertex<numBoundaryVertices; vertex++){
  //   std:: cerr<<boundaryVertices[vertex]<<"  "
  //             <<boundaryNormal[vertex][0]<<"  "
  //             <<boundaryNormal[vertex][1]<<"  "
  //             <<boundaryNormal[vertex][2]<<"  "
  //             <<std::endl;
  //}
  
}


void VertexNoUpdateBoundaryPtemplateStatic::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  // std::cerr<< numBoundaryVertices <<"  " 
   //          << boundaryVertices.size() <<"  "
   //          << boundaryNormal.size()<<std::endl;
  
  size_t dimension = vertexData[0].size();
  
  // std:: cerr<<"derives.................................... "
  for (size_t vertex=0; vertex<numBoundaryVertices; vertex++){//for the boundary vertices
    size_t N=boundaryVertices[vertex];
    double norLength=std::abs(boundaryNormal[vertex][0])+
      std::abs(boundaryNormal[vertex][1])+
      std::abs(boundaryNormal[vertex][2]);
    if(norLength>0.1){
      double temp=vertexDerivs[N][0]*boundaryNormal[vertex][0]
        +vertexDerivs[N][1]*boundaryNormal[vertex][1]
        +vertexDerivs[N][2]*boundaryNormal[vertex][2];    
      // projecting vertex derivs on boundary normal
      
      for (size_t d=0; d<dimension; ++d)
        vertexDerivs[N][d]=temp*boundaryNormal[vertex][d];
    }
    else std::cerr<<"strange normal length....................." <<std::endl;
  }
  
}

VertexNoUpdateBoundaryPtemplateStatic3D::  // BB
VertexNoUpdateBoundaryPtemplateStatic3D(std::vector<double> &paraValue, 
                                std::vector< std::vector<size_t> > 
                                &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()) {
    std::cerr << "VertexNoUpdateBoundaryPtemplateStatic3D::"
	      << "VertexNoUpdateBoundaryPtemplateStatic3D() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=1 ) {
    std::cerr << "VertexNoUpdateBoundaryPtemplateStatic3D::"
	      << "VertexNoUpdateBoundaryPtemplateStatic3D() "
              << "Start of Cell COM indices at first level "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexNoUpdateBoundaryPtemplateStatic3D");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}

void VertexNoUpdateBoundaryPtemplateStatic3D::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
{
  size_t comIndex = variableIndex(0,0);  
  size_t numVertices = T.numVertex();
  size_t numCells = T.numCell();
  size_t dimension = vertexData[0].size();
  size_t neighborIndex=27;
  numBottomCells=0;
  numSideCells=0;
  std::vector<std::vector<double> >  cellNormalsBottom, cellNormalsSide;
  for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){//check all the cells for the boundary 

    if (cellData[cellIndex][neighborIndex]==-5){ //if it is at the side
      numSideCells+=1;
      size_t numCellVertices = T.cell(cellIndex).numVertex();
      std::vector<double> normal(3);
      std::vector<double> com(3);
      com=T.cell(cellIndex).positionFromVertex(vertexData);
      for (size_t verInd=0; verInd<numCellVertices; verInd++){
        size_t verIndPlus=verInd+1;
        if(verIndPlus==numCellVertices)
          verIndPlus=0;
        size_t verGInd=T.cell(cellIndex).vertex(verInd) -> index();
        size_t verGIndPlus=T.cell(cellIndex).vertex(verIndPlus) -> index();
        normal[0]+=
          (vertexData[verGInd][1]-com[1])*(vertexData[verGIndPlus][2]-com[2])-
          (vertexData[verGInd][2]-com[2])*(vertexData[verGIndPlus][1]-com[1]);
        normal[1]+=
          (vertexData[verGInd][2]-com[2])*(vertexData[verGIndPlus][0]-com[0])-
          (vertexData[verGInd][0]-com[0])*(vertexData[verGIndPlus][2]-com[2]);
        normal[2]+=
          (vertexData[verGInd][0]-com[0])*(vertexData[verGIndPlus][1]-com[1])-
          (vertexData[verGInd][1]-com[1])*(vertexData[verGIndPlus][0]-com[0]);

      }
      double tmp=std::sqrt(normal[0]*normal[0]+
                           normal[1]*normal[1]+
                           normal[2]*normal[2]);
      
      sideNormals.resize(numSideCells);
      sideNormals[numSideCells-1].resize(4);
      sideNormals[numSideCells-1][3]=cellIndex;
      
      sideNormals[numSideCells-1][0]=normal[0]/tmp;
      sideNormals[numSideCells-1][1]=normal[1]/tmp;
      sideNormals[numSideCells-1][2]=normal[2]/tmp;

    }

    if (cellData[cellIndex][neighborIndex]==-10){ //if it is at the bottom
      numBottomCells+=1;
      size_t numCellVertices = T.cell(cellIndex).numVertex();
      std::vector<double> normal(3);
      std::vector<double> com(3);
      com=T.cell(cellIndex).positionFromVertex(vertexData);
      for (size_t verInd=0; verInd<numCellVertices; verInd++){
        size_t verIndPlus=verInd+1;
        if(verIndPlus==numCellVertices)
          verIndPlus=0;
        size_t verGInd=T.cell(cellIndex).vertex(verInd) -> index();
        size_t verGIndPlus=T.cell(cellIndex).vertex(verIndPlus) -> index();
        normal[0]+=
          (vertexData[verGInd][1]-com[1])*(vertexData[verGIndPlus][2]-com[2])-
          (vertexData[verGInd][2]-com[2])*(vertexData[verGIndPlus][1]-com[1]);
        normal[1]+=
          (vertexData[verGInd][2]-com[2])*(vertexData[verGIndPlus][0]-com[0])-
          (vertexData[verGInd][0]-com[0])*(vertexData[verGIndPlus][2]-com[2]);
        normal[2]+=
          (vertexData[verGInd][0]-com[0])*(vertexData[verGIndPlus][1]-com[1])-
          (vertexData[verGInd][1]-com[1])*(vertexData[verGIndPlus][0]-com[0]);

      }
      double tmp=std::sqrt(normal[0]*normal[0]+
                           normal[1]*normal[1]+
                           normal[2]*normal[2]);
      
      bottomNormals.resize(numBottomCells);
      bottomNormals[numBottomCells-1].resize(4);
      bottomNormals[numBottomCells-1][3]=cellIndex;
      
      bottomNormals[numBottomCells-1][0]=normal[0]/tmp;
      bottomNormals[numBottomCells-1][1]=normal[1]/tmp;
      bottomNormals[numBottomCells-1][2]=normal[2]/tmp;

    }
  }  
}

void VertexNoUpdateBoundaryPtemplateStatic3D::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  size_t dimension = vertexData[0].size();
  
  for (size_t sideInd=0; sideInd<numSideCells; sideInd++){//for the side cells
    size_t N=T.cell(sideNormals[sideInd][3]).numVertex();
    for(size_t verInd=0; verInd<N; verInd++){
      size_t verGInd=T.cell(sideNormals[sideInd][3]).vertex(verInd) -> index();
      double temp=
        vertexDerivs[verGInd][0]*sideNormals[sideInd][0]+
        vertexDerivs[verGInd][1]*sideNormals[sideInd][1]+
        vertexDerivs[verGInd][2]*sideNormals[sideInd][2];
      std::vector<double> vecTemp(3);
      vecTemp[0]=temp*sideNormals[sideInd][0];
      vecTemp[1]=temp*sideNormals[sideInd][1];
      vecTemp[2]=temp*sideNormals[sideInd][2];

      vertexDerivs[verGInd][0]-=vecTemp[0];
      vertexDerivs[verGInd][1]-=vecTemp[1];
      vertexDerivs[verGInd][2]-=vecTemp[2];    
    }
  }
  
  for (size_t bottomInd=0; bottomInd<numBottomCells; bottomInd++){//for the bottom cells
    size_t N=T.cell(bottomNormals[bottomInd][3]).numVertex();
    for(size_t verInd=0; verInd<N; verInd++){
      size_t verGInd=T.cell(bottomNormals[bottomInd][3]).vertex(verInd) -> index();
      double temp=
        vertexDerivs[verGInd][0]*bottomNormals[bottomInd][0]+
        vertexDerivs[verGInd][1]*bottomNormals[bottomInd][1]+
        vertexDerivs[verGInd][2]*bottomNormals[bottomInd][2];
      std::vector<double> vecTemp(3);
      vecTemp[0]=temp*bottomNormals[bottomInd][0];
      vecTemp[1]=temp*bottomNormals[bottomInd][1];
      vecTemp[2]=temp*bottomNormals[bottomInd][2];

      vertexDerivs[verGInd][0]-=vecTemp[0];
      vertexDerivs[verGInd][1]-=vecTemp[1];
      vertexDerivs[verGInd][2]-=vecTemp[2];
      
    }  
  }
    
}

VertexNoUpdateBoundary3D::  // BB
VertexNoUpdateBoundary3D(std::vector<double> &paraValue, 
                                std::vector< std::vector<size_t> > 
                                &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()) {
    std::cerr << "VertexNoUpdateBoundary3D::"
	      << "VertexNoUpdateBoundary3D() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=1 || indValue[0].size()!=2 ) {
    std::cerr << "VertexNoUpdateBoundary3D::"
	      << "VertexNoUpdateBoundary3D() "
              << "Start of Cell COM indices and cell-neighborIndex at first level "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexNoUpdateBoundary3D");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}

void VertexNoUpdateBoundary3D::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
{
  size_t comIndex = variableIndex(0,0);  
  size_t numVertices = T.numVertex();
  size_t numCells = T.numCell();
  size_t dimension = vertexData[0].size();
  size_t neighborIndex=variableIndex(0,1);
  
  size_t numBottomCells=0;
  size_t numSideCells=0;
  std::vector< double > bottomCells, sideCells;
 
  for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){//check all the cells for the boundary 
    
    if (cellData[cellIndex][neighborIndex]==-5){ //if it is at the side
      numSideCells+=1;
      sideCells.resize(numSideCells);
      sideCells[numSideCells-1]=cellIndex;

    }

    if (cellData[cellIndex][neighborIndex]==-10){ //if it is at the bottom
      numBottomCells+=1;
      bottomCells.resize(numBottomCells);
      bottomCells[numBottomCells-1]=cellIndex;
      
    }
    
  }
  
  size_t numBottomVertices=0;
  size_t numSideVertices=0;
  bool isInVector;

  for (size_t cIndex=0; cIndex<numSideCells; ++cIndex){
    size_t cellIndex=sideCells[cIndex];
    size_t numCellVertices = T.cell(cellIndex).numVertex();  

    for (size_t vIndex=0; vIndex<numCellVertices; ++vIndex){
      isInVector=false;
      for (size_t vI=0; vI<numSideVertices; ++vI)
        if(T.cell(cellIndex).vertex(vIndex) -> index() == sideVertices[vI]){
          isInVector=true;
        }

      if(!isInVector){
        numSideVertices++;
        sideVertices.push_back(T.cell(cellIndex).vertex(vIndex) -> index());
      }
    }
  }

  for (size_t cIndex=0; cIndex<numBottomCells; ++cIndex){
    size_t cellIndex=bottomCells[cIndex];
    size_t numCellVertices = T.cell(cellIndex).numVertex();  

    for (size_t vIndex=0; vIndex<numCellVertices; ++vIndex){
      isInVector=false;
      for (size_t vI=0; vI<numBottomVertices; ++vI)
        if(T.cell(cellIndex).vertex(vIndex) -> index() == bottomVertices[vI]){
          isInVector=true;
        }

      if(!isInVector){
        numBottomVertices++;
        bottomVertices.push_back(T.cell(cellIndex).vertex(vIndex) -> index());
      }
    }
  }

  size_t N=T.numSisterVertex();
    
  std::vector<std::vector<double>> tmpsisters;
  std::vector<std::vector<double>> sisters;
  tmpsisters.resize(N);
  for(size_t i=0; i<N; ++i){
    tmpsisters[i].resize(3);
    tmpsisters[i][0]=T.sisterVertex(i,0);
    tmpsisters[i][1]=T.sisterVertex(i,1);
    tmpsisters[i][2]=0;
  }
  
  size_t counter=1;
  for(size_t i=0; i<N; ++i)
    if(tmpsisters[i][2]==0){
      tmpsisters[i][2]=counter;
      sisters.resize(counter);
      sisters[counter-1].push_back(tmpsisters[i][0]);
      sisters[counter-1].push_back(tmpsisters[i][1]);
      for(size_t j=i+1; j<N; j++)
        if(tmpsisters[j][2]==0){
          size_t M=sisters[counter-1].size();
          size_t ww=0;
          for(size_t k=0; k<M; ++k){
            if(tmpsisters[j][0]==sisters[counter-1][k])
              ww+=1;
            if(tmpsisters[j][1]==sisters[counter-1][k])
              ww+=2;
          }
          if(ww==1)
            sisters[counter-1].push_back(tmpsisters[j][1]);
          if(ww==2)
            sisters[counter-1].push_back(tmpsisters[j][0]);
          if(ww!=0)
            tmpsisters[j][2]=counter;
        }
      counter++;
    }
  
  size_t NN = sisters.size();
  for (size_t vI=0; vI<numSideVertices; ++vI){

    for (size_t vs=0; vs<NN; ++vs)
      for (size_t j=0; j<sisters[vs].size(); ++j)
        if(sideVertices[vI]==sisters[vs][j]){

          for (size_t js=0; js<sisters[vs].size(); ++js)
            sideVertices.push_back(sisters[vs][js]);
        }
  }

}

void VertexNoUpdateBoundary3D::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  size_t dimension = vertexData[0].size();
  size_t numSideVertices=sideVertices.size();
  size_t numBottomVertices=bottomVertices.size();

  // for (size_t sideInd=0; sideInd<numSideVertices; sideInd++){//for the side cells
  //   vertexDerivs[sideVertices[sideInd]][0]=0;
  //   vertexDerivs[sideVertices[sideInd]][1]=0;
  //   vertexDerivs[sideVertices[sideInd]][2]=0;
  // }
  
  for (size_t bottomInd=0; bottomInd<numBottomVertices; bottomInd++){//for the bottom cells
    vertexDerivs[bottomVertices[bottomInd]][0]=0;
    vertexDerivs[bottomVertices[bottomInd]][1]=0;
    vertexDerivs[bottomVertices[bottomInd]][2]=0;
  }
     
}

VertexFromConstStressBoundary::  // BB
VertexFromConstStressBoundary(std::vector<double> &paraValue, 
                                std::vector< std::vector<size_t> > 
                                &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()!=7) {
    std::cerr << "VertexFromConstStressBoundary::"
	      << "VertexFromConstStressBoundary() "
	      << "Uses seven parameters for x-stress, y-stress, vertex sencitivity," 
              << "and initial locations "
              << "of right(x), left(x), top(y) and bottom(y) boundaries."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=0) {
    std::cerr << "VertexFromConstStressBoundary::"
	      << "VertexFromConstStressBoundary() "
              << "no index "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexFromConstStressBoundary");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  
  tmp[0] = "stress_x";
  tmp[1] = "stress_y";
  tmp[2] = "vertex_sencitivity";
  tmp[3] = "right_x";
  tmp[4] = "left_x";
  tmp[5] = "top_y";
  tmp[6] = "bottom_y";
  setParameterId( tmp );
}

void VertexFromConstStressBoundary::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
{

  size_t numVertices = T.numVertex();
  size_t dimension = vertexData[0].size();  
  double epcilon=parameter(2);
 
  double rx=parameter(3);
  double lx=parameter(4);
  double ty=parameter(5);
  double by=parameter(6);
  
  numOldVertices=numVertices;  

  rightVertices.resize(3);
  leftVertices.resize(3);
  topVertices.resize(3);
  bottomVertices.resize(3);
  
  std::vector< std::vector<double> > tmprightVertices, tmpleftVertices, tmptopVertices, tmpbottomVertices;
  tmprightVertices.resize(2);
  tmpleftVertices.resize(2);
  tmptopVertices.resize(2);
  tmpbottomVertices.resize(2);
  
  for (size_t vIndex=0; vIndex<numVertices; ++vIndex){
    if(vertexData[vIndex][0]>rx-epcilon && vertexData[vIndex][0]<rx+epcilon){
      tmprightVertices[0].push_back(vIndex);
      tmprightVertices[1].push_back(vertexData[vIndex][1]);
    }
    if(vertexData[vIndex][0]>lx-epcilon && vertexData[vIndex][0]<lx+epcilon){
      tmpleftVertices[0].push_back(vIndex);
      tmpleftVertices[1].push_back(vertexData[vIndex][1]);
    }
    if(vertexData[vIndex][1]>ty-epcilon && vertexData[vIndex][1]<ty+epcilon){
      tmptopVertices[0].push_back(vIndex);
      tmptopVertices[1].push_back(vertexData[vIndex][0]);
    }
    if(vertexData[vIndex][1]>by-epcilon && vertexData[vIndex][1]<by+epcilon){
      tmpbottomVertices[0].push_back(vIndex);
      tmpbottomVertices[1].push_back(vertexData[vIndex][0]);
    }
  }

  // sorting
  double tmp=0;
   for(size_t i=0; i<tmprightVertices[0].size(); ++i)
     for(size_t j=0; j<tmprightVertices[0].size()-1; ++j)
       if(tmprightVertices[1][j]<tmprightVertices[1][j+1]){
         tmp=tmprightVertices[1][j];
         tmprightVertices[1][j]=tmprightVertices[1][j+1];
         tmprightVertices[1][j+1]=tmp;
         
         tmp=tmprightVertices[0][j];
         tmprightVertices[0][j]=tmprightVertices[0][j+1];
         tmprightVertices[0][j+1]=tmp;
       }
   
   for(size_t i=0; i<tmpleftVertices[0].size(); ++i)
     for(size_t j=0; j<tmpleftVertices[0].size()-1; ++j)
       if(tmpleftVertices[1][j]<tmpleftVertices[1][j+1]){
         tmp=tmpleftVertices[1][j];
         tmpleftVertices[1][j]=tmpleftVertices[1][j+1];
         tmpleftVertices[1][j+1]=tmp;
         
         tmp=tmpleftVertices[0][j];
         tmpleftVertices[0][j]=tmpleftVertices[0][j+1];
         tmpleftVertices[0][j+1]=tmp;
       }
   
   for(size_t i=0; i<tmptopVertices[0].size(); ++i)
     for(size_t j=0; j<tmptopVertices[0].size()-1; ++j)
       if(tmptopVertices[1][j]<tmptopVertices[1][j+1]){
         tmp=tmptopVertices[1][j];
         tmptopVertices[1][j]=tmptopVertices[1][j+1];
         tmptopVertices[1][j+1]=tmp;
         
         tmp=tmptopVertices[0][j];
         tmptopVertices[0][j]=tmptopVertices[0][j+1];
         tmptopVertices[0][j+1]=tmp;
       }
   
   for(size_t i=0; i<tmpbottomVertices[0].size(); ++i)
     for(size_t j=0; j<tmpbottomVertices[0].size()-1; ++j)
       if(tmpbottomVertices[1][j]<tmpbottomVertices[1][j+1]){
         tmp=tmpbottomVertices[1][j];
         tmpbottomVertices[1][j]=tmpbottomVertices[1][j+1];
         tmpbottomVertices[1][j+1]=tmp;
         
         tmp=tmpbottomVertices[0][j];
         tmpbottomVertices[0][j]=tmpbottomVertices[0][j+1];
         tmpbottomVertices[0][j+1]=tmp;
       }

   for(size_t i=0; i<tmprightVertices[0].size(); ++i)
     rightVertices[0].push_back(tmprightVertices[0][i]);

   for(size_t i=0; i<tmpleftVertices[0].size(); ++i)
     leftVertices[0].push_back(tmpleftVertices[0][i]);

   for(size_t i=0; i<tmptopVertices[0].size(); ++i)
     topVertices[0].push_back(tmptopVertices[0][i]);

   for(size_t i=0; i<tmpbottomVertices[0].size(); ++i)
     bottomVertices[0].push_back(tmpbottomVertices[0][i]);

   rightVertices[1].push_back(0);
   leftVertices[1].push_back(0);
   topVertices[1].push_back(0);
   bottomVertices[1].push_back(0);
   
   for(size_t i=1; i<tmprightVertices[0].size(); ++i)
     rightVertices[1].push_back(tmprightVertices[1][i-1]-tmprightVertices[1][i]);
   for(size_t i=1; i<tmpleftVertices[0].size(); ++i)
     leftVertices[1].push_back(tmpleftVertices[1][i-1]-tmpleftVertices[1][i]);
   for(size_t i=1; i<tmptopVertices[0].size(); ++i)
     topVertices[1].push_back(tmptopVertices[1][i-1]-tmptopVertices[1][i]);
   for(size_t i=1; i<tmpbottomVertices[0].size(); ++i)
     bottomVertices[1].push_back(tmpbottomVertices[1][i-1]-tmpbottomVertices[1][i]);

   for(size_t i=0; i<tmprightVertices[0].size()-1; ++i)
     rightVertices[2].push_back(tmprightVertices[1][i]-tmprightVertices[1][i+1]);
   for(size_t i=0; i<tmpleftVertices[0].size()-1; ++i)
     leftVertices[2].push_back(tmpleftVertices[1][i]-tmpleftVertices[1][i+1]);
   for(size_t i=0; i<tmptopVertices[0].size()-1; ++i)
     topVertices[2].push_back(tmptopVertices[1][i]-tmptopVertices[1][i+1]);
   for(size_t i=0; i<tmpbottomVertices[0].size()-1; ++i)
     bottomVertices[2].push_back(tmpbottomVertices[1][i]-tmpbottomVertices[1][i+1]);
   
   rightVertices[2].push_back(0);
   leftVertices[2].push_back(0);
   topVertices[2].push_back(0);
   bottomVertices[2].push_back(0);

   // std::cerr<<std::endl;
   // for(size_t i=0; i<tmprightVertices[0].size(); ++i)
   //   std::cerr<<rightVertices[0][i]<<"  ";
   // std::cerr<<std::endl;   
   // for(size_t i=0; i<tmpleftVertices[0].size(); ++i)
   //   std::cerr<<leftVertices[0][i]<<"  ";
   // std::cerr<<std::endl;
   // for(size_t i=0; i<tmptopVertices[0].size(); ++i)
   //   std::cerr<<topVertices[0][i]<<"  ";
   // std::cerr<<std::endl;
   // for(size_t i=0; i<tmpbottomVertices[0].size(); ++i)
   //   std::cerr<<bottomVertices[0][i]<<"  ";
   // std::cerr<<std::endl;

  // rightVertices[1].push_back(rx);
  // rightVertices[1].push_back(0);
  // leftVertices[1].push_back(lx);
  // leftVertices[1].push_back(0);
  // topVertices[1].push_back(ty);
  // topVertices[1].push_back(0);
  // bottomVertices[1].push_back(by);
  // bottomVertices[1].push_back(0);
}

void VertexFromConstStressBoundary::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //double sx=std::exp(0.022*parameter(0)*parameter(0))+1;  
  double sx=parameter(0);  
  double sy=parameter(1);  

  // double l1=0, l2=0;
  // // fix with the boundary
  // for(size_t i=0; i<rightVertices[0].size(); ++i){
  //   if(i==0)
  //     l1=0;
  //   else
  //     l1=vertexData[rightVertices[0][i-1]][1]-vertexData[rightVertices[0][i]][1];
  //   if(i==rightVertices[0].size()-1)
  //     l2=0;
  //   else
  //     l2=vertexData[rightVertices[0][i]][1]-vertexData[rightVertices[0][i+1]][1];
  //   vertexDerivs[rightVertices[0][i]][1]=(l1-rightVertices[1][i])-(l2-rightVertices[2][i]);
  // }
  // for(size_t i=0; i<leftVertices[0].size(); ++i){
  //   if(i==0)
  //     l1=0;
  //   else
  //     l1=vertexData[leftVertices[0][i-1]][1]-vertexData[leftVertices[0][i]][1];
  //   if(i==leftVertices[0].size()-1)
  //     l2=0;
  //   else
  //     l2=vertexData[leftVertices[0][i]][1]-vertexData[leftVertices[0][i+1]][1];
  //   vertexDerivs[leftVertices[0][i]][1]=(l1-leftVertices[1][i])-(l2-leftVertices[2][i]);
  // }

  // for(size_t i=0; i<topVertices[0].size(); ++i){
  //   if(i==0)
  //     l1=0;
  //   else
  //     l1=vertexData[topVertices[0][i-1]][0]-vertexData[topVertices[0][i]][0];
  //   if(i==topVertices[0].size()-1)
  //     l2=0;
  //   else
  //     l2=vertexData[topVertices[0][i]][0]-vertexData[topVertices[0][i+1]][0];
  //   vertexDerivs[topVertices[0][i]][0]=(l1-topVertices[1][i])-(l2-topVertices[2][i]);
  // }
  // for(size_t i=0; i<bottomVertices[0].size(); ++i){
  //   if(i==0)
  //     l1=0;
  //   else
  //     l1=vertexData[bottomVertices[0][i-1]][0]-vertexData[bottomVertices[0][i]][0];
  //   if(i==bottomVertices[0].size()-1)
  //     l2=0;
  //   else
  //     l2=vertexData[bottomVertices[0][i]][0]-vertexData[bottomVertices[0][i+1]][0];
  //   vertexDerivs[bottomVertices[0][i]][0]=(l1-bottomVertices[1][i])-(l2-bottomVertices[2][i]);
  // }

  double deltaY=vertexData[topVertices[0][0]][1]-vertexData[bottomVertices[0][0]][1];
  double deltaX=vertexData[rightVertices[0][0]][0]-vertexData[leftVertices[0][0]][0];
  
  double Fr= sx*deltaY;
  double Fl=-sx*deltaY;
  double Ft= sy*deltaX;
  double Fb=-sy*deltaX;

  for(size_t i=0; i<rightVertices[0].size(); ++i){
    Fr+=vertexDerivs[rightVertices[0][i]][0];
  }
  for(size_t i=0; i<leftVertices[0].size(); ++i){
    Fl+=vertexDerivs[leftVertices[0][i]][0];
  }
  for(size_t i=0; i<topVertices[0].size(); ++i){
    Ft+=vertexDerivs[topVertices[0][i]][1];
  }
  for(size_t i=0; i<bottomVertices[0].size(); ++i){
    Fb+=vertexDerivs[bottomVertices[0][i]][1];
  }  
  
  //std::cout<<  Fr<<" "<< Fl<<" "<< Ft<<" "<< Fb<<" "<<std::endl;

  Fr/= rightVertices[0].size();
  Fl/=  leftVertices[0].size();
  Ft/=   topVertices[0].size();
  Fb/=bottomVertices[0].size();
  
  //if(totaltime<500){ for growthTec the time after which relax the template
  if(true){
    for(size_t i=0; i<rightVertices[0].size(); ++i){
      vertexDerivs[rightVertices[0][i]][0]=Fr;
    }
    for(size_t i=0; i<leftVertices[0].size(); ++i){
      vertexDerivs[leftVertices[0][i]][0]=Fl;
    }
    for(size_t i=0; i<topVertices[0].size(); ++i){
      vertexDerivs[topVertices[0][i]][1]=Ft;
    }
    for(size_t i=0; i<bottomVertices[0].size(); ++i){
      vertexDerivs[bottomVertices[0][i]][1]=Fb;
    }
        
    // keeping on xy plane
    for(size_t i=0; i<rightVertices[0].size(); ++i){
      vertexDerivs[rightVertices[0][i]][2]=0 ;
    }
    for(size_t i=0; i<leftVertices[0].size(); ++i){
      vertexDerivs[leftVertices[0][i]][2]=0 ;
    }
    for(size_t i=0; i<topVertices[0].size(); ++i){
      vertexDerivs[topVertices[0][i]][2]=0 ;
    }
    for(size_t i=0; i<bottomVertices[0].size(); ++i){
      vertexDerivs[bottomVertices[0][i]][2]=0 ;
    }
    
  }
  
}

void VertexFromConstStressBoundary::
update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h)
{
  size_t numVertices = T.numVertex();
  size_t dimension = vertexData[0].size();  
  double epcilon=parameter(2);
  static double tt=0;
  tt+=h;
  totaltime=tt;

  double rx=vertexData[rightVertices[0][0]][0];
  double lx=vertexData[leftVertices[0][0]][0];
  double ty=vertexData[topVertices[0][0]][1];
  double by=vertexData[bottomVertices[0][0]][1];
  
  // adding new vertices from division  
  if(numVertices>numOldVertices){
    for (size_t vIndex=numOldVertices; vIndex<numVertices; ++vIndex){
      std::cerr<<"vertex added to the boundary";
      if(vertexData[vIndex][0]>rx-epcilon && vertexData[vIndex][0]<rx+epcilon)
        rightVertices[0].push_back(vIndex);
      if(vertexData[vIndex][0]>lx-epcilon && vertexData[vIndex][0]<lx+epcilon)
        leftVertices[0].push_back(vIndex);
      if(vertexData[vIndex][1]>ty-epcilon && vertexData[vIndex][1]<ty+epcilon)
        topVertices[0].push_back(vIndex);
      if(vertexData[vIndex][1]>by-epcilon && vertexData[vIndex][1]<by+epcilon)
        bottomVertices[0].push_back(vIndex);
    }
  }
  
  numOldVertices=numVertices;
}

manipulate::  // BB
manipulate(std::vector<double> &paraValue, 
           std::vector< std::vector<size_t> > 
           &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()!=0) {
    std::cerr << "manipulate::"
	      << "manipulate() "
	      << "Uses no parameter" 
 	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=0) {
    std::cerr << "manipulate::"
	      << "manipulate() "
              << "no index "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("manipulate");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}

void manipulate::
initiate(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  size_t numVertices = T.numVertex();
  size_t numCells = T.numCell();
  size_t numWalls = T.numWall();
  
  size_t dimension = vertexData[0].size();  
  size_t centerIndex=35;
  double xx=5;
  double yy=5;

  for (size_t i=0; i<numCells; ++i){
    if(cellData[i][34]==0){
      cellData[i][0]=1;
      cellData[i][1]=0;
      cellData[i][2]=0;
      
    }
    if(cellData[i][34]==7){
      cellData[i][0]=0;
      cellData[i][1]=1;
      cellData[i][2]=0;
    }
    if(cellData[i][34]!=0 && cellData[i][34]!=7 ){
      cellData[i][0]=cellData[i][centerIndex]-xx;
      cellData[i][1]=cellData[i][centerIndex+1]-yy;
      cellData[i][2]=0;    
      double tmp=std::sqrt(cellData[i][0]*cellData[i][0]+cellData[i][1]*cellData[i][1]);
      cellData[i][0]/=tmp;
      cellData[i][1]/=tmp;
    }
  }
  // for (size_t i=0; i<numWalls; ++i)
  //   std::cout<<wallData[i][0]<<std::endl;
  

  // for (size_t i=0; i<numVertices; ++i){

  //     vertexData[i][1]*=0.3;
  // }

  // for (size_t i=0; i<numCells; ++i){
  //   double Radius=std::sqrt((cellData[i][centerIndex]-xx)*(cellData[i][centerIndex]-xx)
  //                           +(cellData[i][centerIndex+1]-yy)*(cellData[i][centerIndex+1]-yy));
  //   if(Radius<.1)
  //     cellData[i][40]=100;
  // }
  
  // double l=0.1;
  // double r=1;
  // double c=3;

  // for (size_t i=0; i<numCells; ++i){
  //   cellData[i][0]=1;
  //   cellData[i][1]=0;
  //   cellData[i][2]=0;
  //   cellData[i][3]=((r-l)/(2*c))*cellData[i][centerIndex]+((r+l)/2);

  //   cellData[i][4]=0;
  //   cellData[i][5]=1;
  //   cellData[i][6]=0;
  //   cellData[i][7]=((l-r)/(2*c))*cellData[i][centerIndex]+((r+l)/2);
  // }
}

void manipulate::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 

{}

void manipulate::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h)
{
  size_t numVertices = T.numVertex();
  size_t numCells = T.numCell();

  size_t dimension = vertexData[0].size();  
  size_t centerIndex=39;

  // double epcilon=.1;
  // double alpha=3.1415/2;

  // // double A=0.4;
  // double A=12;// fiber strength
  
  // size_t fiberIndex=21;
  // size_t fiberLIndex=16;

  // for (size_t i=0; i<numCells; ++i)
  //   if(cellData[i][centerIndex]>0){
  //     cellData[i][fiberIndex]=A*0.5*(1+std::cos(alpha*cellData[i][centerIndex+1]));
  //     cellData[i][fiberLIndex]=cellData[i][fiberIndex]/2;
  //   }
  //   else{
  //     cellData[i][fiberIndex]=A*0.5*(1-std::cos(alpha*cellData[i][centerIndex+1]));
  //     cellData[i][fiberLIndex]=cellData[i][fiberIndex]/2;
  //   }
  
  // for (size_t vIndex=0; vIndex<numVertices; ++vIndex)
  //   if(std::abs(vertexData[vIndex][0])<epcilon){
     
  //     vertexData[vIndex][0]=A*std::cos(alpha*vertexData[vIndex][1]);
  //   }
}

cellPolarity3D::  // BB
cellPolarity3D(std::vector<double> &paraValue, 
               std::vector< std::vector<size_t> > 
               &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()) {
    std::cerr << "cellPolarity3D::"
	      << "cellPolarity3D() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=2 ) {
    std::cerr << "cellPolarity3D::"
	      << "cellPolarity3D() "
              << "cell_identity_index and start of Cell COM index at first level "
              << "polarity_vector_index and polarity_measure_index at second level "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("cellPolarity3D");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}

void cellPolarity3D::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
{
  size_t cellId = variableIndex(0,0);  
    
  //size_t numVertices = T.numVertex();
  size_t numCells = T.numCell();
  size_t dimension = vertexData[0].size();

  size_t num3dCells=0;
  for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex)//find the number of 3d cells
    if (cellData[cellIndex][cellId]>num3dCells)
      num3dCells=cellData[cellIndex][cellId];  
  
  cellFaces.resize(num3dCells+1);
  cellCentPol.resize(num3dCells+1);

  for (size_t CellInd3d=0; CellInd3d<num3dCells+1; ++CellInd3d){
    for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex)
      if (cellData[cellIndex][cellId]==CellInd3d) 
        cellFaces[CellInd3d].push_back(cellIndex);
    cellCentPol[CellInd3d].resize(6);
  }
  
  // size_t numtemp=cellFaces.size();
  // std::cerr<<"number of real cells "<<numtemp<<std::endl;
  // for(size_t i=0; i< numtemp; ++i){
  //   size_t N=cellFaces[i].size();
  //   for(size_t j=0; j<N ; ++j)
  //     std::cerr<<cellFaces[i][j]<<" ";
  //   std::cerr<<std::endl;
  // }
}

void cellPolarity3D::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // nothing is needed here  
}

void cellPolarity3D::
update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h)
{
  // polarity vector calculation
  // very ad-hoc for marcus , summs the normal vectors of 
  // cell walls scaled with misses stress of the wall. Not good in 
  // general but works if the cell walls are almost symmetric-.
  
  size_t num3dCells=cellCentPol.size();
  size_t comIndex = variableIndex(0,1);  
  size_t wallId=28;
  size_t Lid=29;
  size_t polVecInd = variableIndex(1,0);  
  //size_t polMesInd = 11;
  size_t polMesInd = variableIndex(1,1);
   
  double area=0;
  double totalarea=0;
  std::vector<double> tmpCentPol(6);
  for (size_t CellInd3d=0; CellInd3d<num3dCells; ++CellInd3d) {
    
    totalarea=0;
    for (size_t i=0; i<6; ++i)
      tmpCentPol[i]=0;
    
    if (cellData[cellFaces[CellInd3d][0]][Lid]==1){// if 3dCell is in L1
      size_t numWalls=cellFaces[CellInd3d].size();
      for (size_t wallIndex=0; wallIndex<numWalls; ++wallIndex)
        if (cellData[cellFaces[CellInd3d][wallIndex]][wallId]==0) { // if the wall is anticlinal 
          
          area=T.cell(cellFaces[CellInd3d][wallIndex]).calculateVolume(vertexData);    
          totalarea+=area;          

          tmpCentPol[0]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex  ]*area;
          tmpCentPol[1]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex+1]*area;
          tmpCentPol[2]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex+2]*area;
          double tmp=cellData[cellFaces[CellInd3d][wallIndex]][polMesInd]*area;
         
          tmpCentPol[3]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd  ];
          tmpCentPol[4]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd+1];
          tmpCentPol[5]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd+2];
        }
    }
    
    if (cellData[cellFaces[CellInd3d][0]][Lid]==2 || 
        cellData[cellFaces[CellInd3d][0]][Lid]==3){// if 3dCell is in L2 or L3
      size_t numWalls=cellFaces[CellInd3d].size();
      for (size_t wallIndex=0; wallIndex<numWalls; ++wallIndex) { // for all the cell walls 
         
          area=T.cell(cellFaces[CellInd3d][wallIndex]).calculateVolume(vertexData);    
          totalarea+=area;          

          tmpCentPol[0]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex  ]*area;
          tmpCentPol[1]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex+1]*area;
          tmpCentPol[2]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex+2]*area;
          double tmp=cellData[cellFaces[CellInd3d][wallIndex]][polMesInd]*area;

          tmpCentPol[3]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd  ];
          tmpCentPol[4]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd+1];
          tmpCentPol[5]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd+2];
        }
    }
    
    if(totalarea != 0){
      tmpCentPol[0]/=totalarea;
      tmpCentPol[1]/=totalarea;
      tmpCentPol[2]/=totalarea;
    }
    cellCentPol[CellInd3d]=tmpCentPol;
    
  }
}

// {
//   // polarity vector calculation
//   // very ad-hoc for marcus , summs the normal vectors of 
//   // cell walls scaled with misses stress of the wall. Not good in 
//   // general but works if the cell walls are almost symmetric-.

//   size_t num3dCells=cellCentPol.size();
//   size_t comIndex = variableIndex(0,1);  
//   size_t wallId=28;
//   size_t Lid=29;
//   size_t polVecInd = variableIndex(1,0);  
//   size_t polMesInd = variableIndex(1,1);

//   size_t numSides;
//   std::vector<double> tmpCentPol(6);
//   for (size_t CellInd3d=0; CellInd3d<num3dCells; ++CellInd3d) {
    
//     numSides=0;
//     for (size_t i=0; i<6; ++i)
//       tmpCentPol[i]=0;
    
//     if (cellData[cellFaces[CellInd3d][0]][Lid]==1){// if 3dCell is in L1
//       size_t numWalls=cellFaces[CellInd3d].size();
//       for (size_t wallIndex=0; wallIndex<numWalls; ++wallIndex)
//         if (cellData[cellFaces[CellInd3d][wallIndex]][wallId]==0) { // if the wall is anticlinal 
//           numSides++;        
//           tmpCentPol[0]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex  ];
//           tmpCentPol[1]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex+1];
//           tmpCentPol[2]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex+2];
//           double tmp=cellData[cellFaces[CellInd3d][wallIndex]][polMesInd]*
//             T.cell(cellFaces[CellInd3d][wallIndex]).calculateVolume(vertexData);
//           tmpCentPol[3]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd  ];
//           tmpCentPol[4]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd+1];
//           tmpCentPol[5]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd+2];
//         }
//     }
    
//     if (cellData[cellFaces[CellInd3d][0]][Lid]==2 || 
//         cellData[cellFaces[CellInd3d][0]][Lid]==3){// if 3dCell is in L2 or L3
//       size_t numWalls=cellFaces[CellInd3d].size();
//       for (size_t wallIndex=0; wallIndex<numWalls; ++wallIndex) { // for all the cell walls 
//           numSides++;        
//           tmpCentPol[0]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex  ];
//           tmpCentPol[1]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex+1];
//           tmpCentPol[2]+=cellData[cellFaces[CellInd3d][wallIndex]][comIndex+2];
//           double tmp=cellData[cellFaces[CellInd3d][wallIndex]][polMesInd]*
//             T.cell(cellFaces[CellInd3d][wallIndex]).calculateVolume(vertexData);
//           tmpCentPol[3]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd  ];
//           tmpCentPol[4]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd+1];
//           tmpCentPol[5]+=tmp*cellData[cellFaces[CellInd3d][wallIndex]][polVecInd+2];
//         }
//     }
    
//     if(numSides != 0){
//       tmpCentPol[0]/=numSides;
//       tmpCentPol[1]/=numSides;
//       tmpCentPol[2]/=numSides;
//     }
//     cellCentPol[CellInd3d]=tmpCentPol;
    
//   }

//   // size_t numtemp=cellCentPol.size();
//   // std::cerr<<"centers and pol vec "<<numtemp<<std::endl;
//   // for(size_t i=0; i< numtemp; ++i){
//   //   size_t N=cellCentPol[i].size();
//   //   for(size_t j=0; j<N ; ++j)
//   //     std::cerr<<cellCentPol[i][j]<<" ";
//   //   std::cerr<<std::endl;
//   // }
// }

void cellPolarity3D::
printState(Tissue *T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData, 
	   std::ostream &os)
{ // VTK style
  
  double RR=12; // pol vectors are shown inside this radius(to exclude boundary)
  size_t numPoints = cellCentPol.size();
  size_t dimension=3;
  static size_t index=0;

  std::stringstream name;
  name << "tmp/VTKPolVec" << index  <<".vtu";
  std::ofstream myfile;
  myfile.open (name.str());

  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<numPoints<<"\" NumberOfCells=\"1\">"<< std::endl;

  myfile << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  for (size_t vIndex=0 ; vIndex<numPoints ; ++vIndex){
    
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << cellCentPol[vIndex][d] << " ";
    myfile << std::endl;
  }
  myfile << "</DataArray>"<<std::endl
         << "</Points>"<<std::endl;

  myfile << "<Cells>"<<std::endl;

  //connectivity
  myfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;
  for (size_t cIndex=0 ; cIndex<numPoints ; cIndex++)
    myfile << cIndex <<" ";
  
  myfile << std::endl
         << "</DataArray>"<<std::endl;
  
  //off-sets 
  myfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  myfile << numPoints ;
  
  myfile << std::endl
         << "</DataArray>"<<std::endl;
  
  //types
  myfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  myfile <<2;
  myfile << std::endl
         << "</DataArray>"<<std::endl;
  
  myfile << "</Cells>"<<std::endl;
  myfile << "<PointData>"<<std::endl;
  myfile << "<DataArray type=\"Float64\" Name=\"PolarityVector\" NumberOfComponents=\"3\" format=\"ascii\">"
         <<std::endl;
  for (size_t vIndex=0 ; vIndex<numPoints ; vIndex++){
    double radius=std::sqrt(cellCentPol[vIndex][0]*cellCentPol[vIndex][0]+
                            cellCentPol[vIndex][1]*cellCentPol[vIndex][1]
                            );
    if(radius<RR){
      for( size_t d=0 ; d<dimension ; d++ )
        myfile << cellCentPol[vIndex][d+3] << " ";
      myfile << std::endl;
    }   
    else{
      for( size_t d=0 ; d<dimension ; d++ )
        myfile << 0 << " ";
      myfile << std::endl;
    }
  }
 
  myfile << "</DataArray>"<<std::endl
         << "</PointData>"<<std::endl 
	 << "</Piece>"<<std::endl
	 << "</UnstructuredGrid>"<<std::endl
	 << "</VTKFile>"<<std::endl;
  myfile.close();

  index++;
}

diffusion3D::  // BB
diffusion3D(std::vector<double> &paraValue, 
                                std::vector< std::vector<size_t> > 
                                &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()!=1) {
    std::cerr << "diffusion3D::"
	      << "diffusion3D() "
	      << "Uses one parameter, diffusion constant."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=1 || indValue[0].size()!=3 ) {
    std::cerr << "diffusion3D::"
	      << "diffusion3D() "
              << "needs 3 indices in one level, concentration_index, "
	      << "neighbor_wall_index and 3Dcell_index at the first level "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("diffusion3D");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_0";
  setParameterId( tmp );
}

void diffusion3D::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
{
  size_t concIndex = variableIndex(0,0);  
  size_t neighIndex= variableIndex(0,1);  
  size_t CellIndex3d = variableIndex(0,2);  

  size_t numCells = T.numCell();
  size_t num3dCells = 0;
 
  // for (size_t cellInd=0; cellInd<numCells; ++cellInd)//change cell variable for 3dCell 432
  //   if (cellData[cellInd][CellIndex3d]==432 ){
  //     cellData[cellInd][31]=10;
  //     cellData[cellInd][35]=10;
  //   }

  for (size_t cellInd=0; cellInd<numCells; ++cellInd){//change cell variable for 3dCell 432
    if (cellData[cellInd][30]>=5){
      cellData[cellInd][30]=10;
    }
    if (cellData[cellInd][30]<5 && cellData[cellInd][30]>=2.5){
      cellData[cellInd][30]=5;
    }
    if (cellData[cellInd][30]<2.5 && cellData[cellInd][30]>1){
      cellData[cellInd][30]=2;
    }
  }
  
  for (size_t cellInd=0; cellInd<numCells; ++cellInd)//check all the cells for number of 3d cells
    if (cellData[cellInd][CellIndex3d]>num3dCells )
      num3dCells=cellData[cellInd][CellIndex3d];
  num3dCells++;
  Cells3d.resize(num3dCells);
  //Cells3d(3dCellInd, compInd)(num3Dcells,1+numWalls+1+numNeighbors) 
  //[numwalls,wall(cell)Index1,...,wall(cell)IndexN,
  // numNeghbors,,]

  for (size_t i=0; i<num3dCells; ++i)// reserve first component for number of cell walls  
    Cells3d[i].push_back(0);

  for (size_t cellInd=0; cellInd<numCells; ++cellInd)//putting wall(cell) indices into the vector
    Cells3d[cellData[cellInd][CellIndex3d]].push_back(cellInd);
    
  for (size_t i=0; i<num3dCells; ++i)// storing the cell wall indices  
    Cells3d[i][0]=Cells3d[i].size()-1;

   for (size_t i=0; i<num3dCells; ++i){
    std::vector<size_t> tmpVec;

    for (size_t j=1; j<Cells3d[i][0]+1; ++j){
         
      int walltmp=cellData[Cells3d[i][j]][neighIndex];
      if(walltmp!=-2 && walltmp!=-5 &&      walltmp!=-7 &&      walltmp!=-10)
	tmpVec.push_back(cellData[walltmp][CellIndex3d]);      
      else
        tmpVec.push_back(-1);// elements at the boundary with no neighbor        
    }
    
    //std::sort(tmpVec.begin(), tmpVec.end() );
    //tmpVec.erase( std::unique( tmpVec.begin(), tmpVec.end() ), tmpVec.end() );
    Cells3d[i].push_back(tmpVec.size());    
    
    for (size_t k=0; k<tmpVec.size(); ++k)
      Cells3d[i].push_back(tmpVec[k]);    
   }
  
  faceArea.resize(num3dCells); 
  for (size_t i=0; i<num3dCells; ++i){
    for (size_t j=1; j<Cells3d[i][0]+1; ++j)
      faceArea[i].push_back(T.cell(Cells3d[i][j]).calculateVolume(vertexData));
  }
  
  for (size_t i=0; i<num3dCells; ++i){
    std::cerr<<i<<"               ";
    for (size_t j=0; j<Cells3d[i][0]; ++j)
      std::cerr<<faceArea[i][j]<<" ";
    std::cerr<<std::endl; 
  }
}

void diffusion3D::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  size_t num3Cells = Cells3d.size();
  size_t aI = variableIndex(0,0);//conc index
  assert( aI<cellData[0].size());
  
  for( size_t i=0 ; i<num3Cells ; ++i ) {
    double cellDer=0;
    for( size_t j=Cells3d[i][0]+2 ; j<Cells3d[i].size() ; ++j) {
      size_t neighInd=Cells3d[i][j];
      if(neighInd !=-1)
        cellDer += parameter(0)*
          ( cellData[Cells3d[neighInd][1]][aI] - cellData[Cells3d[i][1]][aI])*
          faceArea[i][j-Cells3d[i][0]-2];
    }
    for( size_t j=1 ; j<Cells3d[i][0]+1 ; ++j) {
      cellDerivs[Cells3d[i][j]][aI]+=cellDer;
    }
  }  
}

VertexTranslateToMax::
VertexTranslateToMax(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexTranslateToMax::VertexTranslateToMax() "
	      << "Uses one parameter, maxPos " << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexTranslateToMax::"
	      << "VertexTranslateToMax() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexTranslateToMax");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "maxPos";
  setParameterId( tmp );
}

void VertexTranslateToMax::
initiate(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  assert( posIndex < vertexData[0].size() );
  
  double max = vertexData[0][posIndex];
  for( size_t i=1 ; i<numVertices ; ++i )
    if (vertexData[i][posIndex]>max)
      max = vertexData[i][posIndex];
  double delta = max-parameter(0);
  for( size_t i=0 ; i<numVertices ; ++i )
    vertexData[i][posIndex] -= delta;
}

void VertexTranslateToMax::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
}

void VertexTranslateToMax::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
			 double h)
{
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  assert( posIndex < vertexData[0].size() );
  
	double max = vertexData[0][posIndex];
  for( size_t i=1 ; i<numVertices ; ++i )
    if (vertexData[i][posIndex]>max)
			max = vertexData[i][posIndex];
	double delta = max-parameter(0);
  for( size_t i=0 ; i<numVertices ; ++i )
		vertexData[i][posIndex] -= delta;
}

CenterCOM::
CenterCOM(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue) {
  //Do some checks on the parameters and variable indeces
  //
  if (paraValue.size() != 0) {
    std::cerr << "CenterCOM::CenterCOM() Uses no parameters.\n";
    std::exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 0) {
    std::cerr << "CenterCOM::CenterCOM() No variable indices used.\n";
    std::exit(EXIT_FAILURE);
  }
  
  //Set the variable values
  //
  setId("CenterCOM");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp(numParameter());
  setParameterId(tmp);
}

void CenterCOM::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs) 
{
  update(T, cellData, wallData, vertexData, 0.0);
}

void CenterCOM::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs) 
{
  
}

void CenterCOM::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  
}

void CenterCOM::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h)
{
  size_t dimension = vertexData[0].size();
  size_t numVertices = T.numVertex();
  
  std::vector<double> com(dimension, 0.0);
  for (size_t i = 0; i < numVertices; ++i) {
    for (size_t d = 0; d < dimension; ++d) {
      com[d] += vertexData[i][d];
    }
  }
  
  for (size_t d = 0; d < dimension; ++d) {
    com[d] /= numVertices;
  }
  
  for (size_t i = 0; i < numVertices; ++i) {
    for (size_t d = 0; d < dimension; ++d) {
      vertexData[i][d] -= com[d];

    }
  }
}

CenterCOMcenterTriangulation::
CenterCOMcenterTriangulation(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue) {
  //Do some checks on the parameters and variable indeces
  //
  if (paraValue.size() != 0) {
    std::cerr << "CenterCOMcenterTriangulation::CenterCOMcenterTriangulation() Uses no parameters.\n";
    std::exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "CenterCOMcenterTriangulation::CenterCOMcenterTriangulation()"
	      << " Uses one variable index for where cell vertices are stored "
	      << "(must match TRBScenterTriangulation)." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  //Set the variable values
  //
  setId("CenterCOMcenterTriangulation");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp(numParameter());
  setParameterId(tmp);
}

void CenterCOMcenterTriangulation::
initiate(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs) 
{
  update(T, cellData, wallData, vertexData, 0.0);
}

void CenterCOMcenterTriangulation::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs) 
{
  
}

void CenterCOMcenterTriangulation::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h)
{
  size_t dimension = vertexData[0].size();
  size_t numVertices = T.numVertex();
  
  std::vector<double> com(dimension, 0.0);
  
  for (size_t i = 0; i < numVertices; ++i) {
    for (size_t d = 0; d < dimension; ++d) {
      com[d] += vertexData[i][d];
    }
  }
  
  for (size_t d = 0; d < dimension; ++d) {
    com[d] /= numVertices;
  }
  
  for (size_t i = 0; i < numVertices; ++i) {
    for (size_t d = 0; d < dimension; ++d) {
      vertexData[i][d] -= com[d];
    }
  }
  // Also update the additional central vertices stored in the cells
  size_t numCells = T.numCell();
  size_t cellVarIndex = variableIndex(0,0);
  for (size_t i = 0; i < numCells; ++i) {
    for (size_t d = 0; d < dimension; ++d) {
      cellData[i][cellVarIndex+d] -= com[d];
    }
  }
  
}

CalculatePCAPlane::
CalculatePCAPlane(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //
  // Do some checks on the parameters and variable indeces
  //
  if (paraValue.size()!=1) {
    std::cerr << "CalculatePCAPlane::"
	      << "CalculatePCAPlane() "
	      << "Uses one parameter: onlyInUpdateFlag "
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=0) {
    std::cerr << "CalculatePCAPlane::"
	      << "CalculatePCAPlane() "
	      << "No indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("CalculatePCAPlane");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "onlyInUpdateFlag";
  setParameterId( tmp );
}

void CalculatePCAPlane::
initiate(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  if (parameter(0)==1.0) {
    size_t numCell = T.numCell();
    for (size_t i=0; i<numCell; ++i)
      T.cell(i).calculatePCAPlane(vertexData);
  }
}

void CalculatePCAPlane::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  if (parameter(0)!=1.0) {
    size_t numCell = T.numCell();
    for (size_t i=0; i<numCell; ++i)
      T.cell(i).calculatePCAPlane(vertexData);
  }
}

void CalculatePCAPlane::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h)
{
  if (parameter(0)==1.0) {
    size_t numCell = T.numCell();
    for (size_t i=0; i<numCell; ++i)
      T.cell(i).calculatePCAPlane(vertexData);
  }
}

InitiateWallLength::
InitiateWallLength(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue ) 
{
  //
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "InitiateWallLength::InitiateWallLength() "
	      << "Uses one parameter, factor. " << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "InitiateWallLength::"
	      << "InitiateWallLength() "
	      << "No variable indices used." << std::endl;
    exit(0);
  }
  //
  //Set the variable values
  //
  setId("InitiateWallLength");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "factor";
  setParameterId( tmp );
}

void InitiateWallLength::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
	
{
  size_t dimension = vertexData[0].size();
  size_t numWall = T.numWall();
  
  for (size_t i=0; i<numWall; ++i) {
    double distance=0.0;
    size_t v1I=T.wall(i).vertex1()->index();
    size_t v2I=T.wall(i).vertex2()->index();
    for (size_t d=0; d<dimension; ++d )
      distance += (vertexData[v2I][d]-vertexData[v1I][d])*(vertexData[v2I][d]-vertexData[v1I][d]);
    distance = std::sqrt(distance);
    wallData[i][0] = parameter(0)*distance;
  }
}

void InitiateWallLength::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
}

InitiateWallVariableConstant::
InitiateWallVariableConstant(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue ) 
{
  //
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "InitiateWallVariableConstant::InitiateWallVariableConstant() "
	      << "Uses one parameter, the value. " << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "InitiateWallVariableConstant::"
	      << "InitiateWallVariableConstant() "
	      << "One variable index at first level used, the variable to be initiated"
	      << std::endl;
    exit(0);
  }
  //
  //Set the variable values
  //
  setId("InitiateWallVariableConstant");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "initValue";
  setParameterId( tmp );
}

void InitiateWallVariableConstant::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
	
{
  size_t dimension = vertexData[0].size();
  size_t numWall = T.numWall();
  
  for (size_t i=0; i<numWall; ++i) {
    wallData[i][variableIndex(0,0)] = parameter(0);
  }
}

void InitiateWallVariableConstant::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
}

InitiateWallMesh::
InitiateWallMesh(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > 
		 &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "InitiateWallMesh::InitiateWallMesh() "
	      << "Uses one parameter, #of added vertices per wall." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "InitiateWallMesh::"
	      << "InitiateWallMesh() "
	      << "No variable indices used." << std::endl;
    exit(0);
  }
  //
  //Set the variable values
  //
  setId("InitiateWallMesh");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "numAddedVerticesPerWall";
  setParameterId( tmp );
}

void InitiateWallMesh::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  size_t numWallOld = T.numWall();
  size_t numWall = numWallOld;
  size_t numVertex = T.numVertex();
  size_t numMesh = (size_t) parameter(0);
  if (!numMesh) 
    return;
  double scaleFactor = 1.0/(1.0+numMesh); 
  
  for (size_t i=0; i<numWallOld; ++i) {		
    // Add the new walls and vertices, and set indices, 
    // and wall lengths and vertex positions (into the data matrices)
    wallData[i][0] *= scaleFactor;
    size_t v1I = T.wall(i).vertex1()->index();
    size_t v2I = T.wall(i).vertex2()->index();
    std::vector<double> n(dimension);
    for (size_t d=0; d<dimension; ++d)
      n[d] = scaleFactor*(vertexData[v2I][d]-vertexData[v1I][d]);
    for (size_t k=0; k<numMesh; ++k) {
      T.addWall(T.wall(i));
      T.wall(numWall+k).setIndex(numWall+k);
      wallData.push_back(wallData[i]);//with the new correct length
      T.addVertex(T.vertex(v1I));
      T.vertex(numVertex+k).setIndex(numVertex+k);
      vertexData.push_back(vertexData[v1I]);
      for (size_t d=0; d<dimension; ++d)
	vertexData[numVertex+k][d] += (k+1)*n[d];
    }
    // Adjust cell connections
    size_t c1I = T.wall(i).cell1()->index();
    size_t c2I = T.wall(i).cell2()->index();
    for (size_t k=0; k<numMesh; ++k) {
      T.cell(c1I).addWall(T.wallP(numWall+k));
      T.cell(c2I).addWall(T.wallP(numWall+k));
      T.cell(c1I).addVertex(T.vertexP(numVertex+k));
      T.cell(c2I).addVertex(T.vertexP(numVertex+k));
    }
    // Adjust wall connections
    T.wall(i).setVertex2(T.vertexP(numVertex));
    for (size_t k=0; k<numMesh-1; ++k) {
      T.wall(numWall+k).setVertex1(T.vertexP(numVertex+k));
      T.wall(numWall+k).setVertex2(T.vertexP(numVertex+k+1));
    }
    T.wall(numWall+numMesh-1).setVertex1(T.vertexP(numVertex+numMesh-1));
    // Adjust vertex connections
    std::vector<Cell*> tmpCell;
    if (T.wall(i).cell1()!=T.background())
      tmpCell.push_back(T.wall(i).cell1());
    if (T.wall(i).cell2()!=T.background())
      tmpCell.push_back(T.wall(i).cell2());
    std::vector<Wall*> tmpWall(2);
    tmpWall[1] = T.wallP(i);
    for (size_t k=0; k<numMesh; ++k) {
      T.vertex(numVertex+k).setCell(tmpCell);
      tmpWall[0] = tmpWall[1];
      tmpWall[1] = T.wallP(numWall+k);
      T.vertex(numVertex+k).setWall(tmpWall);
    }
    size_t numVertexWall = T.vertex(v2I).numWall();
    size_t updatedFlag=0;
    for (size_t k=0; k<numVertexWall; ++k) {
      if (T.vertex(v2I).wall(k)==T.wallP(i)) {
	T.vertex(v2I).setWall(k,T.wallP(numWall+numMesh-1));
	++updatedFlag;
      }
    }
    assert(updatedFlag==1);
    numWall += numMesh;
    numVertex += numMesh;
  }
  // Sort the walls and vertices in the cells again since things have been added
  T.sortCellWallAndCellVertex();
  T.checkConnectivity(1);
}

void InitiateWallMesh::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
}

StrainTest::
StrainTest(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) 
{
  //
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()<2 ) {
    std::cerr << "StrainTest::StrainTest() "
	      << "Uses at least two parameters, version, angle and optional parameters." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "StrainTest::"
	      << "StrainTest() "
	      << "No variable indices used." << std::endl;
    exit(0);
  }
  //
  //Set the variable values
  //
  setId("StrainTest");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "version";
  tmp[1] = "angle";
  setParameterId( tmp );
}

void StrainTest::
initiate(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
}

void StrainTest::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  size_t dimension=vertexData[0].size();
  double angleRad = 2*3.14159*parameter(1)/360;
  assert(dimension==2);
  std::vector<double> n(dimension);
  n[0] = std::cos(angleRad);
  n[1] = std::sin(angleRad);
  if (parameter(0)==0.0) {
    if (numParameter()!=2) {
      std::cerr << "StrainTest::derivs() version " << parameter(0) << " requires " << 2
		<< " parameters." << std::endl;
      exit(-1);
    }
    for (size_t d=0; d<dimension; ++d) {
      double sign = vertexData[3][d]>0.0 ? 1.0 : -1.0;
      vertexDerivs[3][d] = n[d]*sign;
    }
  }
  else if (parameter(0)==1.0) {
    if (numParameter()!=2) {
      std::cerr << "StrainTest::derivs() version " << parameter(0) << " requires " << 2
		<< " parameters." << std::endl;
      exit(-1);
    }
    for (size_t vI=0; vI<vertexData.size(); ++vI) {
      for (size_t d=0; d<dimension; ++d) {
	double sign = vertexData[vI][d]>0.0 ? 1.0 : -1.0;
	vertexDerivs[vI][d] = n[d]*sign;
      }
    }
  }	
}

CalculateVertexStressDirection::CalculateVertexStressDirection(std::vector<double> &paraValue, 
							       std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "CalculateVertexStressDirection::CalculateVertexStressDirection() "
	      << "Uses one parameter: onlyInUpdateFlag\n";
    exit(0);
  }
  
  if (indValue.size() != 1) {
    std::cerr << "CalculateVertexStressDirection::CalculateVertexStressDirection() "
	      << "Wall force indexes.\n";
    exit(0);
  }
	
  setId("CalculateVertexStressDirection");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "onlyInUpdateFlag";
  setParameterId(tmp);
  
  for (size_t i = 0; i < indValue[0].size(); ++i) {
    wallForceIndexes_.push_back(indValue[0][i]);
  }
}

void CalculateVertexStressDirection::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs) 
{
  if (parameter(0) == 1.0) {
    size_t numVertex = T.numVertex();
    
    for (size_t i = 0; i < numVertex; ++i) {
      T.vertex(i).calculateStressDirection(vertexData, wallData, wallForceIndexes_);
    }
  }
}

void CalculateVertexStressDirection::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs) 
{
  if (parameter(0) != 1.0) {
    size_t numVertex = T.numVertex();
    for (size_t i=0; i<numVertex; ++i) {
      T.vertex(i).calculateStressDirection(vertexData, wallData, wallForceIndexes_);
    }
  }
}

void CalculateVertexStressDirection::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h)
{
  if (parameter(0) == 1.0) {
    size_t numVertex = T.numVertex();
    for (size_t i=0; i<numVertex; ++i) {
      T.vertex(i).calculateStressDirection(vertexData, wallData, wallForceIndexes_);
    }
  }
}

MoveVerticesRandomlyCapCylinder::MoveVerticesRandomlyCapCylinder(std::vector<double> &paraValue, 
							       std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "MoveVerticesRandomlyCapCylinder::MoveVerticesRandomlyCapCylinder() "
	      << "Uses one parameter: onlyInUpdateFlag\n";
    exit(0);
  }
  
  if (indValue.size() != 0) {
    std::cerr << "MoveVerticesRandomlyCapCylinder::MoveVerticesRandomlyCapCylinder() "
	      << "no parameter index.\n";
    exit(0);
  }
	
  setId("MoveVerticesRandomlyCapCylinder");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "onlyInUpdateFlag";
  setParameterId(tmp);
  
  
}

void MoveVerticesRandomlyCapCylinder::initiate(Tissue &T,
					       DataMatrix &cellData,
					       DataMatrix &wallData,
					       DataMatrix &vertexData,
					       DataMatrix &cellDerivs,
					       DataMatrix &wallDerivs,
					       DataMatrix &vertexDerivs)
{


  size_t numVertices = T.numVertex();
  double fac=parameter(0);
  
  double PI=3.14159265;
  double R=10;  
  double zmax=15;
  double zmin=-15;
  
  // Move vertices
  for( size_t VertexIndex=0 ; VertexIndex<numVertices ; ++VertexIndex ) {
    double a= rand();
    double b= rand();
    double c= RAND_MAX ;
    double d=a/c;
    double f=b/c;
    double x=vertexData[ VertexIndex][0];
    double y=vertexData[ VertexIndex][1];
    double z=vertexData[ VertexIndex][2];
   
    // if (z>zmin && z<zmax){
    //   vertexData[ VertexIndex][2]+=fac*(d-0.5);
    // }
   


    if ((z<zmin && z>-24.8) || ( z>zmax && z<24.8)){
      double teta=0;
      if (z<zmin){
    	teta=std::atan((std::sqrt(x*x+y*y))/(z-zmin));  
      }   
      if (z>zmax){
    	teta=std::atan((std::sqrt(x*x+y*y))/(z-zmax));  
      }   
      double phi=std::atan(y/x);     

      if (x<0 ){phi +=PI;}
      if (z<0 ){teta +=PI;}
 
      teta +=fac*(d-0.5);
      phi  +=0.3*fac*(f-0.5);
      if ( z<zmin && z>-24.8 ){
    	vertexData[ VertexIndex][0]=R*(std::sin(teta))*(std::cos(phi));
    	vertexData[ VertexIndex][1]=R*(std::sin(teta))*(std::sin(phi));
    	vertexData[ VertexIndex][2]=-15+R*(std::cos(teta));
      }
      if ( z>zmax && z< 24.8 ){
    	vertexData[ VertexIndex][0]=R*(std::sin(teta))*(std::cos(phi));
    	vertexData[ VertexIndex][1]=R*(std::sin(teta))*(std::sin(phi));
    	vertexData[ VertexIndex][2]=15+(R*(std::cos(teta)));
      }
    } 
   
    // if (z>zmin && z<zmax){
       
    //   double phi=std::atan(y/x);     
    //   if (x<0 ){phi +=PI;}
    //   phi  +=fac*(f-0.5);
    //   vertexData[ VertexIndex][0]=R*(std::cos(phi));
    //   vertexData[ VertexIndex][1]=R*(std::sin(phi));
      

    // }
    
    std::cerr<<"                                "<<a/c<<"  "<<b/c << std::endl;
    
  }
  
}






void MoveVerticesRandomlyCapCylinder::derivs(Tissue &T,
					    DataMatrix &cellData,
					    DataMatrix &wallData,
					    DataMatrix &vertexData,
					    DataMatrix &cellDerivs,
					    DataMatrix &wallDerivs,
					    DataMatrix &vertexDerivs) 
{
   
  
}



void MoveVerticesRandomlyCapCylinder::update(Tissue &T,
					    DataMatrix &cellData,
					    DataMatrix &wallData,
					    DataMatrix &vertexData,
					    double h)
{
  
}




scaleTemplate::scaleTemplate(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "scaleTemplate::scaleTemplate() "
	      << "Uses one parameter: scaling factor\n";
    exit(0);
  }
  
  if (indValue.size() != 0) {
    std::cerr << "scaleTemplate::scaleTemplate() "
	      << "no parameter index.\n";
    exit(0);
  }
	
  setId("scaleTemplate");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "onlyInUpdateFlag";
  setParameterId(tmp);
  
  
}

void scaleTemplate::initiate(Tissue &T,
					       DataMatrix &cellData,
					       DataMatrix &wallData,
					       DataMatrix &vertexData,
					       DataMatrix &cellDerivs,
					       DataMatrix &wallDerivs,
					       DataMatrix &vertexDerivs)
{

  size_t numVertices = T.numVertex();
  size_t numWall = T.numWall();
  double fac=parameter(0);
  
  // Move vertices
  for( size_t VertexIndex=0 ; VertexIndex<numVertices ; ++VertexIndex ) {
    // ad-hoc
    // if(vertexData[ VertexIndex][1]>=0)
    //   vertexData[ VertexIndex][1]=std::pow(vertexData[ VertexIndex][1],0.75);
       
    // if(vertexData[ VertexIndex][1]<0)
    //   vertexData[ VertexIndex][1]=-std::pow(-vertexData[ VertexIndex][1],0.75);
    vertexData[ VertexIndex][0]=fac*vertexData[ VertexIndex][0];
    vertexData[ VertexIndex][1]=fac*vertexData[ VertexIndex][1];
    vertexData[ VertexIndex][2]=fac*vertexData[ VertexIndex][2];     
  } 

  for (size_t i=0; i<numWall; ++i) {
    wallData[i][0] = fac*wallData[i][0];
  }
  // for (size_t i=0; i<numWall; ++i) {
  //   double distance=0.0;
  //   size_t v1I=T.wall(i).vertex1()->index();
  //   size_t v2I=T.wall(i).vertex2()->index();
  //   for (size_t d=0; d<dimension; ++d )
  //     distance += (vertexData[v2I][d]-vertexData[v1I][d])*(vertexData[v2I][d]-vertexData[v1I][d]);
  //   distance = std::sqrt(distance);
  //   wallData[i][0] = parameter(0)*distance;
  // } 
}

void scaleTemplate::derivs(Tissue &T,
			   DataMatrix &cellData,
			   DataMatrix &wallData,
			   DataMatrix &vertexData,
			   DataMatrix &cellDerivs,
			   DataMatrix &wallDerivs,
			   DataMatrix &vertexDerivs) 
{
}

void scaleTemplate::update(Tissue &T,
			   DataMatrix &cellData,
			   DataMatrix &wallData,
			   DataMatrix &vertexData,
			   double h)
{  
}

copyCellVector::copyCellVector(std::vector<double> &paraValue, 
                               std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 0) {
    std::cerr << "copyCellVector::copyCellVector() "
	      << "Uses no parameter\n";
    exit(0);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 2) {
    std::cerr << "scaleTemplate::scaleTemplate() "
	      << "one index level with 2 indices.\n";
    exit(0);
  }
	
  setId("copyCellVector");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // std::vector<std::string> tmp(numParameter());
  // tmp[0] = "onlyInUpdateFlag";
  // setParameterId(tmp);  
}

void copyCellVector::initiate(Tissue &T,
			      DataMatrix &cellData,
			      DataMatrix &wallData,
			      DataMatrix &vertexData,
			      DataMatrix &cellDerivs,
			      DataMatrix &wallDerivs,
			      DataMatrix &vertexDerivs)
{
  size_t numCells = T.numCell();
  // copy vectors
  for( size_t cellInd=0 ; cellInd<numCells ; ++cellInd ) {
    cellData[cellInd][variableIndex(0,1)  ]=cellData[cellInd][variableIndex(0,0)  ];
    cellData[cellInd][variableIndex(0,1)+1]=cellData[cellInd][variableIndex(0,0)+1];
    cellData[cellInd][variableIndex(0,1)+2]=cellData[cellInd][variableIndex(0,0)+2];
    cellData[cellInd][variableIndex(0,1)+3]=cellData[cellInd][variableIndex(0,0)+3];
  } 
}

void copyCellVector::derivs(Tissue &T,
                            DataMatrix &cellData,
                            DataMatrix &wallData,
                            DataMatrix &vertexData,
                            DataMatrix &cellDerivs,
                            DataMatrix &wallDerivs,
                            DataMatrix &vertexDerivs) 
{
}

void copyCellVector::update(Tissue &T,
                            DataMatrix &cellData,
                            DataMatrix &wallData,
                            DataMatrix &vertexData,
                            double h)
{
}

limitZdis::limitZdis(std::vector<double> &paraValue, 
                     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 0) {
    std::cerr << "limitZdis::limitZdis() "
	      << "Uses no parameter\n";
    exit(0);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 2) {
    std::cerr << "limitZdis::limitZdis() "
	      << "one index level with 2 indices.\n";
    exit(0);
  }
	
  setId("limitZdis");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // std::vector<std::string> tmp(numParameter());
  // tmp[0] = "onlyInUpdateFlag";
  // setParameterId(tmp);  
}

void limitZdis::initiate(Tissue &T,
                         DataMatrix &cellData,
                         DataMatrix &wallData,
                         DataMatrix &vertexData,
                         DataMatrix &cellDerivs,
                         DataMatrix &wallDerivs,
                         DataMatrix &vertexDerivs)
{
 size_t numCells = T.numCell();
 for (size_t cellIndex= 0; cellIndex< numCells; ++cellIndex)
   if(cellData[cellIndex][38]==-1 && (cellData[cellIndex][37]==-4 ||cellData[cellIndex][37]==0) ) // for hypocotyl
     { size_t numWalls = T.cell(cellIndex).numWall();      
       for (size_t wallindex=0; wallindex<numWalls; ++wallindex) { 
         size_t vInd= T.cell(cellIndex).vertex(wallindex)->index();
         if(vertexData[vInd][2]>-50)
           topVertices.push_back(vInd);
         if(vertexData[vInd][2]<-50)
           bottomVertices.push_back(vInd);
       }
     }
}

void limitZdis::derivs(Tissue &T,
                       DataMatrix &cellData,
                       DataMatrix &wallData,
                       DataMatrix &vertexData,
                       DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs,
                       DataMatrix &vertexDerivs) 
{
  size_t numTop=topVertices.size();
  size_t numBottom=bottomVertices.size();
  
  double tmpZforce=0;
  for (size_t i= 0; i< numTop; ++i)
    tmpZforce+=vertexDerivs[topVertices[i]][2];
  tmpZforce/=numTop;
  for (size_t i= 0; i< numTop; ++i)
    vertexDerivs[topVertices[i]][2]=tmpZforce;
  
  tmpZforce=0;
  for (size_t i= 0; i< numBottom; ++i)
    tmpZforce+=vertexDerivs[bottomVertices[i]][2];
  tmpZforce/=numBottom;
  for (size_t i= 0; i< numBottom; ++i)
    vertexDerivs[bottomVertices[i]][2]=tmpZforce;
}

void limitZdis::update(Tissue &T,
                       DataMatrix &cellData,
                       DataMatrix &wallData,
                       DataMatrix &vertexData,
                       double h)
{
}

randomizeMT::randomizeMT(std::vector<double> &paraValue, 
                         std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 4) {
    std::cerr << "randomizeMT::randomizeMT"
	      << "Uses two parameters parameter(0) for MT(0/1) and parameter(1) "
              << "for concentration(0/1) randomization in the range of [parameter(2), parameter(3)]"
              <<std::endl;
    exit(0);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 2) {
    std::cerr << "randomizeMT::randomizeMT()"
	      << "one index level with 2 indices for MT and concentration.\n";
    exit(0);
  }
	
  setId("randomizeMT");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // std::vector<std::string> tmp(numParameter());
  // tmp[0] = "onlyInUpdateFlag";
  // setParameterId(tmp);  
}

void randomizeMT::initiate(Tissue &T,
                           DataMatrix &cellData,
                           DataMatrix &wallData,
                           DataMatrix &vertexData,
                           DataMatrix &cellDerivs,
                           DataMatrix &wallDerivs,
                           DataMatrix &vertexDerivs)
{
  size_t numCells = T.numCell();
  double Min, Max;
  Min=parameter(2);
  Max=parameter(3);
  size_t MTInd=variableIndex(0,0);
  size_t conInd=variableIndex(0,1);

  //std::srand (time(NULL));

  for( size_t cellInd=0 ; cellInd<numCells ; ++cellInd ) {
    double ttmp=Min+(Max-Min)*((double) rand() / (RAND_MAX));
    if (parameter(1)==1)
      cellData[cellInd][conInd]=ttmp;
    // calculate the average normal vector to the cell plane
    size_t numVer=T.cell(cellInd).numVertex();
    std::vector<std::vector<double> > verticesPosition;
    verticesPosition.resize(numVer);
  
    // storing the vertex positions
    for(size_t k=0; k< numVer; ++k) {
      verticesPosition[k].resize(3);
      size_t Vind=T.cell(cellInd).vertex(k) -> index();
      verticesPosition[k][0]=vertexData[Vind][0];
      verticesPosition[k][1]=vertexData[Vind][1];
      verticesPosition[k][2]=vertexData[Vind][2];
    }    
    // Finding the average normal vector to the cell plane
    std::vector<double> normal(3);
    normal[0]=0;
    normal[1]=0;
    normal[2]=0;

    for (size_t k=1; k< numVer-1; ++k) {
      std::vector<double> x0(3),x1(3),x2(3);
      size_t ind1=T.cell(cellInd).vertex(0)   -> index();
      size_t ind2=T.cell(cellInd).vertex(k)   -> index();
      size_t ind3=T.cell(cellInd).vertex(k+1) -> index();
      
      //normal to the element
      std::vector<double> temp(3);
      temp[0]=(vertexData[ind2][1]-vertexData[ind1][1])*(vertexData[ind3][2]-vertexData[ind1][2])
        -(vertexData[ind2][2]-vertexData[ind1][2])*(vertexData[ind3][1]-vertexData[ind1][1]);
      temp[1]=(vertexData[ind2][2]-vertexData[ind1][2])*(vertexData[ind3][0]-vertexData[ind1][0])
        -(vertexData[ind2][0]-vertexData[ind1][0])*(vertexData[ind3][2]-vertexData[ind1][2]);
      temp[2]=(vertexData[ind2][0]-vertexData[ind1][0])*(vertexData[ind3][1]-vertexData[ind1][1])
        -(vertexData[ind2][1]-vertexData[ind1][1])*(vertexData[ind3][0]-vertexData[ind1][0]);
    
      // it is area weighted because of the size of outer product
      normal[0]+=temp[0];
      normal[1]+=temp[1];
      normal[2]+=temp[2];
    }
    double norm=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    // normalizing the normal vector
    normal[0]/=norm;
    normal[1]/=norm;
    normal[2]/=norm;
  
    // choose a rondom direction 
    //double teta=(parameter(0)*3.1415/180)*(1-2*((double) rand() / (RAND_MAX)));
    double randVec[3]={1-2*((double) rand() / (RAND_MAX)),
                       1-2*((double) rand() / (RAND_MAX)),
                       1-2*((double) rand() / (RAND_MAX))};
    
    double randMT[3];
    //outer product between random direction and cell normal
    randMT[0]=randVec[1]*normal[2]-randVec[2]*normal[1];
    randMT[1]=randVec[2]*normal[0]-randVec[0]*normal[2];
    randMT[2]=randVec[0]*normal[1]-randVec[1]*normal[0];

    double tmp=std::sqrt(randMT[0]*randMT[0]+randMT[1]*randMT[1]+randMT[2]*randMT[2]);
    randMT[0]/=tmp;
    randMT[1]/=tmp;
    randMT[2]/=tmp;
    //insert the random direction to the MT index of the cell vector 
    if(parameter(0)==1){
      cellData[cellInd][variableIndex(0,0)  ]=randMT[0];
      cellData[cellInd][variableIndex(0,0)+1]=randMT[1];
      cellData[cellInd][variableIndex(0,0)+2]=randMT[2];
      cellData[cellInd][variableIndex(0,0)+3]=1;
    }  
  } 
}

void randomizeMT::derivs(Tissue &T,
                         DataMatrix &cellData,
                         DataMatrix &wallData,
                         DataMatrix &vertexData,
                         DataMatrix &cellDerivs,
                         DataMatrix &wallDerivs,
                         DataMatrix &vertexDerivs) 
{
}

void randomizeMT::update(Tissue &T,
                         DataMatrix &cellData,
                         DataMatrix &wallData,
                         DataMatrix &vertexData,
                         double h)
{
}

restrictVertexRadially::restrictVertexRadially(std::vector<double> &paraValue, 
                                               std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 0) {
    std::cerr << "restrictVertexRadially::restrictVertexRadially() "
	      << "Uses no parameter\n";
    exit(0);
  }
  
  if (indValue.size() != 0) {
    std::cerr << "restrictVertexRadially::restrictVertexRadially() "
	      << "no index.\n";
    exit(0);
  }
  setId("restrictVertexRadially");
  setParameter(paraValue);  
  setVariableIndex(indValue);  
}

void restrictVertexRadially::derivs(Tissue &T,
                                    DataMatrix &cellData,
                                    DataMatrix &wallData,
                                    DataMatrix &vertexData,
                                    DataMatrix &cellDerivs,
                                    DataMatrix &wallDerivs,
                                    DataMatrix &vertexDerivs) 
{
  size_t numVertices = T.numVertex();
    
  // restrict vertices radially in xy plane
  for( size_t VertexIndex=0 ; VertexIndex<numVertices ; ++VertexIndex ) {
    double tmp=std::sqrt(vertexData[VertexIndex][0]*vertexData[VertexIndex][0]+
                         vertexData[VertexIndex][1]*vertexData[VertexIndex][1]);
    double tmp1=vertexData[VertexIndex][0]/tmp;
    double tmp2=vertexData[VertexIndex][1]/tmp;
    tmp=(vertexDerivs[VertexIndex][0]*tmp1+vertexDerivs[VertexIndex][1]*tmp2);
    vertexDerivs[VertexIndex][0]=tmp*tmp1;
    vertexDerivs[VertexIndex][1]=tmp*tmp2;
  }
}

void restrictVertexRadially::update(Tissue &T,
                                    DataMatrix &cellData,
                                    DataMatrix &wallData,
                                    DataMatrix &vertexData,
                                    double h)
{
}

CreationPrimordiaTime::
CreationPrimordiaTime(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
{  

  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=5 ) {
    std::cerr << "CreationPrimordiaTime::CreationPrimordiaTime() "
	      << "Uses five parameter(s) " 
              << "k_c, the constant production rate," 
              << "tp, time intervall for creating primordia,"
              << "teta, angle between primordia, "
              << "z, vertical distance from max for new primordium."
              << "and a flag(equal to 1) if constant -amount- creation is needed"
              <<std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 2  ) {
    std::cerr << "CreationPrimordiaTime::"
	      << "CreationPrimordiaTime() "
	      << "one levels of indices used: index for variable to be updated given "
              << "as first indexand COM index in the second level" 
              << std::endl;
    exit(0);
  }
  	
  // Set the variable values
  setId("CreationPrimordiaTime");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_c";
  tmp[1] = "delta_time";
  tmp[2] = "delta_teta";
  tmp[3] = "z_distanceTip";
  tmp[4] = "production_flag";
 
  setParameterId(tmp);
}

void CreationPrimordiaTime::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //Do the update for each cell
  size_t numCells = T.numCell();
  size_t cIndex = variableIndex(0,0);
  
  if(parameter(4)==1)
    //For the cells in the list
    for (size_t cellI = 0; cellI < proCells.size(); ++cellI) 
      cellDerivs[proCells[cellI]][cIndex] += parameter(0)/T.cell(cellI).calculateVolume(vertexData);
  
  if(parameter(4)==0)
    for (size_t cellI = 0; cellI < proCells.size(); ++cellI) 
      cellDerivs[proCells[cellI]][cIndex] += parameter(0);  
}

void CreationPrimordiaTime::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h)
{

  size_t nCells=T.numCell();
  size_t comInd=variableIndex(0,1);
  static double timeT=0;
  static size_t n=1;
  double deltaz=20 ;
  double R=40;
 
  timeT+=h;
  double tmp=0;       
  if(timeT>n*parameter(1)) {
    //make a primordium

    // determining the tip as maximum z
    double maxZ=cellData[0][comInd+2];
    for(size_t i=1; i< nCells; ++i)
      if(cellData[i][comInd+2]>maxZ)
        maxZ=cellData[i][comInd+2]; 
    
    size_t pInd=0;
    double distance=100000;
    
    for(size_t i=1; i< nCells; ++i)
      if(cellData[i][comInd+2]<maxZ-parameter(3) &&  cellData[i][comInd+2]>maxZ-parameter(3)-deltaz ){
        tmp=3.1415*fmod(n*parameter(2),360)/180;       
        double tmpDistance=(cellData[i][comInd+2]-maxZ+parameter(3))*(cellData[i][comInd+2]-maxZ+parameter(3))
          +(cellData[i][comInd+1]-R*std::sin(tmp))*(cellData[i][comInd+1]-R*std::sin(tmp))
          +(cellData[i][comInd]-R*std::cos(tmp))*(cellData[i][comInd]-R*std::cos(tmp));
        if(tmpDistance<distance){
          pInd=i;
          distance=tmpDistance;
        }       
      }
   
    if (pInd!=0)   
      proCells.push_back(pInd);
    
    n++;
    // std::cerr<<tmp*180/3.1415<<" "<<std::cos(tmp)<<" "<<std::sin(tmp)<<std::endl;
    // std::cerr<<"end primordium ";
    // for(size_t i=0; i< proCells.size(); ++i)
    //   std::cerr<<proCells[i]<<"  "; 
    // std::cerr<<std::endl; 
  }
}

VertexFromRotationalForceLinear::
VertexFromRotationalForceLinear(std::vector<double> &paraValue, 
                                std::vector< std::vector<size_t> > 
                                &indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()<2 || paraValue.size()>4 ) {
    std::cerr << "VertexFromRotationalForceLinear::"
	      << "VertexFromRotationalForceLinear() "
	      << "Uses a force vector that should be in one (x), two (x,y) or three (x,y,z) "
	      << "dimensions plus a deltaT that sets the time the linear increase is applied." 
	      << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() < 1 ) {
    std::cerr << "VertexFromRotationalForceLinear::"
	      << "VertexFromRotationalForceLinear() "
	      << "List of vertex indices given in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexFromRotationalForceLinear");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "F_x";
  if (numParameter()==2)
    tmp[1] = "deltaT";
  else if (numParameter()==3) {
    tmp[1] = "F_y";
    tmp[2] = "deltaT";
  }
  else {
    tmp[1] = "F_y";
    tmp[2] = "F_z";
    tmp[3] = "deltaT";
  }
  timeFactor_=0.0;
  setParameterId( tmp );
}






void VertexFromRotationalForceLinear::
initiate(Tissue &T,
         DataMatrix &cellData,
         DataMatrix &wallData,
         DataMatrix &vertexData,
         DataMatrix &cellDerivs,
         DataMatrix &wallDerivs,
         DataMatrix &vertexDerivs ){


  // find the boundary vertices and store the positions and indices
  size_t numVertices =T.numVertex();
  double UpT=-50;
  double DnT=-250;
  for(size_t vIndex=0; vIndex<numVertices; ++vIndex){
    if(vertexData[vIndex][2]>UpT){
      size_t n=boundVerticesUp.size()+1;
      boundVerticesUp.resize(n);
      boundVerticesUp[n-1].push_back(vertexData[vIndex][0]);
      boundVerticesUp[n-1].push_back(vertexData[vIndex][1]);
      boundVerticesUp[n-1].push_back(vertexData[vIndex][2]);
      boundVerticesUp[n-1].push_back(vIndex);
    }
    if(vertexData[vIndex][2]<DnT){
      size_t n=boundVerticesDn.size()+1;
      boundVerticesDn.resize(n);
      boundVerticesDn[n-1].push_back(vertexData[vIndex][0]);
      boundVerticesDn[n-1].push_back(vertexData[vIndex][1]);
      boundVerticesDn[n-1].push_back(vertexData[vIndex][2]);
      boundVerticesDn[n-1].push_back(vIndex);
    }
  }
}

void VertexFromRotationalForceLinear::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  // set derivetives of boundary to zero
  size_t UpN, DnN;
  UpN=boundVerticesUp.size();
  DnN=boundVerticesDn.size();
  for(size_t v=0; v<UpN; ++v)
    for(size_t d=0; d<3; ++d)
      vertexDerivs[boundVerticesUp[v][3]][d]=0;
  for(size_t v=0; v<DnN; ++v)
    for(size_t d=0; d<3; ++d)
      vertexDerivs[boundVerticesDn[v][3]][d]=0;
  
}

void VertexFromRotationalForceLinear::update(Tissue &T,
				   DataMatrix &cellData,
				   DataMatrix &wallData,
				   DataMatrix &vertexData,
				   double h)
{
  // rotate the position plane of boundaries
  size_t UpN, DnN;
  UpN=boundVerticesUp.size();
  DnN=boundVerticesDn.size();
  double deltaTet,tet;
  deltaTet=0.2;
  tet=h*deltaTet;  
  double centerUp[3]={0,0,-107.50};
  double centerDn[3]={0,0,-207.50};
  double rot[3][3]= { { 0 , 0 , 0 } ,
                      { 0 , 1 , 0 } ,
                      { 0 , 0 , 0 } };

  rot[0][0]= std::cos(tet);
  rot[0][2]= std::sin(tet);
  rot[2][0]=-std::sin(tet);
  rot[2][2]= std::cos(tet);

  // // rotate

  for(size_t v=0; v<UpN; ++v)
    for(size_t d=0; d<3; ++d)
      vertexData[boundVerticesUp[v][3]][d]= 
        rot[d][0]*(boundVerticesUp[v][0]-centerUp[0])+
        rot[d][1]*(boundVerticesUp[v][1]-centerUp[1])+
        rot[d][2]*(boundVerticesUp[v][2]-centerUp[2])+
        +centerUp[d];
  
  rot[0][2]*=-1;
  rot[2][0]*=-1;
  
  for(size_t v=0; v<DnN; ++v)
    for(size_t d=0; d<3; ++d)
      vertexData[boundVerticesDn[v][3]][d]= 
        rot[d][0]*(boundVerticesDn[v][0]-centerDn[0])+
        rot[d][1]*(boundVerticesDn[v][1]-centerDn[1])+
        rot[d][2]*(boundVerticesDn[v][2]-centerDn[2])+
        +centerDn[d];
  
  for(size_t v=0; v<UpN; ++v)
    for(size_t d=0; d<3; ++d)
      boundVerticesUp[v][d]=vertexData[boundVerticesUp[v][3]][d];
  for(size_t v=0; v<DnN; ++v)
    for(size_t d=0; d<3; ++d)
      boundVerticesDn[v][d]=vertexData[boundVerticesDn[v][3]][d];
  


}




ThresholdSwitch::ThresholdSwitch(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "ThresholdSwitch::ThresholdSwitch()  parameter used "
        << "Threshold" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=1 ) {
    std::cerr << "ThresholdSwitch::ThresholdSwitch() "
              << "Two levels of variable indices are used, "
              << "one for the threshold variable (single index)"
              << ", and one for a list of variables to change" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "const";
  tmp[1] = "switchtype";
  setParameterId( tmp );
}

void ThresholdSwitch::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void ThresholdSwitch::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void ThresholdSwitch::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex = variableIndex(0,0);
  size_t c2Index = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex]>=parameter(0)  ) {
          cellData[cellI][c2Index]=1;
      }
    else if ( cellData[cellI][cIndex]<parameter(0) && parameter(1)==0 )
  {cellData[cellI][c2Index]=0;}

}
}




AndGate::AndGate(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "AndGate::AndGate()  parameter used "
        << "gatetype" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=1 ) {
    std::cerr << "AndGate::AndGate() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are two indices "
              << ", and one for the output variables" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "gatetype";
  setParameterId( tmp );
}

void AndGate::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void AndGate::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void AndGate::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    size_t gate_condition=0;

    if (cellData[cellI][cIndex_input1]==1 && cellData[cellI][cIndex_input2]==1)
        {gate_condition=1;}

    if (gate_condition==1) {
          cellData[cellI][cIndex_output]=1;
      }
    else if (gate_condition==0 && parameter(0)==0)
     {cellData[cellI][cIndex_output]=0;}

}
}



AndNotGate::AndNotGate(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "AndNotGate::AndNotGate()  parameter used "
        << "gatetype" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=1 ) {
    std::cerr << "AndNotGate::AndNotGate() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are two indices "
              << ", and one for the output variables" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "gatetype";
  setParameterId( tmp );
}

void AndNotGate::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void AndNotGate::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void AndNotGate::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    size_t gate_condition=0;

    if (cellData[cellI][cIndex_input1]==1 && cellData[cellI][cIndex_input2]==0)
        {gate_condition=1;}

    if (gate_condition==1) {
          cellData[cellI][cIndex_output]=1;
      }
    else if (gate_condition==0 && parameter(0)==0)
     {cellData[cellI][cIndex_output]=0;}

}
}



AndSpecialGate::AndSpecialGate(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "AndSpecialGate::AndSpecialGate()  does not use any parameters." << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=3 || indValue[1].size()!=1 ) {
    std::cerr << "AndSpecialGate::AndSpecialGate() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are three indices "
              << ", and one for the output variables" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue); 
  setVariableIndex(indValue);
  //
}

void AndSpecialGate::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void AndSpecialGate::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void AndSpecialGate::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_input3 = variableIndex(0,2);
  size_t cIndex_output = variableIndex(1,0);

  //size_t cond=0;
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  
    cellData[cellI][cIndex_output]=0;
   // if (cellData[cellI][cIndex_input2]==1 || cellData[cellI][cIndex_input2]==0){cond=1;}
   // if (cellData[cellI][cIndex_input2]==0){cond=1;}
   //    if (cellData[cellI][cIndex_input1]==1 && cond==1 && cellData[cellI][cIndex_input3]==0)

    if (cellData[cellI][cIndex_input1]==1 && cellData[cellI][cIndex_input2]==0 && cellData[cellI][cIndex_input3]==0)
        {cellData[cellI][cIndex_output]=1;}


  }
}




AndSpecialGate2::AndSpecialGate2(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "AndSpecialGate2::AndSpecialGate2()  does not use any parameters." << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=3 || indValue[1].size()!=1 ) {
    std::cerr << "AndSpecialGate2::AndSpecialGate2() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are three indices "
              << ", and one for the output variables" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue); 
  setVariableIndex(indValue);
  //
}


void AndSpecialGate2::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void AndSpecialGate2::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void AndSpecialGate2::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_input3 = variableIndex(0,2);
  size_t cIndex_output = variableIndex(1,0);

  //size_t cond=0;
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  


    if (cellData[cellI][cIndex_input1]==1 && cellData[cellI][cIndex_input2]==1 && cellData[cellI][cIndex_input3]==0)
        {cellData[cellI][cIndex_output]=1;}


  }
}




AndSpecialGate3::AndSpecialGate3(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "AndSpecialGate3::AndSpecialGate3()  needs one parameter." << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=3 || indValue[1].size()!=1 ) {
    std::cerr << "AndSpecialGate3::AndSpecialGate3() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are three indices "
              << ", and one for the output variables" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue); 
  setVariableIndex(indValue);
  //
}


void AndSpecialGate3::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void AndSpecialGate3::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void AndSpecialGate3::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_input3 = variableIndex(0,2);
  size_t cIndex_output = variableIndex(1,0);

  //size_t cond=0;
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

   // if (cellData[cellI][cIndex_input2]==1 || cellData[cellI][cIndex_input2]==0){cond=1;}
   // if (cellData[cellI][cIndex_input2]==0){cond=1;}
   //    if (cellData[cellI][cIndex_input1]==1 && cond==1 && cellData[cellI][cIndex_input3]==0)

    if (cellData[cellI][cIndex_input1]>parameter(0) && cellData[cellI][cIndex_input2]==1 && cellData[cellI][cIndex_input3]==0)
        {cellData[cellI][cIndex_output]=1;}
   //   else
  //      {cellData[cellI][cIndex_output]=0;}

    //if (cellData[cellI][cIndex_input1]==1 && cellData[cellI][cIndex_input2]==1)
        //{if(cellData[cellI][cIndex_input3]==0)
         // {cellData[cellI][cIndex_output]=1;}
        //else
         // {cellData[cellI][cIndex_output]=0;}
        //}
       // else
      //  {cellData[cellI][cIndex_output]=0;}

  }
}

AndGateCount::AndGateCount(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "No parameters are used in AndGateCount" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=1 ) {
    std::cerr << "AndGateCount::AndGateCount() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are two indices "
              << ", and one for the output variable" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}

void AndGateCount::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void AndGateCount::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void AndGateCount::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex_input1]==1 && cellData[cellI][cIndex_input2]==1)
        {cellData[cellI][cIndex_output]+=1;}

}
}



OrGateCount::OrGateCount(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "No parameters are used in OrGateCount" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=1 ) {
    std::cerr << "OrGateCount::OrGateCount() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are two indices "
              << ", and one for the output variable" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}

void OrGateCount::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void OrGateCount::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void OrGateCount::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex_input1]==1 || cellData[cellI][cIndex_input2]==1)
        {cellData[cellI][cIndex_output]+=1;}

}
}


OrSpecialGateCount::OrSpecialGateCount(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "No parameters are used in OrSpecialGateCount" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=1 ) {
    std::cerr << "OrSpecialGateCount::OrSpecialGateCount() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are two indices "
              << ", and one for the output variable" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}

void OrSpecialGateCount::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void OrSpecialGateCount::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void OrSpecialGateCount::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex_input1]==1 || cellData[cellI][cIndex_input2]>0)
        {cellData[cellI][cIndex_output]+=1;}

}
}





AndThresholdsGate::AndThresholdsGate(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "AndThresholdsGate::AndThresholdsGate() 2 parameter used "
        << "which are the two thresholds." << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=1 ) {
    std::cerr << "AndThresholdsGate::AndThresholdsGate() "
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are two indices "
              << ", and one for the output variables" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "gatetype";
  setParameterId( tmp );
}

void AndThresholdsGate::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void AndThresholdsGate::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void AndThresholdsGate::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex_input1]>parameter(0) && cellData[cellI][cIndex_input2]>parameter(1))
        {cellData[cellI][cIndex_output]=1;}

}
}

Count::Count(std::vector<double> &paraValue,
  std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "No parameters are used in Count" << std::endl;
    exit(0);
  }
  if( indValue.size()!=1 || indValue[0].size()!=1 ) {
    std::cerr << "Count::Count() "
              << "One level of variable index is used "
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
}


void Count::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void Count::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void Count::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_output = variableIndex(0,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

   cellData[cellI][cIndex_output]+=1.0;
          

}
}





FlagCount::FlagCount(std::vector<double> &paraValue,
  std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "No parameters are used in FlagCount" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=1 ) {
    std::cerr << "FlagCount::FlagCount() "
              << "Two levels of variable indices are used, "
              << "One for the input variable, "
              << ", and one for the output variable" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}


void FlagCount::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void FlagCount::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void FlagCount::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input = variableIndex(0,0);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex_input]==1)
        {cellData[cellI][cIndex_output]+=1.0;
          }

}
}


ThresholdReset::ThresholdReset(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "ThresholdReset::ThresholdReset()  parameter used "
        << "Threshold" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=1 ) {
    std::cerr << "ThresholdReset::ThresholdReset() "
              << "Two levels of variable indices are used, "
              << "one for the threshold variable (single index)"
              << ", and one for a list of variables to change" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "const";
  tmp[1] = "switchtype";
  setParameterId( tmp );
}

void ThresholdReset::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void ThresholdReset::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void ThresholdReset::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex = variableIndex(0,0);
  size_t c2Index = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex]>=parameter(0)) {
          cellData[cellI][c2Index]=0;
      }
    else if ( cellData[cellI][cIndex]<parameter(0) && parameter(1)==0 )
  {cellData[cellI][c2Index]=0;}

}
}




ThresholdNoisyReset::ThresholdNoisyReset(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "ThresholdNoisyReset::ThresholdNoisyReset()  parameter used "
        << "Threshold and noise amplitude" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=1 ) {
    std::cerr << "ThresholdNoisyReset::ThresholdNoisyReset() "
              << "Two levels of variable indices are used, "
              << "one for the threshold variable (single index)"
              << ", and one for a list of variables to change" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "const";
  tmp[1] = "amplitude";
  setParameterId( tmp );
}

void ThresholdNoisyReset::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void ThresholdNoisyReset::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void ThresholdNoisyReset::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex = variableIndex(0,0);
  size_t c2Index = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex]>=parameter(0)  ) {
          double rrr;
          //rrr= (myRandom::Rnd()-0.5);
          rrr= (myRandom::Rnd())*parameter(1);

          cellData[cellI][c2Index]=rrr;
      }

}
}




ThresholdResetAndCount::ThresholdResetAndCount(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indexes
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "ThresholdResetAndCount::ThresholdResetAndCount()  parameter used "
        << "Threshold" << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=2 ) {
    std::cerr << "ThresholdResetAndCount::ThresholdResetAndCount() "
              << "Two levels of variable indices are used, "
              << "one for the threshold variable (single index)"
              << ", and one for a list of variables to change" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "const";
  tmp[1] = "switchtype";
  setParameterId( tmp );
}

void ThresholdResetAndCount::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void ThresholdResetAndCount::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void ThresholdResetAndCount::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex = variableIndex(0,0);
  size_t c2Index = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex]>=parameter(0)  ) {
          double rrr;
          rrr= (myRandom::Rnd()-0.5);
          cellData[cellI][c2Index]=rrr;cellData[cellI][13]=cellData[cellI][13]+1;
      }


}
}



FlagNoisyReset::FlagNoisyReset(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "FlagNoisyReset::FlagNoisyReset() 2 parameters used "
        << "Threshold and  noise intensity. " << std::endl;
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=1 ) {
    std::cerr << "FlagNoisyReset::FlagNoisyReset() "
              << "Two levels of variable indices are used, "
              << "one for the threshold variable (single index)"
              << ", and one for a list of variables to change" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "const";
  tmp[1] = "switchtype";
  setParameterId( tmp );
}

void FlagNoisyReset::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void FlagNoisyReset::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void FlagNoisyReset::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();

    size_t cIndex = variableIndex(0,0);
    size_t c2Index = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex]==parameter(0)) {
          double rrr;
          rrr= (myRandom::Rnd()-0.5)*parameter(1);
          cellData[cellI][c2Index]=rrr;
      }

}
}



ThresholdAndFlagNoisyReset::ThresholdAndFlagNoisyReset(std::vector<double> &paraValue, 
                                                       std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=3 ) {
    std::cerr << "ThresholdAndFlagNoisyReset::ThresholdAndFlagNoisyReset()  requires 3 parameters"
        << "Threshold, flag value, and noise amplitude for resetting" << std::endl;
    exit(0);
  }


  if( indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=1 ) {
    std::cerr << "ThresholdAndFlagNoisyReset::ThresholdAndFlagNoisyReset()"
              << "Two levels of variable indices are used, "
              << "One for the input variables, which are two indices "
              << ", and one for the output variable" 
              << std::endl;
    exit(0);
  }

  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "const";
  tmp[1] = "switchtype";
  setParameterId( tmp );
}

void ThresholdAndFlagNoisyReset::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void ThresholdAndFlagNoisyReset::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void ThresholdAndFlagNoisyReset::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  

  size_t cIndex_input1 = variableIndex(0,0);
  size_t cIndex_input2 = variableIndex(0,1);
  size_t cIndex_output = variableIndex(1,0);

  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if ( cellData[cellI][cIndex_input1]>parameter(0) && cellData[cellI][cIndex_input2]==parameter(1)) {
          double rrr;
          rrr= (myRandom::Rnd()-0.5)*parameter(2);
          cellData[cellI][cIndex_output]=rrr;
      }

}
}



FlagAddValue::FlagAddValue(std::vector<double> &paraValue,
  std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "FlagAddValue::FlagAddValue  "
        << "Uses one parameter: add_value\n";
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=1 ) {
    std::cerr << "FlagAddValue::FlagAddValue() "
              << "Two levels of variable indices are used, "
              << "One for the input variable, "
              << ", and one for the output variable" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}


void FlagAddValue::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void FlagAddValue::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void FlagAddValue::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input = variableIndex(0,0);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  

    if (cellData[cellI][cIndex_input]==1)
        {cellData[cellI][cIndex_output]+=parameter(0);
          }

}
}



CopyVariable::CopyVariable(std::vector<double> &paraValue,
  std::vector< std::vector<size_t> > &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "CopyVariable::CopyVariable  "
        << "has no parameters. \n";
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=1 ) {
    std::cerr << "CopyVariable::CopyVariable() "
              << "Two levels of variable indices are used, "
              << "One for the input variable, "
              << ", and one for the output variable" 
              << std::endl;
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("add");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}


void CopyVariable::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  // Nothing to be done for the derivative function.
}

  
void CopyVariable::
derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex)
{
  // Nothing to be done for the derivative function.
}


void CopyVariable::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       double h) 
{

    size_t numCells = T.numCell();
  
  size_t cIndex_input = variableIndex(0,0);
  size_t cIndex_output = variableIndex(1,0);


  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  
    cellData[cellI][cIndex_output]=cellData[cellI][cIndex_input];
}

}
