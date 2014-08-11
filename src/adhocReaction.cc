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



//--------------------------------------------


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
      size_t numVertexWalls=T.vertex(vertexIndex).numWall();        //take the number of walls connected to it       //std::cerr<<"boundary vertex "<< vertexIndex  <<" has "<< numVertexWalls << " walls"<< std::endl;
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
  for (size_t vertex=0; vertex<numBoundaryVertices-1; vertex++){//for the boundary vertices
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


