 

/////program to fill loops according to shortes-loop-paths combining graph_BiConnectedComponents_Boost remove_graph-junctions_01 shortest-loop-path_03 and  vtkContourTriangulator
//02: recursive version
//03: debugging

#include <vtkSmartPointer.h>


#include <vtkXMLPolyDataReader.h>//for vtp-files (cannot contain 3D cells?)

//#include <vtkCleanPolyData.h>
//#include <vtkTriangleFilter.h>
//#include <vtkExtractEdges.h>

#include <vtkPolyDataToGraph.h>
#include <vtkBoostBiconnectedComponents.h>

#include <vtkThresholdGraph.h>
// #include <vtkThreshold.h>
// #include <vtkThresholdPoints.h>
#include <vtkSelectionNode.h>
#include <vtkSelectionSource.h>
#include <vtkExtractSelectedGraph.h>
#include <vtkRemoveIsolatedVertices.h>
#include <vtkBoostConnectedComponents.h>
#include <vtkIntArray.h>

#include <vtkIdList.h>
#include <vtkDijkstraGraphGeodesicPath.h>

#include <vtkContourTriangulator.h>
//#include <.h>

#include <vtkGraphToPolyData.h>

#include <vtkXMLPolyDataWriter.h>//for vtp-files 

//for isNestedLoop(
#include <vtkVertexListIterator.h>
#include <vtkAdjacentVertexIterator.h>

//for fill_loops(
#include <vtkAppendPolyData.h>
#include <vtkDataSetAttributes.h>
#include <vtkIdTypeArray.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkInEdgeIterator.h>
//#include <.h>

//for ProgressFunction(
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

#define VERBOSE 0
#define P_VERBOSE 1
//#define COM_INT= -1 //none biconnected components


#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()



void ProgressFunction( vtkObject* caller, long unsigned int eventId, void* clientData, void* callData ){
    vtkAlgorithm *d= static_cast<vtkAlgorithm*>(caller);
    if(P_VERBOSE) fprintf(stderr, "\rFilter progress: %5.1f%%", 100.0 * d->GetProgress());
    if(P_VERBOSE) std::cerr.flush();
    }


bool isNestedLoop(vtkSmartPointer<vtkGraph> graph){

  if (graph->GetNumberOfVertices() <= 2){
    return 0;
  }
  if (graph->GetNumberOfEdges() <= 2){
    return 0;
  }

  VTK_CREATE(vtkVertexListIterator, it);
  graph->GetVertices(it);

  bool ans;
  while(ans= it->HasNext()){ 
    vtkIdType nextVertex= it->Next();

    VTK_CREATE(vtkAdjacentVertexIterator, itc); 
    graph->GetAdjacentVertices(nextVertex, itc);
	
    int num_of_adjacent_vertices= 0;
    while (itc->HasNext()){
      vtkIdType u = itc->Next();
      num_of_adjacent_vertices++;
    }

    if (num_of_adjacent_vertices > 2){
      return 1;
    }
  }
  return 0;
}

vtkIdType getEndpointId(vtkSmartPointer<vtkPolyData> mesh){

  VTK_CREATE(vtkIdList, cellIdList);
  for(vtkIdType i= 0; i < mesh->GetNumberOfPoints(); i++){
    mesh->GetPointCells(i, cellIdList); //Make sure that routine BuildLinks() has been called. 

    if(cellIdList->GetNumberOfIds() == 1){
      return i;
    }
  }
    
  return -1;
}


vtkSmartPointer<vtkPolyData> fill_loops(vtkSmartPointer<vtkPolyData> in_mesh, vtkSmartPointer<vtkPolyData> connected_filled_loops){


  int COM_INT= -1; //none biconnected components
  VTK_CREATE(vtkCallbackCommand, progressCallback);
  progressCallback->SetCallback(ProgressFunction);


  vtkSmartPointer<vtkPolyData> sub_mesh;

  VTK_CREATE(vtkAppendPolyData, append);
  append->AddInputData(connected_filled_loops);

  VTK_CREATE(vtkContourTriangulator, ct);
  ct->TriangulationErrorDisplayOn();

  /////////////////Identifying loops

  VTK_CREATE(vtkPolyDataToGraph, polyDataToGraphFilter);
  polyDataToGraphFilter->SetInputData(in_mesh);

  if(P_VERBOSE) std::cout << "Executing vtkPolyDataToGraph..." << std::endl;
  polyDataToGraphFilter->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  polyDataToGraphFilter->Update();
  if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;



  //vtkGraph* graph= polyDataToGraphFilter->GetOutput();


  VTK_CREATE(vtkBoostBiconnectedComponents, bbc);
  bbc->SetInputConnection(polyDataToGraphFilter->GetOutputPort());
  bbc->SetOutputArrayName("biconnected_component");

  if(P_VERBOSE) std::cout << "Executing vtkBoostBiconnectedComponents..." << std::endl;
  bbc->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  bbc->Update();
  if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;



  /////////////////Identifying loops done.


  /////////////////Removing none-loop vertices



  /////////////////Removing only vertices connected to purely the same kind done.

    
  VTK_CREATE(vtkSelectionSource, sel);
  //sel->SetContentType(vtkSelectionNode::THRESHOLDS); //Thresholds
  sel->SetContentType(vtkSelectionNode::INDICES);
  sel->SetFieldType(vtkSelectionNode::VERTEX); //Edge
  sel->SetArrayName("biconnected_component");
  //sel->AddThreshold(0,bbc->GetOutput()->GetVertexData()->GetArray("biconnected_component")->GetMaxId());


  VTK_CREATE(vtkMutableUndirectedGraph, mdgraph0);
  mdgraph0->DeepCopy(bbc->GetOutput());
       

  vtkIntArray* array0= vtkIntArray::SafeDownCast(mdgraph0->GetVertexData()->GetArray("biconnected_component"));

  VTK_CREATE(vtkIdTypeArray, vr_array0);
  VTK_CREATE(vtkIdTypeArray, er_array0);


  VTK_CREATE(vtkVertexListIterator,it0);
  mdgraph0->GetVertices(it0);

  bool ans;
  bool remove_v;
  int n= 0;
  while(ans= it0->HasNext()){ 
    n++;
    vtkIdType nextVertex= it0->Next();

    //std::cout << "Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << mdgraph->GetNumberOfEdgePoints(nextVertex) << " adjacent vertices" << std::endl;  

    VTK_CREATE(vtkAdjacentVertexIterator, itc00);
    mdgraph0->GetAdjacentVertices(nextVertex, itc00);
	
    remove_v= true;
    while (itc00->HasNext()){
      vtkIdType u = itc00->Next();
      int vv= array0->GetValue(u);
      //std::cout << "Checking vertex " << nextCVertex << " with value: " << vv << std::endl;  
      if (vv != COM_INT){
  	remove_v= false;
      }
    }
    if(remove_v){
      //vr_array0->InsertNextValue(nextVertex);
      sel->AddID(0, nextVertex);
      if(VERBOSE) std::cout << "Added vertex " << nextVertex << " to the vertex-remove-list!"<< std::endl; 
    }
  }

  ////coud also use vtkExtractSelectionGraph
  // mdgraph0->RemoveEdges(er_array); //Removes a collection of vertices from t
  // std::cout << "Removed edges!"<< std::endl; 
  //mdgraph0->RemoveVertices(vr_array0); //Removes a collection of vertices from the graph along with any connected edges. AND often screws up edge relations!!!
  //std::cout << "Removed vertices and their edges!"<< std::endl; 

  /////////////////Removing only vertices connected to purely the same kind done.

  sel->SetInverse(1);
  sel->Update();

  VTK_CREATE(vtkExtractSelectedGraph, extract_graph);
  extract_graph->AddInputConnection(bbc->GetOutputPort());
  extract_graph->SetSelectionConnection(sel->GetOutputPort());
  extract_graph->RemoveIsolatedVerticesOff();//If set, removes vertices with no adjacent edges in an edge selection. A vertex selection ignores this flag and always returns the full set of selected vertices. Default is on. 
  extract_graph->Update();


  VTK_CREATE(vtkRemoveIsolatedVertices, riv);
  riv->SetInputConnection(extract_graph->GetOutputPort());
  //  riv->SetInputData(extract_graph->GetOutputData());
  riv->Update();

  // VTK_CREATE(vtkThresholdGraph, bbc_gth);
  // //VTK_CREATE(vtkThreshold, bbc_gth);
  // //VTK_CREATE(vtkThresholdPoints, bbc_gth);
  // bbc_gth->SetInputConnection(bbc->GetOutputPort());
  // bbc_gth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_VERTICES, "biconnected_component");
  // //bbc_gth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_EDGES, "biconnected_component");
  // //bbc_gth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "biconnected_component");
  // //bbc_gth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "biconnected_component");
  // bbc_gth->SetLowerThreshold(0); //all none biconnected components
  // //bbc_gth->SetUpperThreshold(bbc->GetOutput()->GetVertexData()->GetArray("biconnected_component")->GetMaxId());//or GetDataTypeMax()
  // bbc_gth->SetUpperThreshold(bbc->GetOutput()->GetVertexData()->GetArray("biconnected_component")->GetDataTypeMax());//or GetDataTypeMax()
  // //bbc_gth->ThresholdBetween(0, bbc->GetOutput()->GetVertexData()->GetArray("biconnected_component")->GetMaxId());
  // //bbc_gth->ThresholdByUpper(0);
  // //bbc_gth->AllScalarsOff();//does not exist for vtkThresholdGraph?
  // //// because AllScalarsOff() does not exist for vtkThresholdGraph needs to be done manually only removeing vertices that are purely connected to -1
  // bbc_gth->Update();

  // //if (bbc_gth->GetOutput()->GetNumberOfVertices() < 1){
  // // std::cerr  << "bbc_gth->GetOutput()GetNumberOfVertices():" << bbc_gth->GetOutput()->GetNumberOfVertices() << std::endl;
  // // std::cerr  << "bbc_gth->GetOutput()GetNumberOfEdges():" << bbc_gth->GetOutput()->GetNumberOfEdges() << std::endl;
  // //}
  /////////////////Removing none-loop vertices done.


  /////////////////Label each connected-component

  VTK_CREATE(vtkBoostConnectedComponents,bcc);

  //bcc->SetInputConnection(bbc_gth->GetOutputPort());
  //bcc->SetInputData(mdgraph0);
  //bcc->SetInputConnection(extract_graph->GetOutputPort());
  bcc->SetInputConnection(riv->GetOutputPort());

  if(P_VERBOSE) std::cout << "Executing vtkBoostConnectedComponents..." << std::endl;
  bcc->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  //bcc->SetOutputArrayName("component");
  bcc->Update();
  if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

 
  if(!bcc->GetOutput()->GetVertexData()->GetArray("component")){
    std::cerr  << "No array named \"component\" in bcc->GetOutput()!"  << std::endl;
    exit(1);
    //return sub_mesh= NULL;
  }

  /////////////////Label each connected-component done.

  int n_cc;

  // std::cout << VTK_UNSIGNED_CHAR << " unsigned char" << std::endl;
  // std::cout << VTK_UNSIGNED_INT << " unsigned int" << std::endl;
  // std::cout << VTK_FLOAT << " float" << std::endl;
  // std::cout << VTK_DOUBLE << " double" << std::endl;

  // std::cout  << "component array type string: " << bcc->GetOutput()->GetVertexData()->GetArray("component")->GetDataTypeAsString() << std::endl;
  // std::cout  << "component array type: " << bcc->GetOutput()->GetVertexData()->GetArray("component")->GetDataType() << std::endl;
  
  //vtkSmartPointer<vtkDataArray> cc_array= bcc->GetOutput()->GetVertexData()->GetArray("component");
  //vtkSmartPointer<vtkIdTypeArray> cc_array= vtkIdTypeArray::SafeDownCast(bcc->GetOutput()->GetVertexData()->GetArray("component"));
  vtkSmartPointer<vtkIntArray> cc_array= vtkIntArray::SafeDownCast(bcc->GetOutput()->GetVertexData()->GetArray("component"));

  if(cc_array){
    if(VERBOSE) std::cout  << " GetNumberOfComponents() " << cc_array->GetNumberOfComponents() << " cc_array->GetNumberOfTuples() " << cc_array->GetNumberOfTuples() << std::endl;
  }
  else{
    std::cerr  << "No vtkIntArray named \"component\" in bcc->GetOutput()!"  << std::endl;
    n_cc= 0;
    exit(1);
    //return sub_mesh= NULL;
  }

  n_cc= 0;
  for(vtkIdType i = 0; i < cc_array->GetNumberOfTuples(); i++){
    //int v= int(cc_array->GetComponent(i,0));
    int v= cc_array->GetValue(i);
    if(VERBOSE) 
      std::cout  << "entry " << i << " has value: " << v  << std::endl;
    if (v > n_cc){
      n_cc++;
      if(v != n_cc){
	std::cerr  << "entry " << i << " has value: " << v << " != n_cc " << n_cc  << std::endl;
        // ////write test output 
	// vtkSmartPointer<vtkGraphToPolyData> tgraphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
	// tgraphToPolyData->SetInputData(bcc->GetOutput());
	// tgraphToPolyData->Update();
        // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        // twriter->SetFileName("testc.vtp");
        // twriter->SetInputConnection(tgraphToPolyData->GetOutputPort());
        // twriter->Write();
        // std::cout   << "Wrote testc.vtp" << std::endl;
        // ////write test output done.

	exit(1);
      }	
      n_cc++;
    }
  }
  
  if(VERBOSE) std::cout  << "n_cc " << n_cc << std::endl;


  for (vtkIdType i = 0; i <= n_cc; i++){//for each cc

    if(VERBOSE) std::cout  << "doing " << i << std::endl;
    
    VTK_CREATE(vtkThresholdGraph, cc_gth);
    cc_gth->SetInputConnection(bcc->GetOutputPort());
    cc_gth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_VERTICES, "component");
    cc_gth->SetLowerThreshold(i); //all none biconnected components
    cc_gth->SetUpperThreshold(i);
    cc_gth->Update();

    if(VERBOSE) std::cout  << "Extracted cc # " << i << std::endl;

    //////check if cc is a nested loop
    if (isNestedLoop(cc_gth->GetOutput())){/////////////////Fill one loop of a group of nested loops

      if(VERBOSE) std::cout  << "doning " << i << std::endl;

      /////////////////Removing junctions

      vtkSmartPointer<vtkMutableUndirectedGraph> mdgraph=  vtkSmartPointer<vtkMutableUndirectedGraph>::New();
      mdgraph->DeepCopy(cc_gth->GetOutput());
       

      vtkIntArray* array= vtkIntArray::SafeDownCast(mdgraph->GetVertexData()->GetArray("biconnected_component"));

      VTK_CREATE(vtkIdTypeArray, vr_array);
      VTK_CREATE(vtkIdTypeArray, er_array);


      vtkSmartPointer<vtkVertexListIterator> it= vtkSmartPointer<vtkVertexListIterator>::New();
      mdgraph->GetVertices(it);

      bool ans;
      int n= 0;
      while(ans= it->HasNext()){ 
	n++;
	vtkIdType nextVertex= it->Next();
	//std::cout << "Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << mdgraph->GetNumberOfEdgePoints(nextVertex) << " adjacent vertices" << std::endl;  

	vtkSmartPointer<vtkAdjacentVertexIterator> itc= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
	mdgraph->GetAdjacentVertices(nextVertex, itc);
	
	int num_of_adjacent_vertices= 0;
	while (itc->HasNext()){
	  vtkIdType u = itc->Next();
	  num_of_adjacent_vertices++;
	}

	//if (num_of_adjacent_vertices > 2){
	if (array->GetValue(nextVertex) == COM_INT){
	  if(VERBOSE) std::cout << "Pre-Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << num_of_adjacent_vertices << " adjacent vertices" << std::endl;  
	  //if (num_of_adjacent_vertices > 2){ //does not work in cases of -1 not belonging to a triangle, have to check if adacent vertex has one identical adjacent vertex of vertex with v==-1

	  //std::cout << "Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << num_of_adjacent_vertices << " adjacent vertices" << std::endl;  
	  
	  ////checking first "corner"
	  vtkSmartPointer<vtkAdjacentVertexIterator> itc0= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
	  mdgraph->GetAdjacentVertices(nextVertex, itc0);
	  while (itc0->HasNext()){
	    vtkSmartPointer<vtkAdjacentVertexIterator> itc1= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
	    mdgraph->GetAdjacentVertices(itc0->Next(), itc1);
	    while (itc1->HasNext()){
	      vtkSmartPointer<vtkAdjacentVertexIterator> itc2= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
	      mdgraph->GetAdjacentVertices(itc1->Next(), itc2);
	      while (itc2->HasNext()){
		vtkIdType nextCVertex= itc2->Next();
		if (nextCVertex == nextVertex){
		  int vv= array->GetValue(nextCVertex);
		  if(VERBOSE) std::cout << "Checking vertex " << nextCVertex << " with value: " << vv << std::endl;  
		  if (vv == COM_INT){
		    
		    int value_found= 0;
		    for(vtkIdType i= 0; i < vr_array->GetNumberOfTuples(); i++)
		      if (vr_array->GetValue(i) == nextCVertex)
			value_found= 1;
		    if(!value_found){
		      vr_array->InsertNextValue(nextCVertex);
		      if(VERBOSE) std::cout << "Added vertex " << nextCVertex << " with value: " << vv << " to the vertex-remove-list!"<< std::endl; 
		      VTK_CREATE(vtkInEdgeIterator, iei);
		      mdgraph->GetInEdges(nextCVertex, iei);
		      while (iei->HasNext()){
			vtkInEdgeType next_Edge= iei->Next();
			er_array->InsertNextValue(next_Edge.Id);
			if(VERBOSE) std::cout << "Added edge " << next_Edge.Id << " to the edge-remove-list!"<< std::endl; 
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      ////coud also use vtkExtractSelectionGraph
      mdgraph->RemoveEdges(er_array); //Removes a collection of vertices from t
      if(VERBOSE) std::cout << "Removed edges!"<< std::endl; 
      // mdgraph->RemoveVertices(vr_array); //Removes a collection of vertices from the graph along with any connected edges. 
      // std::cout << "Removed vertices and their edges!"<< std::endl; 

      /////////////////Removing junctions done.

      /////////////////changing form vtkGraph to vtkPolyData
      VTK_CREATE(vtkGraphToPolyData, graphToPolyData);
      graphToPolyData->SetInputData(mdgraph);
      
      if(P_VERBOSE) std::cout << "Executing vtkGraphToPolyData..." << std::endl;
      graphToPolyData->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      graphToPolyData->Update();
      if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

      vtkPolyData* mesh= graphToPolyData->GetOutput();

      /////////////////find shortes-loop-path
      ////is there a need to check if sVertex is not a junction???
      vtkIdType sVertex= 0; //try Id==0
      vtkIdType eVertex;

      ////get connected vertices and the cells connecting them
      vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();

      //get all cells that vertex 'id' is a part of
      VTK_CREATE(vtkIdList, cellIdList);
      mesh->BuildLinks();
      
      //sVertex= getEndpointId(mesh);
      if(getEndpointId(mesh)>0){
	std::cout << "mesh still contains end-points! Aborting! " << std::endl;

	vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	twriter->SetFileName("test.vtp");
	twriter->SetInputData(mesh);
	twriter->Write();
	std::cout   << "Wrote test.vtp" << std::endl;
	////write test output done.

	exit(1);
      }
      //      if(sVertex < 0){

      sVertex= -1; 
      do{
	sVertex++;
	mesh->GetPointCells(sVertex, cellIdList); //Make sure that routine BuildLinks() has been called. 
	//} while ((cellIdList->GetNumberOfIds() != 2) || sVertex >= mesh->GetNumberOfPoints());
      } while ((cellIdList->GetNumberOfIds() != 2) || sVertex >= mesh->GetNumberOfPoints());


      if(VERBOSE) std::cout << "sVertex " << sVertex << std::endl;


      //mesh->GetPointCells(sVertex, cellIdList); //Make sure that routine BuildLinks() has been called. 

      if (cellIdList->GetNumberOfIds() != 2){
	std::cerr << "Point not connected to just one or two lines (cells)! # of connected cells: "<< cellIdList->GetNumberOfIds() << " loop splitting embigous!"<< std::endl;
	exit(1);
      }
      /*
	cout << "Vertex 0 is used in cells ";
	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
	cout << cellIdList->GetId(i) << ", ";
	}
	cout << endl;
      */
 
      for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++){
	//cout << "id " << i << " : " << cellIdList->GetId(i) << endl;
 
	vtkSmartPointer<vtkIdList> pointIdList =  vtkSmartPointer<vtkIdList>::New();
	//mesh-> BuildCells(); //not needed
	mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);
	  
	if(VERBOSE) cout << "End-points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;
 
	if(pointIdList->GetId(0) != sVertex){
	  //cout << "Connected to " << pointIdList->GetId(0) << endl;
	  connectedVertices->InsertNextId(pointIdList->GetId(0));
	}
	else{
	  //cout << "Connected to " << pointIdList->GetId(1) << endl;
	  connectedVertices->InsertNextId(pointIdList->GetId(1));
	}
      }

      // eVertex= connectedVertices->GetId(0);

      // if (connectedVertices->GetNumberOfIds() == 2){

      // 	eVertex= connectedVertices->GetId(0);
      // 	std::cout << "Endpoint unembigous! Using: "<< eVertex << std::endl;
    
      // 	////split open line connection between sVertex and eVertex
      // 	mesh->DeleteCell(cellIdList->GetId(0)); //has to correspond to eVertex= connectedVertices->GetId(0);
      // 	mesh->RemoveDeletedCells();
      // 	//mesh->Update();
      // 	std::cout << "Deleted cell between vertex " << sVertex << " and " << eVertex << " : "<< cellIdList->GetId(0) << std::endl;
      // }
      // // }
      // // else{
      // // 	std::cerr << "More than just one connected point (" << connectedVertices->GetNumberOfIds() << "), endpoint embigous!"<< std::endl;
      // // 	exit(1);
      // // }

      eVertex= connectedVertices->GetId(0);
      if(VERBOSE) std::cout << "Endpoint unembigous! Using: "<< eVertex << std::endl;
    
      ////split open line connection between sVertex and eVertex
      mesh->DeleteCell(cellIdList->GetId(0)); //has to correspond to eVertex= connectedVertices->GetId(0);
      mesh->RemoveDeletedCells();
      //mesh->Update();
      if(VERBOSE) std::cout << "Deleted cell between vertex " << sVertex << " and " << eVertex << " : "<< cellIdList->GetId(0) << std::endl;
    

      // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      // twriter->SetFileName("test.vtp");
      // twriter->SetInputData(mesh);
      // twriter->Write();
      // std::cout   << "Wrote test.vtp" << std::endl;
      // ////write test output done.


      VTK_CREATE(vtkDijkstraGraphGeodesicPath, dggp);

      dggp->SetInputData(mesh);
      //dggp->SetInputConnection(extractEdges->GetOutputPort());
      dggp->SetStartVertex(sVertex);
      dggp->SetEndVertex(eVertex);
      dggp->StopWhenEndReachedOn();
      // UseScalarWeightsOn ()
      // RepelPathFromVerticesOn ()
      // SetRepelVertices (vtkPoints *)
  
      if(P_VERBOSE) std::cout << "Executing vtkDijkstraGraphGeodesicPath..." << std::endl;
      dggp->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      dggp->Update();
      if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;


      /////////////////find shortes-loop-path done.
	

      /////////////////fill shortes-loop-path 
      ct->SetInputConnection(dggp->GetOutputPort());
      /////////////////fill shortes-loop-path done.

      // ////write test output 
      // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      // twriter->SetFileName("test_mesh.vtp");
      // twriter->SetInputData(mesh);
      // twriter->Write();
      // std::cout   << "Wrote test_mesh.vtp" << std::endl;
      // ////write test output done.

      std::cout  << "Filled a nested loop." << std::endl;
      sub_mesh= mesh;

    }/////////////////Fill one loop of a group of nested loops done.


    else{/////////////////Fill not nested loop 
  
      VTK_CREATE(vtkGraphToPolyData, graphToPolyData);
      graphToPolyData->SetInputConnection(cc_gth->GetOutputPort());
      
      if(P_VERBOSE) std::cout << "Executing vtkGraphToPolyData..." << std::endl;
      graphToPolyData->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      graphToPolyData->Update();
      if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

      // ////write test output 
      // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      // twriter->SetFileName("test_ccgth.vtp");
      // twriter->SetInputConnection(graphToPolyData->GetOutputPort());
      // twriter->Write();
      // std::cout   << "Wrote test_ccgth.vtp" << std::endl;
      // ////write test output done.

      /////////////////fill loop-path 
      ct->SetInputConnection(graphToPolyData->GetOutputPort());
      /////////////////fill loop-path done.

      std::cout  << "Filled a not nested loop." << std::endl;
      sub_mesh= NULL;

    }/////////////////Fill not nested loop done.
    

    if(P_VERBOSE) std::cout << "Executing vtkContourTriangulator..." << std::endl;
    ct->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    ct->Update();
    if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

    //append->AddInputData(ct->GetOutput());
    append->AddInputConnection(ct->GetOutputPort());

    if(P_VERBOSE) std::cout << "Executing vtkAppendPolydata..." << std::endl;
    append->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    append->Update();
    if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

    connected_filled_loops->DeepCopy(append->GetOutput());

    // ////write test output 
    // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    // twriter->SetFileName("test.vtp");
    // twriter->SetInputData(connected_filled_loops);
    // twriter->Write();
    // std::cout   << "Wrote test.vtp" << std::endl;
    // //exit(1);
    // ////write test output done.


    if(VERBOSE) std::cout << "Filled a loop of cc " << i << std::endl;

  }//for each cc

  return sub_mesh;

}

 
  int main(int argc, char* argv[]){

    if( argc <= 2 ) {
      std::cerr << "Usage: " << argv[0];
      std::cerr << " inputGraph";
      std::cerr << " outputGraph";
      std::cerr << std::endl;  
      return EXIT_FAILURE;
    }

    if(!(strcasestr(argv[1],".vtp"))) {
      std::cout << "The input should end with .vtp" << std::endl; 
      return -1;
    }

    if(!(strcasestr(argv[2],".vtp"))) {
      std::cout << "The output should end with .vtp" << std::endl; 
      return -1;
    }

    VTK_CREATE(vtkCallbackCommand, progressCallback);
    progressCallback->SetCallback(ProgressFunction);


    //vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New(); //*.vtp
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();  

    reader->SetFileName(argv[1]);
    reader->Update();

    ////cleaning should be done before this program is called to make sure the rootID does not change!!!

    // VTK_CREATE(vtkCleanPolyData, cleanFilter);
    // cleanFilter->SetInputConnection(reader->GetOutputPort());
    // cleanFilter->Update();

    // vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
    // triangleFilter->SetInputConnection(cleanFilter->GetOutputPort());
    // triangleFilter->Update();
    
    // vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
    // extractEdges->SetInputConnection(triangleFilter->GetOutputPort());
    // extractEdges->Update();
   
    // std::cout << "Prepared mesh!" << std::endl;

    VTK_CREATE(vtkPolyData, connected_filled_loops);

    vtkSmartPointer<vtkPolyData> sub_loop_mesh;
    sub_loop_mesh= reader->GetOutput();

    int nol= 0;
    while(sub_loop_mesh){
      sub_loop_mesh= fill_loops(sub_loop_mesh, connected_filled_loops);
      nol++;
      std::cout << "fill_loops executed " << nol << " times." << std::endl;
    }//while(sub_loop_mesh)



    vtkSmartPointer<vtkXMLPolyDataWriter> writer= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(argv[2]);
    //writer->SetInputConnection(graphToPolyData->GetOutputPort());
    //writer->SetInputConnection(dggp->GetOutputPort());
    writer->SetInputData(connected_filled_loops);
    writer->Write();
        

    return EXIT_SUCCESS;
  }
