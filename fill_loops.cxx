 

/////program to fill loops according to shortes-loop-paths combining graph_BiConnectedComponents_Boost remove_graph-junctions_01 shortest-loop-path_03 and  vtkContourTriangulator
//02: recursive version

#include <vtkSmartPointer.h>


#include <vtkXMLPolyDataReader.h>//for vtp-files (cannot contain 3D cells?)

//#include <vtkCleanPolyData.h>
//#include <vtkTriangleFilter.h>
//#include <vtkExtractEdges.h>

#include <vtkPolyDataToGraph.h>
#include <vtkBoostBiconnectedComponents.h>

#include <vtkThresholdGraph.h>

#include <vtkBoostConnectedComponents.h>

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

#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()



void ProgressFunction( vtkObject* caller, long unsigned int eventId, void* clientData, void* callData ){
    vtkAlgorithm *d= static_cast<vtkAlgorithm*>(caller);
    fprintf(stderr, "\rFilter progress: %5.1f%%", 100.0 * d->GetProgress());
    std::cerr.flush();
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

  std::cout << "Executing vtkPolyDataToGraph..." << std::endl;
  polyDataToGraphFilter->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  polyDataToGraphFilter->Update();
  std::cout  << std::endl << "done." << std::endl;



  //vtkGraph* graph= polyDataToGraphFilter->GetOutput();


  VTK_CREATE(vtkBoostBiconnectedComponents, bbc);
  bbc->SetInputConnection(polyDataToGraphFilter->GetOutputPort());
  bbc->SetOutputArrayName("biconnected_component");

  std::cout << "Executing vtkBoostBiconnectedComponents..." << std::endl;
  bbc->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  bbc->Update();
  std::cout  << std::endl << "done." << std::endl;



  /////////////////Identifying loops done.


  /////////////////Removing none-loop vertices

  VTK_CREATE(vtkThresholdGraph, bbc_gth);
  bbc_gth->SetInputConnection(bbc->GetOutputPort());
  bbc_gth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_EDGES, "biconnected_component");
  bbc_gth->SetLowerThreshold(0); //all none biconnected components
  bbc_gth->SetUpperThreshold(bbc->GetOutput()->GetVertexData()->GetArray("biconnected_component")->GetMaxId());//or GetDataTypeMax()
  //bbc_gth->AllScalarsOff();//does not exist for vtkThresholdGraph?
  bbc_gth->Update();

  //if (bbc_gth->GetOutput()->GetNumberOfVertices() < 1){
  std::cerr  << "bbc_gth->GetOutput()GetNumberOfVertices():" << bbc_gth->GetOutput()->GetNumberOfVertices() << std::endl;
  std::cerr  << "bbc_gth->GetOutput()GetNumberOfEdges():" << bbc_gth->GetOutput()->GetNumberOfEdges() << std::endl;
  //}
    

  /////////////////Removing none-loop vertices done.


  /////////////////Label each connected-component

  VTK_CREATE(vtkBoostConnectedComponents,bcc);

  bcc->SetInputConnection(bbc_gth->GetOutputPort());

  std::cout << "Executing vtkBoostConnectedComponents..." << std::endl;
  bcc->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  //bcc->SetOutputArrayName("component");
  bcc->Update();
  std::cout  << std::endl << "done." << std::endl;

  // ////write test output
  // vtkSmartPointer<vtkGraphToPolyData> tgraphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
  // tgraphToPolyData->SetInputData(bcc->GetOutput());
  // tgraphToPolyData->Update();
      
  // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  // twriter->SetFileName("test.vtp");
  // twriter->SetInputConnection(tgraphToPolyData->GetOutputPort());
  // twriter->Write();
  // ////write test output done.

  if(!bcc->GetOutput()->GetVertexData()->GetArray("component")){
    std::cerr  << "No array named \"component\" in bcc->GetOutput()!"  << std::endl;
    exit(1);
    //return sub_mesh= NULL;
  }

  /////////////////Label each connected-component done.

  int n_cc;

  vtkSmartPointer<vtkDataArray> cc_array= bcc->GetOutput()->GetVertexData()->GetArray("component");
  //vtkSmartPointer<vtkIdTypeArray> cc_array= vtkIdTypeArray::SafeDownCast(bcc->GetOutput()->GetVertexData()->GetArray("component"));

  if(cc_array){
    std::cerr  << " GetNumberOfComponents() " << cc_array->GetNumberOfComponents() << " cc_array->GetNumberOfTuples() " << cc_array->GetNumberOfTuples() << std::endl;
  }
  else{
    std::cerr  << "No vtkIdTypeArray array named \"component\" in bcc->GetOutput()!"  << std::endl;
    n_cc= 0;
    exit(1);
    //return sub_mesh= NULL;
  }

  for(vtkIdType i = 0; i < cc_array->GetNumberOfTuples(); i++){
    int v= int(cc_array->GetComponent(i,0));
    if (v >= n_cc){
      if(v != n_cc){
	std::cerr  << "v (" << v <<") != n_cc " << n_cc  << std::endl;
	exit(1);
      }	
      n_cc++;
    }
  }
  
  std::cerr  << "n_cc " << n_cc << std::endl;


  for (vtkIdType i = 0; i < n_cc; i++){//for each cc

    std::cerr  << "doning " << i << std::endl;
    
    VTK_CREATE(vtkThresholdGraph, cc_gth);
    cc_gth->SetInputConnection(bcc->GetOutputPort());
    cc_gth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_VERTICES, "component");
    cc_gth->SetLowerThreshold(i); //all none biconnected components
    cc_gth->SetUpperThreshold(i);
    cc_gth->Update();

    std::cerr  << "Extracted cc # " << i << std::endl;

    //////check if cc is a nested loop
    if (isNestedLoop(cc_gth->GetOutput())){/////////////////Fill one loop of a group of nested loops

      std::cerr  << "doning " << i << std::endl;

      /////////////////Removing junctions

      vtkSmartPointer<vtkMutableUndirectedGraph> mdgraph=  vtkSmartPointer<vtkMutableUndirectedGraph>::New();
      mdgraph->DeepCopy(cc_gth->GetOutput());
       

      vtkIntArray* array= vtkIntArray::SafeDownCast(mdgraph->GetVertexData()->GetArray("biconnected_component"));

      VTK_CREATE(vtkIdTypeArray, vr_array);
      VTK_CREATE(vtkIdTypeArray, er_array);


      int com_int= -1; //none biconnected components
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
	if (array->GetValue(nextVertex) == com_int){
	  std::cout << "Pre-Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << num_of_adjacent_vertices << " adjacent vertices" << std::endl;  
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
		  std::cout << "Checking vertex " << nextCVertex << " with value: " << vv << std::endl;  
		  if (vv == com_int){
		    
		    int value_found= 0;
		    for(vtkIdType i= 0; i < vr_array->GetNumberOfTuples(); i++)
		      if (vr_array->GetValue(i) == nextCVertex)
			value_found= 1;
		    if(!value_found){
		      vr_array->InsertNextValue(nextCVertex);
		      std::cout << "Added vertex " << nextCVertex << " with value: " << vv << " to the vertex-remove-list!"<< std::endl; 
		      VTK_CREATE(vtkInEdgeIterator, iei);
		      mdgraph->GetInEdges(nextCVertex, iei);
		      while (iei->HasNext()){
			vtkInEdgeType next_Edge= iei->Next();
			er_array->InsertNextValue(next_Edge.Id);
			std::cout << "Added edge " << next_Edge.Id << " to the edge-remove-list!"<< std::endl; 
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
      std::cout << "Removed edges!"<< std::endl; 
      // mdgraph->RemoveVertices(vr_array); //Removes a collection of vertices from the graph along with any connected edges. 
      // std::cout << "Removed vertices and their edges!"<< std::endl; 

      /////////////////Removing junctions done.

      /////////////////changing form vtkGraph to vtkPolyData
      VTK_CREATE(vtkGraphToPolyData, graphToPolyData);
      graphToPolyData->SetInputData(mdgraph);
      
      std::cout << "Executing vtkGraphToPolyData..." << std::endl;
      graphToPolyData->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      graphToPolyData->Update();
      std::cout  << std::endl << "done." << std::endl;

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
	exit(1);
      }
      //      if(sVertex < 0){

      sVertex= -1; 
      do{
	sVertex++;
	mesh->GetPointCells(sVertex, cellIdList); //Make sure that routine BuildLinks() has been called. 
	//} while ((cellIdList->GetNumberOfIds() != 2) || sVertex >= mesh->GetNumberOfPoints());
      } while ((cellIdList->GetNumberOfIds() != 2) || sVertex >= mesh->GetNumberOfPoints());


      std::cout << "sVertex " << sVertex << std::endl;


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
	  
	if(VERBOSE)
	  cout << "End-points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;
 
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
      std::cout << "Endpoint unembigous! Using: "<< eVertex << std::endl;
    
      ////split open line connection between sVertex and eVertex
      mesh->DeleteCell(cellIdList->GetId(0)); //has to correspond to eVertex= connectedVertices->GetId(0);
      mesh->RemoveDeletedCells();
      //mesh->Update();
      std::cout << "Deleted cell between vertex " << sVertex << " and " << eVertex << " : "<< cellIdList->GetId(0) << std::endl;
    

      vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      twriter->SetFileName("test.vtp");
      twriter->SetInputData(mesh);
      twriter->Write();
      std::cout   << "Wrote test.vtp" << std::endl;
      ////write test output done.


      VTK_CREATE(vtkDijkstraGraphGeodesicPath, dggp);

      dggp->SetInputData(mesh);
      //dggp->SetInputConnection(extractEdges->GetOutputPort());
      dggp->SetStartVertex(sVertex);
      dggp->SetEndVertex(eVertex);
      dggp->StopWhenEndReachedOn();
      // UseScalarWeightsOn ()
      // RepelPathFromVerticesOn ()
      // SetRepelVertices (vtkPoints *)
  
      std::cout << "Executing vtkDijkstraGraphGeodesicPath..." << std::endl;
      dggp->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      dggp->Update();
      std::cout  << std::endl << "done." << std::endl;


      /////////////////find shortes-loop-path done.
	

      /////////////////fill shortes-loop-path 
      ct->SetInputConnection(dggp->GetOutputPort());
      /////////////////fill shortes-loop-path done.

      ////write test output 
      // vtkSmartPointer<vtkGraphToPolyData> tgraphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
      // tgraphToPolyData->SetInputData(bcc->GetOutput());
      // tgraphToPolyData->Update();
      // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      // twriter->SetFileName("test.vtp");
      // twriter->SetInputConnection(dggp->GetOutputPort());
      // twriter->Write();
      // std::cout   << "Wrote test.vtp" << std::endl;
      // ////write test output done.

	sub_mesh= mesh;
    }/////////////////Fill one loop of a group of nested loops done.


    else{/////////////////Fill not nested loop 
  
      VTK_CREATE(vtkGraphToPolyData, graphToPolyData);
      graphToPolyData->SetInputConnection(cc_gth->GetOutputPort());
      
      std::cout << "Executing vtkGraphToPolyData..." << std::endl;
      graphToPolyData->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      graphToPolyData->Update();
      std::cout  << std::endl << "done." << std::endl;

      //   ////write test output 
      // // vtkSmartPointer<vtkGraphToPolyData> tgraphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
      // // tgraphToPolyData->SetInputData(bcc->GetOutput());
      // // tgraphToPolyData->Update();
      //   vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      //   twriter->SetFileName("test.vtp");
      //   twriter->SetInputConnection(graphToPolyData->GetOutputPort());
      //   twriter->Write();
      //   std::cout   << "Wrote test.vtp" << std::endl;
      //   ////write test output done.

      /////////////////fill loop-path 
      ct->SetInputConnection(graphToPolyData->GetOutputPort());
      /////////////////fill loop-path done.
      sub_mesh= NULL;
    }/////////////////Fill not nested loop done.
    

    std::cout << "Executing vtkContourTriangulator..." << std::endl;
    ct->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    ct->Update();
    std::cout  << std::endl << "done." << std::endl;

    //append->AddInputData(ct->GetOutput());
    append->AddInputConnection(ct->GetOutputPort());

    std::cout << "Executing vtkAppendPolydata..." << std::endl;
    append->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    append->Update();
    std::cout  << std::endl << "done." << std::endl;

    connected_filled_loops= append->GetOutput();

    std::cout << "Filled a loop of cc " << i << std::endl;

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

    vtkSmartPointer<vtkPolyData> connected_filled_loops= vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPolyData> sub_loop_mesh;
    sub_loop_mesh= reader->GetOutput();

    int nol= 0;
    while(sub_loop_mesh){
      sub_loop_mesh= fill_loops(sub_loop_mesh, connected_filled_loops);
      nol++;
      std::cout << "Filled " << nol << " loops" << std::endl;
    }//while(sub_loop_mesh)



    vtkSmartPointer<vtkXMLPolyDataWriter> writer= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(argv[2]);
    //writer->SetInputConnection(graphToPolyData->GetOutputPort());
    //writer->SetInputConnection(dggp->GetOutputPort());
    writer->SetInputData(connected_filled_loops);
    writer->Write();
        

    return EXIT_SUCCESS;
  }
