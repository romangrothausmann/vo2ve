 

/////program to fill loops according to shortes-loop-paths combining graph_BiConnectedComponents_Boost remove_graph-junctions_01 shortest-loop-path_03 and  vtkContourTriangulator
//02: recursive version
//03: debugging
//04: replace vtkContourTriangulator by vtkTrianglulate of an n-gon
//05: switch from VERTEX ID-selection to EDGE threshold-selections
//{06: append performance improvements (not working yet)}

//todo: 
// - compute Euler Characteristics to optain total numer of loops before head and for total progress report;-)
// - save skippted graphes too for manual checking
// - use openMP for main for-loop

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

//#include <vtkContourTriangulator.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>
#include <vtkStripper.h> //to convert vtkLines to vtkPolyLine
#include <vtkPolyLine.h>
#include <vtkTriangleFilter.h> 
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
#define P_VERBOSE 0
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

  mesh->BuildLinks();
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
  //VTK_CREATE(vtkPolyData, sub_mesh);

  VTK_CREATE(vtkAppendPolyData, append_fl);
  append_fl->AddInputData(connected_filled_loops);

  //VTK_CREATE(vtkAppendPolyData, append_rm);

  //vtkSmartPointer<vtkPolyData> mesh;
  // VTK_CREATE(vtkPolyData, pmesh);
  // VTK_CREATE(vtkPolygon, polygon);
  // VTK_CREATE(vtkCellArray, outCells);

  //VTK_CREATE(vtkContourTriangulator, ct);
  //ct->TriangulationErrorDisplayOn();
  //VTK_CREATE(vtkTriangleFilter, tf);

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



  /////////////////Removing only vertices connected to purely the same kind, i.e. bbc==-1

    
  VTK_CREATE(vtkSelectionSource, sel);
  sel->SetContentType(vtkSelectionNode::THRESHOLDS); 
  //sel->SetContentType(vtkSelectionNode::INDICES);
  //sel->SetFieldType(vtkSelectionNode::VERTEX); 
  sel->SetFieldType(vtkSelectionNode::EDGE); 
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

    int nextVertex_v= array0->GetValue(nextVertex);
    //std::cout << "Checking vertex " << nextVertex << " with value: " << nextVertex_v << " which has " << mdgraph->GetNumberOfEdgePoints(nextVertex) << " adjacent vertices" << std::endl;  

    VTK_CREATE(vtkAdjacentVertexIterator, itc00);
    mdgraph0->GetAdjacentVertices(nextVertex, itc00);

    int non= 0; 
    while (itc00->HasNext()){
      non++;
      vtkIdType u = itc00->Next();
    }
    //std::cout << "Vertex " << nextVertex << " with value: " << nextVertex_v << " has a non of: " << non <<std::endl;  
    // if(non < 2){
    //   remove_v= true;
    //   if(VERBOSE) std::cout << "Removing vertex " << nextVertex << " with value: " << nextVertex_v << " which has less than 2 non." <<std::endl;  
    // }

    if ((nextVertex_v != COM_INT)&&(non > 1)){
      //sel->AddID(0, nextVertex);
      sel->AddThreshold(nextVertex_v, nextVertex_v);
      if(VERBOSE) std::cout << "Added " << nextVertex_v << " to the edge-selection thresholds!"<< std::endl; 
    }


    // VTK_CREATE(vtkAdjacentVertexIterator, itc00);
    // mdgraph0->GetAdjacentVertices(nextVertex, itc00);
    
    // remove_v= false;
    // if(nextVertex_v == COM_INT){
    //   remove_v= true;
    // }
    // int non= 0; 
    // while (itc00->HasNext()){
    //   non++;
    //   vtkIdType u = itc00->Next();
    //   int vv= array0->GetValue(u);
    //   //std::cout << "Checking vertex " << u << " with value: " << vv << std::endl;  
    //   if ((nextVertex_v == COM_INT) && (vv != COM_INT)){//only apply for nextVertex_v == -1
    // 	remove_v= false;
    //   }
    // }
    // //std::cout << "Vertex " << nextVertex << " with value: " << nextVertex_v << " has a non of: " << non <<std::endl;  
    // if(non < 2){
    //   remove_v= true;
    //   if(VERBOSE) std::cout << "Removing vertex " << nextVertex << " with value: " << nextVertex_v << " which has less than 2 non." <<std::endl;  
    // }
    // if(remove_v){
    //   //vr_array0->InsertNextValue(nextVertex);
    //   sel->AddID(0, nextVertex);
    //   if(VERBOSE) std::cout << "Added vertex " << nextVertex << " to the vertex-remove-list!"<< std::endl; 
    // }
  }

  ////coud also use vtkExtractSelectionGraph
  // mdgraph0->RemoveEdges(er_array); //Removes a collection of vertices from t
  // std::cout << "Removed edges!"<< std::endl; 
  //mdgraph0->RemoveVertices(vr_array0); //Removes a collection of vertices from the graph along with any connected edges. AND often screws up edge relations!!!
  //std::cout << "Removed vertices and their edges!"<< std::endl; 



  //sel->SetInverse(1);
  sel->Update();

  VTK_CREATE(vtkExtractSelectedGraph, extract_graph);
  extract_graph->AddInputConnection(bbc->GetOutputPort());
  extract_graph->SetSelectionConnection(sel->GetOutputPort());
  extract_graph->RemoveIsolatedVerticesOff();//If set, removes vertices with no adjacent edges in an edge selection. A vertex selection ignores this flag and always returns the full set of selected vertices. Default is on. 
  //if(P_VERBOSE) 
  std::cout << "Executing vtkExtractSelectedGraph..." << std::endl;
  //extract_graph->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  extract_graph->Update();
  //if(P_VERBOSE) 
  std::cout  << std::endl << "done." << std::endl;

  // ////write test output 
  // vtkSmartPointer<vtkGraphToPolyData> tgraphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
  // tgraphToPolyData->SetInputData(extract_graph->GetOutput());
  // tgraphToPolyData->Update();
  // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  // twriter->SetFileName("teseg.vtp");
  // twriter->SetInputConnection(tgraphToPolyData->GetOutputPort());
  // twriter->Write();
  // std::cout   << "Wrote testeg.vtp" << std::endl;
  // ////write test output done.


  /////////////////Removing only vertices connected to purely the same kind done.

  VTK_CREATE(vtkRemoveIsolatedVertices, riv);
  riv->SetInputConnection(extract_graph->GetOutputPort());
  //  riv->SetInputData(extract_graph->GetOutputData());
  riv->Update();


  // ////write test output 
  // vtkSmartPointer<vtkGraphToPolyData> tgraphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
  // tgraphToPolyData->SetInputData(riv->GetOutput());
  // tgraphToPolyData->Update();
  // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  // twriter->SetFileName("testriv.vtp");
  // twriter->SetInputConnection(tgraphToPolyData->GetOutputPort());
  // twriter->Write();
  // std::cout   << "Wrote testriv.vtp" << std::endl;
  // ////write test output done.


  /////////////////Removing none-loop vertices done.


  /////////////////Label each connected-component

  VTK_CREATE(vtkBoostConnectedComponents,bcc);

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
    if(VERBOSE) std::cout  << "entry " << i << " has value: " << v  << std::endl;
    if (v > n_cc){
      n_cc++;
      if(v != n_cc){
	std::cerr  << "entry " << i << " has value: " << v << " != n_cc " << n_cc  << std::endl;
        ////write test output 
	vtkSmartPointer<vtkGraphToPolyData> tgraphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
	tgraphToPolyData->SetInputData(bcc->GetOutput());
	tgraphToPolyData->Update();
        vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        twriter->SetFileName("testc.vtp");
        twriter->SetInputConnection(tgraphToPolyData->GetOutputPort());
        twriter->Write();
        std::cout   << "Wrote testc.vtp" << std::endl;
        ////write test output done.
	exit(1);
      }	
      n_cc++;
    }
  }
  
  if(VERBOSE) std::cout  << "n_cc " << n_cc << std::endl;


  for (vtkIdType i = 0; i <= n_cc; i++){//for each cc

    //if(VERBOSE) 
    //std::cout  << "doing cc: " << i << " of " << n_cc << "[" << i*100/n_cc << "%]" std::endl;
    printf("Doing cc: %d of %d [%5.1f%%]\n", i+1 , n_cc+1,  (i+1)*100.0/(n_cc+1));
 
    VTK_CREATE(vtkPolyData, pmesh);
    VTK_CREATE(vtkPolygon, polygon);
    VTK_CREATE(vtkCellArray, outCells);
    VTK_CREATE(vtkTriangleFilter, tf);


   
    VTK_CREATE(vtkThresholdGraph, cc_gth);
    cc_gth->SetInputConnection(bcc->GetOutputPort());
    cc_gth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_VERTICES, "component");
    cc_gth->SetLowerThreshold(i); //all none biconnected components
    cc_gth->SetUpperThreshold(i);
    cc_gth->Update();

    if(VERBOSE) std::cout  << "Extracted cc # " << i << std::endl;


    VTK_CREATE(vtkMutableUndirectedGraph, mdgraph);
    mdgraph->DeepCopy(cc_gth->GetOutput());
       

    vtkSmartPointer<vtkIntArray> array= vtkIntArray::SafeDownCast(mdgraph->GetVertexData()->GetArray("biconnected_component"));

    VTK_CREATE(vtkVertexListIterator, it);
    mdgraph->GetVertices(it);

    bool ans;
    int n= 0;
    vtkIdType ccv;
    while(ans= it->HasNext()){ 
      vtkIdType nextVertex= it->Next();

      if(array->GetValue(nextVertex) != COM_INT){
	ccv= array->GetValue(nextVertex);
	n++;
      }
      //std::cout << "Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << mdgraph->GetNumberOfEdgePoints(nextVertex) << " adjacent vertices" << std::endl;  
    }
    //if(n < 3){
    //if(n < 6){ //6 to avoid cross in square graph
    if(n < 10){ //10 to avoid small but odd structures such as in 13-621_CO_LW-seg_03_11_01.png
      std::cout << "Skipping cc " << i << " which has less than 3 vertices with value: " << ccv << std::endl;  
      continue;
    }

    //////check if cc is a nested loop
    if (isNestedLoop(cc_gth->GetOutput())){/////////////////Fill one loop of a group of nested loops

      if(VERBOSE) std::cout  << "doning " << i << std::endl;

      // /////////////////Removing junctions

      // vtkSmartPointer<vtkMutableUndirectedGraph> mdgraph=  vtkSmartPointer<vtkMutableUndirectedGraph>::New();
      // mdgraph->DeepCopy(cc_gth->GetOutput());
       

      // vtkIntArray* array= vtkIntArray::SafeDownCast(mdgraph->GetVertexData()->GetArray("biconnected_component"));

      // VTK_CREATE(vtkIdTypeArray, vr_array);
      // VTK_CREATE(vtkIdTypeArray, er_array);


      // vtkSmartPointer<vtkVertexListIterator> it= vtkSmartPointer<vtkVertexListIterator>::New();
      // mdgraph->GetVertices(it);

      // bool ans;
      // int n= 0;
      // while(ans= it->HasNext()){ 
      // 	n++;
      // 	vtkIdType nextVertex= it->Next();
      // 	//std::cout << "Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << mdgraph->GetNumberOfEdgePoints(nextVertex) << " adjacent vertices" << std::endl;  

      // 	vtkSmartPointer<vtkAdjacentVertexIterator> itc= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
      // 	mdgraph->GetAdjacentVertices(nextVertex, itc);
	
      // 	int num_of_adjacent_vertices= 0;
      // 	while (itc->HasNext()){
      // 	  vtkIdType u = itc->Next();
      // 	  num_of_adjacent_vertices++;
      // 	}

      // 	// //if(nextVertex == 622)
      // 	// std::cout << "Vertex " << nextVertex << ", num_of_adjacent_vertices " << num_of_adjacent_vertices << std::endl; 
      // 	// if (num_of_adjacent_vertices < 2){
	    
      // 	//   vr_array->InsertNextValue(nextVertex);
      // 	//   //if(VERBOSE) 
      // 	//   std::cout << "Added vertex " << nextVertex << " (which has less than two neighbours) to the vertex-remove-list!"<< std::endl; 
      // 	//   VTK_CREATE(vtkInEdgeIterator, iei);
      // 	//   mdgraph->GetInEdges(nextVertex, iei);
      // 	//   while (iei->HasNext()){
      // 	//     vtkInEdgeType next_Edge= iei->Next();
      // 	//     er_array->InsertNextValue(next_Edge.Id);
      // 	//     //if(VERBOSE) 
      // 	//     std::cout << "Added edge " << next_Edge.Id << " to the edge-remove-list!"<< std::endl; 
      // 	//   }
      // 	// }

      // 	if (array->GetValue(nextVertex) == COM_INT){
      // 	  if(VERBOSE) std::cout << "Pre-Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << num_of_adjacent_vertices << " adjacent vertices" << std::endl;  
      // 	  //if (num_of_adjacent_vertices > 2){ //does not work in cases of -1 not belonging to a triangle, have to check if adacent vertex has one identical adjacent vertex of vertex with v==-1

      // 	  //std::cout << "Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << num_of_adjacent_vertices << " adjacent vertices" << std::endl;  
	  
      // 	  if (num_of_adjacent_vertices < 2){
	    
      // 	    vr_array->InsertNextValue(nextVertex);
      // 	    //if(VERBOSE) 
      // 	    std::cout << "Added vertex (which has less than two neighbours) " << nextVertex << " to the vertex-remove-list!"<< std::endl; 
      // 	    VTK_CREATE(vtkInEdgeIterator, iei);
      // 	    mdgraph->GetInEdges(nextVertex, iei);
      // 	    while (iei->HasNext()){
      // 	      vtkInEdgeType next_Edge= iei->Next();
      // 	      er_array->InsertNextValue(next_Edge.Id);
      // 	      //if(VERBOSE) 
      // 	      std::cout << "Added edge " << next_Edge.Id << " to the edge-remove-list!"<< std::endl; 
      // 	    }
      // 	  }
      // 	  else if (num_of_adjacent_vertices == 2){
      // 	    vtkSmartPointer<vtkAdjacentVertexIterator> itc0= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
      // 	    mdgraph->GetAdjacentVertices(nextVertex, itc0);
      // 	    int av= array->GetValue(itc0->Next());
      // 	    int bv= array->GetValue(itc0->Next());
      // 	    if (itc0->HasNext()){
      // 	      std::cerr << "Vertex " << nextVertex << " has more than 2 neighbours here. Aborting!"<< std::endl;
      // 	      exit(1);
      // 	    }
      // 	    //std::cout << "Vertex " << nextVertex << " has 2 neighbours, av " << av << ", bv " << bv << std::endl;
      // 	    if (av != bv){
      // 	      vr_array->InsertNextValue(nextVertex);
      // 	      //if(VERBOSE) 
      // 	      std::cout << "Added vertex " << nextVertex << " to the vertex-remove-list!"<< std::endl; 
      // 	      VTK_CREATE(vtkInEdgeIterator, iei);
      // 	      mdgraph->GetInEdges(nextVertex, iei);
      // 	      while (iei->HasNext()){
      // 		vtkInEdgeType next_Edge= iei->Next();
      // 		er_array->InsertNextValue(next_Edge.Id);
      // 		if(VERBOSE) std::cout << "Added edge " << next_Edge.Id << " to the edge-remove-list!"<< std::endl; 
      // 	      }
      // 	    }
      // 	  }//else if (num_of_adjacent_vertices == 2){
      // 	  else{

      // 	    ////checking first "corner"
      // 	    vtkSmartPointer<vtkAdjacentVertexIterator> itc0= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
      // 	    mdgraph->GetAdjacentVertices(nextVertex, itc0);
      // 	    while (itc0->HasNext()){
      // 	      vtkSmartPointer<vtkAdjacentVertexIterator> itc1= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
      // 	      mdgraph->GetAdjacentVertices(itc0->Next(), itc1);
      // 	      while (itc1->HasNext()){
      // 		vtkSmartPointer<vtkAdjacentVertexIterator> itc2= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
      // 		mdgraph->GetAdjacentVertices(itc1->Next(), itc2);
      // 		while (itc2->HasNext()){
      // 		  vtkIdType nextCVertex= itc2->Next();
      // 		  if (nextCVertex == nextVertex){
      // 		    int vv= array->GetValue(nextCVertex);
      // 		    if(VERBOSE) std::cout << "Checking vertex " << nextCVertex << " with value: " << vv << std::endl;  
      // 		    if (vv == COM_INT){
		    
      // 		      int value_found= 0;
      // 		      for(vtkIdType i= 0; i < vr_array->GetNumberOfTuples(); i++)
      // 			if (vr_array->GetValue(i) == nextCVertex)
      // 			  value_found= 1;
      // 		      if(!value_found){
      // 			vr_array->InsertNextValue(nextCVertex);
      // 			if(VERBOSE) std::cout << "Added vertex " << nextCVertex << " with value: " << vv << " to the vertex-remove-list!"<< std::endl; 
      // 			VTK_CREATE(vtkInEdgeIterator, iei);
      // 			mdgraph->GetInEdges(nextCVertex, iei);
      // 			while (iei->HasNext()){
      // 			  vtkInEdgeType next_Edge= iei->Next();
      // 			  er_array->InsertNextValue(next_Edge.Id);
      // 			  if(VERBOSE) std::cout << "Added edge " << next_Edge.Id << " to the edge-remove-list!"<< std::endl; 
      // 			}
      // 		      }
      // 		    }
      // 		  }
      // 		}
      // 	      }
      // 	    }
      // 	  }//if (num_of_adjacent_vertices < 2){}else{
      // 	}
      // }

      // ////coud also use vtkExtractSelectionGraph
      // mdgraph->RemoveEdges(er_array); //Removes a collection of vertices from t
      // if(VERBOSE) std::cout << "Removed edges!"<< std::endl; 
      // // mdgraph->RemoveVertices(vr_array); //Removes a collection of vertices from the graph along with any connected edges. 
      // // std::cout << "Removed vertices and their edges!"<< std::endl; 

      // /////////////////Removing junctions done.

      // VTK_CREATE(vtkRemoveIsolatedVertices, riv2);
      // riv2->SetInputData(mdgraph);
      // //  riv->SetInputData(extract_graph->GetOutputData());
      // riv2->Update();

      /////////////////changing form vtkGraph to vtkPolyData
      VTK_CREATE(vtkGraphToPolyData, graphToPolyData);
      //graphToPolyData->SetInputData(mdgraph);
      //graphToPolyData->SetInputConnection(riv2->GetOutputPort());
      graphToPolyData->SetInputConnection(cc_gth->GetOutputPort());

      if(P_VERBOSE) std::cout << "Executing vtkGraphToPolyData for nested loop..." << std::endl;
      graphToPolyData->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      graphToPolyData->Update();
      if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

      vtkSmartPointer<vtkPolyData> rmesh= graphToPolyData->GetOutput();

      // ////write test output 
      // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      // twriter->SetFileName("test_rj.vtp");
      // twriter->SetInputData(rmesh);
      // twriter->Write();
      // std::cout   << "Wrote test_rj.vtp" << std::endl;
      // ////write test output done.

      /////////////////find shortes-loop-path
      ////is there a need to check if sVertex is not a junction???
      vtkIdType sVertex= 0; //try Id==0
      vtkIdType eVertex;

      ////get connected vertices and the cells connecting them
      VTK_CREATE(vtkIdList, connectedVertices);

      //get all cells that vertex 'id' is a part of
      VTK_CREATE(vtkIdList, cellIdList);
      if(VERBOSE) std::cout << "Building Links... " << std::flush << std::endl;
      rmesh->BuildLinks();
      
      vtkIdType eid= getEndpointId(rmesh);
      if(eid>0){
	std::cerr << "rmesh still contains end-points! Aborting! eid: " << eid << std::endl;

	vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	twriter->SetFileName("test.vtp");
	twriter->SetInputData(rmesh);
	twriter->Write();
	std::cout   << "Wrote test.vtp" << std::endl;
	////write test output done.

	exit(1);
      }
      //      if(sVertex < 0){

      sVertex= -1; 
      do{
	sVertex++;
	if(VERBOSE) std::cout << "sVertex " << sVertex << std::endl;
	rmesh->GetPointCells(sVertex, cellIdList); //Make sure that routine BuildLinks() has been called. 
	//} while ((cellIdList->GetNumberOfIds() != 2) || sVertex >= rmesh->GetNumberOfPoints());
      } while ((cellIdList->GetNumberOfIds() != 2) || sVertex >= rmesh->GetNumberOfPoints());


      if(VERBOSE) std::cout << "sVertex " << sVertex << std::endl;


      //rmesh->GetPointCells(sVertex, cellIdList); //Make sure that routine BuildLinks() has been called. 

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
	//rmesh-> BuildCells(); //not needed
	rmesh->GetCellPoints(cellIdList->GetId(i), pointIdList);
	  
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
      // 	rmesh->DeleteCell(cellIdList->GetId(0)); //has to correspond to eVertex= connectedVertices->GetId(0);
      // 	rmesh->RemoveDeletedCells();
      // 	//rmesh->Update();
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
      rmesh->DeleteCell(cellIdList->GetId(0)); //has to correspond to eVertex= connectedVertices->GetId(0);
      rmesh->RemoveDeletedCells();
      //rmesh->Update();
      if(VERBOSE) std::cout << "Deleted cell between vertex " << sVertex << " and " << eVertex << " : "<< cellIdList->GetId(0) << std::endl;
    
      // ////write test output 
      // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      // twriter->SetFileName("testrmesh.vtp");
      // twriter->SetInputData(rmesh);
      // twriter->Write();
      // std::cout   << "Wrote testrmesh.vtp" << std::endl;
      // ////write test output done.


      VTK_CREATE(vtkDijkstraGraphGeodesicPath, dggp);

      dggp->SetInputData(rmesh);
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

      /////////////////create n-gon

      // VTK_CREATE(vtkPolyData, pmesh);
      // VTK_CREATE(vtkPolygon, polygon);
      // VTK_CREATE(vtkCellArray, outCells);

      VTK_CREATE(vtkCleanPolyData, cpd);
      cpd->SetInputConnection(dggp->GetOutputPort());
      cpd->PointMergingOn();
      VTK_CREATE(vtkStripper, stripper);
      //stripper->SetInputConnection(dggp->GetOutputPort());
      stripper->SetInputConnection(cpd->GetOutputPort());

      if(P_VERBOSE) std::cout << "Executing vtkStripper..." << std::endl;
      stripper->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      stripper->Update();
      if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;


      //pmesh= dggp->GetOutput();
      //pmesh->DeepCopy(graphToPolyData->GetOutput());
      pmesh->DeepCopy(stripper->GetOutput());

      // VTK_CREATE(vtkPolygon, polygon);
      //polygon->GetPoints()->DeepCopy(mesh); //only good if n-gon is use without assignment to vtkPolyData: http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolygonIntersection
      // polygon->SetPoints(mesh->GetPoints());
      // polygon->SetPointIds(mesh->GetPointIds());
      vtkIdType mnop= pmesh->GetPoints()->GetNumberOfPoints();

      vtkSmartPointer<vtkCell> cell= pmesh->GetCell(0);
      if(VERBOSE) std::cout << "Cell 0 is of type: " << vtkCellTypes::GetClassNameFromTypeId(cell->GetCellType()) << std::endl;  
      if (cell->GetCellType() != VTK_POLY_LINE){ //from VTK/Filtering/vtkCellType.h, VTK_POLY_LINE: set of 1D lines
	std::cerr << "Cell 0 is not a VTK_POLY_LINE. Aborting!"<< std::endl; 
	////write test output 
	vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	twriter->SetFileName("test.vtp");
	twriter->SetInputData(pmesh);
	twriter->Write();
	std::cout   << "Wrote test.vtp" << std::endl;
	////write test output done.
	exit(1);
      }

      vtkSmartPointer<vtkPolyLine> pl= vtkPolyLine::SafeDownCast(cell);
      polygon->GetPointIds()->SetNumberOfIds(mnop);
      if(VERBOSE) std::cout << "Polygon will contain " << mnop << " points2."<< std::endl; 

      for (vtkIdType k= 0; k < mnop; k++){
	//polygon->GetPointIds()->InsertNextId(k);
	// outCell->GetPoints()->SetPoint(k, outPoints->GetPoint(k));//necessary???
	polygon->GetPointIds()->SetId(k, pl->GetPointIds()->GetId(k)); //works but may yield a "folded" polygon -> need convex hull
	//polygon->SetPointIds(pl->GetPointIds()); //works but may yield a "folded" polygon -> need convex hull
      }
 
      //  VTK_CREATE(vtkCellArray, outCells);
      outCells->InsertNextCell(polygon);
      pmesh->SetPolys(outCells);//SetPolys

      /////////////////create n-gon done.

      // //ct->SetInputConnection(dggp->GetOutputPort());
      tf->SetInputData(pmesh);
      /////////////////fill shortes-loop-path done.

      // ////write test output 
      // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      // twriter->SetFileName("test_mesh.vtp");
      // twriter->SetInputData(mesh);
      // twriter->Write();
      // std::cout   << "Wrote test_mesh.vtp" << std::endl;
      // ////write test output done.

      std::cout  << "Filled a nested loop." << std::endl;

      if(rmesh){
	//if(rmesh->GetPoints()->GetNumberOfPoints() > 0){
	VTK_CREATE(vtkAppendPolyData, append2);
      	append2->AddInputData(rmesh);
      	//if(sub_mesh->GetPoints())
      	//if(sub_mesh->GetPoints()->GetNumberOfPoints() > 0)
	append2->AddInputData(sub_mesh);

      	if(P_VERBOSE) 
	  std::cout << "Executing vtkAppendPolydata (append2)..." << std::endl;
	//append2->AddObserver(vtkCommand::ProgressEvent, progressCallback);
	append2->Update();
      	if(P_VERBOSE) 
	  std::cout  << std::endl << "done." << std::endl;

	//sub_mesh->DeepCopy(append2->GetOutput());
	sub_mesh= append2->GetOutput();
      }
      // else{
      // 	//sub_mesh= NULL;
      // 	sub_mesh= rmesh;
      // }
      //sub_mesh= rmesh;

    }/////////////////Fill one loop of a group of nested loops done.


    else{/////////////////Fill not nested loop 
  
      VTK_CREATE(vtkGraphToPolyData, graphToPolyData);
      graphToPolyData->SetInputConnection(cc_gth->GetOutputPort());
      
      if(P_VERBOSE) std::cout << "Executing vtkGraphToPolyData for not nested loop..." << std::endl;
      graphToPolyData->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      graphToPolyData->Update();
      if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

      // ////write test output 
      // vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      // twriter->SetFileName("test_ccgth.vtp");
      // //twriter->SetInputConnection(graphToPolyData->GetOutputPort());
      // twriter->Write();
      // std::cout   << "Wrote test_ccgth.vtp" << std::endl;
      // ////write test output done.


      /////////////////fill loop-path 

      /////////////////create n-gon

      // VTK_CREATE(vtkPolyData, pmesh);
      // VTK_CREATE(vtkPolygon, polygon);
      // VTK_CREATE(vtkCellArray, outCells);

      VTK_CREATE(vtkStripper, stripper);
      stripper->SetInputConnection(graphToPolyData->GetOutputPort());

      if(P_VERBOSE) std::cout << "Executing vtkStripper..." << std::endl;
      stripper->AddObserver(vtkCommand::ProgressEvent, progressCallback);
      stripper->Update();
      if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;


      //pmesh= dggp->GetOutput();
      //pmesh->DeepCopy(graphToPolyData->GetOutput());
      pmesh->DeepCopy(stripper->GetOutput());

      // VTK_CREATE(vtkPolygon, polygon);
      //polygon->GetPoints()->DeepCopy(mesh); //only good if n-gon is use without assignment to vtkPolyData: http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolygonIntersection
      // polygon->SetPoints(mesh->GetPoints());
      // polygon->SetPointIds(mesh->GetPointIds());
      vtkIdType mnop= pmesh->GetPoints()->GetNumberOfPoints();

      vtkSmartPointer<vtkCell> cell= pmesh->GetCell(0);
      if(VERBOSE) std::cout << "Cell 0 is of type: " << vtkCellTypes::GetClassNameFromTypeId(cell->GetCellType()) << std::endl;  
      if (cell->GetCellType() != VTK_POLY_LINE){ //from VTK/Filtering/vtkCellType.h, VTK_POLY_LINE: set of 1D lines
	std::cerr << "Cell 0 is not a VTK_POLY_LINE. Aborting!"<< std::endl; 
	////write test output 
	vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	twriter->SetFileName("testp.vtp");
	twriter->SetInputData(pmesh);
	twriter->Write();
	std::cout   << "Wrote test.vtp" << std::endl;
	////write test output done.
	exit(1);
      }

      vtkSmartPointer<vtkPolyLine> pl= vtkPolyLine::SafeDownCast(cell);
      polygon->GetPointIds()->SetNumberOfIds(mnop);
      if(VERBOSE) std::cout << "Polygon will contain " << mnop << " points2."<< std::endl; 

      for (vtkIdType k= 0; k < mnop; k++){
	//polygon->GetPointIds()->InsertNextId(k);
	// outCell->GetPoints()->SetPoint(k, outPoints->GetPoint(k));//necessary???
	polygon->GetPointIds()->SetId(k, pl->GetPointIds()->GetId(k)); //works but may yield a "folded" polygon -> need convex hull
	//polygon->SetPointIds(pl->GetPointIds()); //works but may yield a "folded" polygon -> need convex hull
      }
 
      // if(VERBOSE) std::cout << "Triangulating polygon...";
      // polygon->Triangulate(0,polygon->GetPointIds(),pmesh->GetPoints());
      // if(VERBOSE) std::cout << " done." << std::endl;

      //  VTK_CREATE(vtkCellArray, outCells);
      outCells->InsertNextCell(polygon);
      pmesh->SetPolys(outCells);//SetPolys

      /////////////////create n-gon done.

      // ct->SetInputConnection(graphToPolyData->GetOutputPort());
      tf->SetInputData(pmesh);
      /////////////////fill loop-path done.

      std::cout  << "Filled a not nested loop." << std::endl;
      //sub_mesh= NULL;

    }/////////////////Fill not nested loop done.
    

    // if(P_VERBOSE) std::cout << "Executing vtkContourTriangulator..." << std::endl;
    // ct->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    // ct->Update();
    // if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

    //append->AddInputConnection(ct->GetOutputPort());

    // //tf->PassLinesOff();//no lines will exist in ouput!	
    // if(P_VERBOSE) std::cout << "Executing vtkTriangleFilter..." << std::endl;
    // tf->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    // tf->Update();
    // if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

    // append->AddInputConnection(tf->GetOutputPort());
    append_fl->AddInputData(pmesh);

    if(P_VERBOSE) std::cout << "Executing vtkAppendPolydata..." << std::endl;
    //append_fl->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    append_fl->Update();
    if(P_VERBOSE) std::cout  << std::endl << "done." << std::endl;

    connected_filled_loops->DeepCopy(append_fl->GetOutput());

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

  // //if(P_VERBOSE) 
  // std::cout << "Executing vtkAppendPolydata..." << std::endl;
  // append->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  // append->Update();
  // //if(P_VERBOSE) 
  // std::cout  << std::endl << "done." << std::endl;

  // connected_filled_loops->DeepCopy(append->GetOutput());

  // //if(P_VERBOSE) 
  // std::cout << "Executing vtkAppendPolydata (append2)..." << std::endl;
  // append2->AddObserver(vtkCommand::ProgressEvent, progressCallback);
  // append2->Update();
  // //if(P_VERBOSE) 
  // std::cout  << std::endl << "done." << std::endl;
  // sub_mesh->DeepCopy(append2->GetOutput());

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
    //writer->SetCompression(true);
    writer->Write();
    std::cout << "Wrote " << argv[2] << std::endl;

     std::cout << "All done!" << std::endl;


    return EXIT_SUCCESS;
  }
