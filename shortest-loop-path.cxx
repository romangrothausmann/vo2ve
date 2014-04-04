 

/////program to find the shortes loop-path in a graph
////neede to remove triangular cc contained in the output of vo-skel2poly run through graph_ConnectedComponents_Boost



#include <vtkXMLPolyDataReader.h>//for vtp-files (cannot contain 3D cells?)

#include <vtkPolyDataToGraph.h>
#include <vtkDataSetAttributes.h>

#include <vtkVertexListIterator.h>
#include <vtkAdjacentVertexIterator.h>
#include <vtkInEdgeIterator.h>
#include <vtkOutEdgeIterator.h>
#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkMutableUndirectedGraph.h>

#include <vtkGraphToPolyData.h>

#include <vtkXMLPolyDataWriter.h>//for vtp-files 

#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


int searchVertex(vtkIdType vertex, vtkIdType target_vertex, vtkIdTypeArray slp_array, vtkGraph graph){

  VTK_CREATE(vtkAdjacentVertexIterator, it);
  graph->GetAdjacentVertices(vertex, it);
	
  if(!found){
    while (it->HasNext()){
      vtkIdType u = it->Next();
      if (u == target_vertex){
	found= 1;
      }
      else{
	searchVertex(u, graph);
      }
    }
    else{
      vtkIdType nv= slp->AddVertex();
      slp->AddEdge(lv, nv);
      std::cout << "Added vertex " << nv << " and an edge to the shortes-loop-path!"<< std::endl; 
  }
}


int main(int argc, char* argv[]){

    if( argc <= 3 ) {
        std::cerr << "Usage: " << argv[0];
        std::cerr << " inputGraph";
        std::cerr << " outputGraph";
        std::cerr << " startVertex";
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

    //vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New(); //*.vtp
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();  

    reader->SetFileName(argv[1]);
    reader->Update();

    ////cleaning should be done before this program is called to make sure the rootID does not change!!!


    vtkSmartPointer<vtkPolyDataToGraph> polyDataToGraphFilter= vtkSmartPointer<vtkPolyDataToGraph>::New();
    polyDataToGraphFilter->SetInputConnection(reader->GetOutputPort());
    polyDataToGraphFilter->Update();


    vtkGraph* graph= polyDataToGraphFilter->GetOutput();
    // vtkSmartPointer<vtkMutableUndirectedGraph> mdgraph=  vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    // mdgraph->DeepCopy(graph);

    //mdgraph->ShallowCopy(graph);
    // if (!graph->ShallowCopy(graph)){
    //     std::cerr << "1. MutableUndirectedGraph creation NOT successful! " << std::endl;
    //     return -1;
    //     }
 
//     vtkSmartPointer<vtkGraphToPolyData> tgraphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
// #if VTK_MAJOR_VERSION <= 5
//     tgraphToPolyData->SetInput(mdgraph);
// #else
//     tgraphToPolyData->SetInputData(mdgraph);
// #endif
//     tgraphToPolyData->Update();
      
//     vtkSmartPointer<vtkXMLPolyDataWriter> twriter= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//     twriter->SetFileName("test.vtp");
//     twriter->SetInputConnection(tgraphToPolyData->GetOutputPort());
//     twriter->Write();
        



    VTK_CREATE(vtkIdTypeArray, slp_array);


    int sv= atoi(argv[3]);
    vtkSmartPointer<vtkVertexListIterator> it= vtkSmartPointer<vtkVertexListIterator>::New();
    graph->GetVertices(it);

    bool ans;
    int n= 0;
    while(ans= it->HasNext()){ 
        n++;
        vtkIdType nextVertex= it->Next();
	//std::cout << "Checking vertex " << nextVertex << " with value: " << array->GetValue(nextVertex) << " which has " << graph->GetNumberOfEdgePoints(nextVertex) << " adjacent vertices" << std::endl;  

	vtkSmartPointer<vtkAdjacentVertexIterator> itc= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
	graph->GetAdjacentVertices(nextVertex, itc);
	
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
	  graph->GetAdjacentVertices(nextVertex, itc0);
	  while (itc0->HasNext()){
	    vtkSmartPointer<vtkAdjacentVertexIterator> itc1= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
	    graph->GetAdjacentVertices(itc0->Next(), itc1);
	    while (itc1->HasNext()){
	      vtkSmartPointer<vtkAdjacentVertexIterator> itc2= vtkSmartPointer<vtkAdjacentVertexIterator>::New(); //VTK_CREATE(vtkAdjacentVertexIterator, it);
	      graph->GetAdjacentVertices(itc1->Next(), itc2);
	      while (itc2->HasNext()){
		vtkIdType nextCVertex= itc2->Next();
		if (nextCVertex == nextVertex){
		  int vv= array->GetValue(nextCVertex);
		  std::cout << "Checking vertex " << nextCVertex << " with value: " << vv << std::endl;  
		  if (vv == com_int){
		    // VTK_CREATE(vtkInEdgeIterator, iei);
		    // graph->GetInEdges(nextCVertex, iei);
		    // while (iei->HasNext()){
		    //   vtkInEdgeType next_Edge= iei->Next();
		    //   graph->RemoveEdge(next_Edge.Id);
		    //   std::cout << "Removed in-edge " << next_Edge.Id << std::endl; 
		    //   }
		    // VTK_CREATE(vtkOutEdgeIterator, oei);
		    // graph->GetOutEdges(nextCVertex, oei);
		    // while (oei->HasNext()){
		    //   //vtkOutEdgeType next_Edge= ei->Next();
		    //   //graph->RemoveEdge(next_Edge);
		    //   graph->RemoveEdge(oei->Next().Id);
		    //   std::cout << "Removed out-edge " << oei->Next().Id << std::endl; 
		    //   }
		    // graph->RemoveVertex(nextCVertex); //Removes the vertex from the graph along with any connected edges. Note: This invalidates the last vertex index, which is reassigned to v. 
		    // std::cout << "Removed vertex " << nextCVertex << " with value: " << vv << " and its edges!"<< std::endl; 
		    // array->RemoveTuple(nextCVertex);
		    
		    int value_found= 0;
		    for(vtkIdType i= 0; i < vr_array->GetNumberOfTuples(); i++)
		      if (vr_array->GetValue(i) == nextCVertex)
			value_found= 1;
		    if(!value_found){
		      vr_array->InsertNextValue(nextCVertex);
		      std::cout << "Added vertex " << nextCVertex << " with value: " << vv << " to the vertex-remove-list!"<< std::endl; 
		      VTK_CREATE(vtkInEdgeIterator, iei);
		      graph->GetInEdges(nextCVertex, iei);
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

    graph->RemoveEdges(er_array); //Removes a collection of vertices from t
    std::cout << "Removed edges!"<< std::endl; 
    // graph->RemoveVertices(vr_array); //Removes a collection of vertices from the graph along with any connected edges. 
    // std::cout << "Removed vertices and their edges!"<< std::endl; 


    vtkSmartPointer<vtkGraphToPolyData> graphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    graphToPolyData->SetInput(graph);
#else
    graphToPolyData->SetInputData(graph);
#endif
    graphToPolyData->Update();
      
    vtkSmartPointer<vtkXMLPolyDataWriter> writer= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(argv[2]);
    writer->SetInputConnection(graphToPolyData->GetOutputPort());
    writer->Write();
        


    return EXIT_SUCCESS;
    }
