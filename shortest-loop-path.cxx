 

/////program to find the shortes loop-path sub-graph in a graph (vtkPolyData) using vtkDijkstraGraphGeodesicPath


#include <vtkSmartPointer.h>


#include <vtkXMLPolyDataReader.h>//for vtp-files (cannot contain 3D cells?)

#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkExtractEdges.h>
#include <vtkIdList.h>

#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkExtractSelectedGraph.h>

//#include <.h>

#include <vtkXMLPolyDataWriter.h>//for vtp-files 

#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


////from http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/VertexConnectivity
//vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);

 
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


    vtkIdType sVertex= atoi(argv[3]);
    vtkIdType eVertex;

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
 
    //vtkSmartPointer<vtkPolyData> mesh = extractEdges->GetOutput();
    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
   
    std::cout << "Prepared mesh!" << std::endl;

    vtkIdType id= sVertex;

    ////get connected vertices and the cells connecting them
    vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();

    //get all cells that vertex 'id' is a part of
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    mesh->BuildLinks();
    mesh->GetPointCells(id, cellIdList); //Make sure that routine BuildLinks() has been called. 
 
    if (cellIdList->GetNumberOfIds() != 2){
      std::cerr << "Point not connected to just two lines (cells)! # of connected cells: "<< cellIdList->GetNumberOfIds() << "loop splitting ambiguous!"<< std::endl;
      return EXIT_FAILURE;
    }
    /*
    cout << "Vertex 0 is used in cells ";
    for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
      {
	cout << cellIdList->GetId(i) << ", ";
      }
    cout << endl;
    */
 
    for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
      {
	//cout << "id " << i << " : " << cellIdList->GetId(i) << endl;
 
	vtkSmartPointer<vtkIdList> pointIdList =  vtkSmartPointer<vtkIdList>::New();
	//mesh-> BuildCells(); //not needed
	mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);
 
	cout << "End points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;
 
	if(pointIdList->GetId(0) != id)
	  {
	    //cout << "Connected to " << pointIdList->GetId(0) << endl;
	    connectedVertices->InsertNextId(pointIdList->GetId(0));
	  }
	else
	  {
	    //cout << "Connected to " << pointIdList->GetId(1) << endl;
	    connectedVertices->InsertNextId(pointIdList->GetId(1));
	  }
      }


    if (connectedVertices->GetNumberOfIds() > 2){
      std::cerr << "More than just one connected point (" << connectedVertices->GetNumberOfIds() << "), endpoint ambiguous!"<< std::endl;
      return EXIT_FAILURE;
    }

    eVertex= connectedVertices->GetId(0);
    std::cout << "Endpoint unambiguous! Using: "<< eVertex << std::endl;
    
    ////split open line connection between sVertex and eVertex
    mesh->DeleteCell(cellIdList->GetId(0)); //has to correspond to eVertex= connectedVertices->GetId(0);
    //mesh->RemoveDeletedCells(); 	
    std::cout << "Deleted cell: "<< cellIdList->GetId(0) << std::endl;


    VTK_CREATE(vtkDijkstraGraphGeodesicPath, dggp);

    //dggp->SetInputConnection(polyDataToGraphFilter->GetOutputPort());
    dggp->SetInputConnection(reader->GetOutputPort());
    //dggp->SetInputConnection(extractEdges->GetOutputPort());
    dggp->SetStartVertex(sVertex);
    dggp->SetEndVertex(eVertex);
// StopWhenEndReachedOn ()
// UseScalarWeightsOn ()
// RepelPathFromVerticesOn ()
// SetRepelVertices (vtkPoints *)

      
//     VTK_CREATE(vtkExtractSelectedGraph, esg);
//     esg->AddInputConnection(polyDataToGraphFilter->GetOutputPort());
//     esg->SetSelectionConnection(dggp->GetOutputPort());

//     vtkSmartPointer<vtkGraphToPolyData> graphToPolyData= vtkSmartPointer<vtkGraphToPolyData>::New();
// // #if VTK_MAJOR_VERSION <= 5
// //     graphToPolyData->SetInput(graph);
// // #else
// //     graphToPolyData->SetInputData(graph);
// // #endif
// //     graphToPolyData->Update();
//     graphToPolyData->SetInputConnection(esg->GetOutputPort());

    vtkSmartPointer<vtkXMLPolyDataWriter> writer= vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(argv[2]);
    //writer->SetInputConnection(graphToPolyData->GetOutputPort());
    writer->SetInputConnection(dggp->GetOutputPort());
    writer->Write();
        


    return EXIT_SUCCESS;
    }
