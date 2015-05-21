/////convert a vector (vtkPolyData) skeleton to a graph (vtkGraph)


#include <vtkXMLPolyDataReader.h>//for vtp-files (cannot contain 3D cells?)
//#include <vtkPolyDataReader.h>//for vtk-files (cannot contain 3D cells?)


#include <vtkCircularLayoutStrategy.h>
#include <vtkGraphLayoutView.h>
//#include <vtkIntArray.h>
//#include <vtkMutableUndirectedGraph.h>
#include <vtkRenderWindowInteractor.h>


#include "vtkPolyDataToGraph.h"

int main(int argc, char* argv[]){

  if( argc <= 1 ) {
      std::cerr << "Usage: " << argv[0];
      std::cerr << " inputGraph";
      //std::cerr << " outputPoints";
      std::cerr << std::endl;  
      return EXIT_FAILURE;
    }

  if(!(strcasestr(argv[1],".vtp"))) {
    std::cout << "The input should end with .vtp" << std::endl; 
    return -1;
  }

  //vtkPolyDataReader *reader = vtkPolyDataReader::New(); //*.vtk
  vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New(); //*.vtp
 
  reader->SetFileName(argv[1]);
  reader->Update();


  vtkSmartPointer<vtkPolyDataToGraph> polyDataToGraphFilter= vtkSmartPointer<vtkPolyDataToGraph>::New();
  polyDataToGraphFilter->SetInputConnection(reader->GetOutputPort());
  polyDataToGraphFilter->Update();

vtkSmartPointer<vtkCircularLayoutStrategy> circularLayoutStrategy =  vtkSmartPointer<vtkCircularLayoutStrategy>::New();

  vtkSmartPointer<vtkGraphLayoutView> graphLayoutView = vtkSmartPointer<vtkGraphLayoutView>::New();
  graphLayoutView->AddRepresentationFromInput(polyDataToGraphFilter->GetOutput());
  //graphLayoutView->AddRepresentationFromInput(graphToPolyDataFilter->GetOutputPort());

  //graphLayoutView->SetLayoutStrategy(circularLayoutStrategy);
  //graphLayoutView->SetLayoutStrategy("Simple 2D");

  //graphLayoutView->SetVertexLabelVisibility(true);
  graphLayoutView->SetEdgeLabelVisibility(true);
  graphLayoutView->SetEdgeVisibility(true);
  //graphLayoutView->SetEdgeLabelArrayName("Weights"); //default is "labels"
  //graphLayoutView->SetVertexLabelArrayName("VertexIDs"); //default is "labels"
  graphLayoutView->ResetCamera();
  graphLayoutView->Render();
  graphLayoutView->GetInteractor()->Start();
 
  return EXIT_SUCCESS;
}
