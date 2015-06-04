/////program to convert a 3D voxel-skeleton into a vtkPolyData (better to save than vtkGraph with points)
/////using Examples/Iterators/NeighborhoodIterators5.cxx and Examples/Iterators/NeighborhoodIterators6.cxx (for jumping)
//01: iterating over whole image (i.e. not jumping adjacent fg-pixels), creating coincident points (i.e. lines are not connected to each other)
//02: points then merged by vtkCleanPolyData
//03: clean up

//ToDo:
//If a path of length pathLength already exists, an additional line does not necessarily create a loop of length pathLength+1 !!!! see e.g. extract_01_skel-ana_01_vskel_pruned-ends_ana_x136y78z47_p+2.mha
//- code needs to be modified to actually check for the cylcle length AFTER the new line as been added and if it creates a loop below maxLoopLengh remove the line again
//- only run iterator with half of the radius filled with 1, track keeping not needed than any more?

#include <itkImageFileReader.h>
#include <itkNeighborhoodIterator.h>

#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkMergePoints.h>
#include <vtkXMLPolyDataWriter.h>

#define VERBOSE 0

const unsigned int Dimension = 3;

typedef  unsigned char InputPixelType;

typedef itk::Image<InputPixelType,  Dimension>   InputImageType;


void i2dA(InputImageType::IndexType iArray, double dArray[]){
    for (int i=0; i<Dimension; ++i) dArray[i] = (double)iArray[i];
}

//// could be implemented with vtkDijkstraGraphGeodesicPath returning loop length (shortest path + 1) only if below maxLoopLength
//// however as the search only needs to go up to a loop length of maxLoopLength, vtkDijkstraGraphGeodesicPath is not used
unsigned char pathExists(const vtkIdType pId, const vtkIdType pIdEnd, vtkPolyData* mesh, unsigned char pathLength, unsigned char maxPathLength){//does this need a char overflow check?

    if(VERBOSE) std::cerr << std::setw (pathLength+1) << +pathLength << ": " << std::endl;
    if(VERBOSE) std::cerr << std::setw (pathLength+1) << +pathLength << ": " << "pId " << +pId << std::endl;
    if(VERBOSE) std::cerr << std::setw (pathLength+1) << +pathLength << ": " << "pIdEnd " << +pIdEnd << std::endl;
    //get all cells that vertex 'id' is a part of
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    //mesh->BuildLinks();//ought to be called before pathExists as mesh is not changed during recursion
    mesh->GetPointCells(pId, cellIdList); //make sure that BuildLinks() has been called.
    if(VERBOSE) std::cerr << std::setw (pathLength+1) << +pathLength << ": " << "GetNumberOfIds " << +cellIdList->GetNumberOfIds() << std::endl;

    //loop through all the cells that use the seed point
    for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++){//so far does not avoid going back to previous visited vertices
        //if(cell->GetNumberOfEdges() <= 0)//check if cell is a line; can be skipped as pd is constructed in this program

        if(VERBOSE) std::cerr << std::setw (pathLength+1) << +pathLength << ": " << "i " << +i << " GetNumberOfIds " << +cellIdList->GetNumberOfIds() << std::endl;

        vtkLine* line = vtkLine::SafeDownCast(mesh->GetCell(cellIdList->GetId(i)));
        if(line){
            vtkIdType p0, p1, pIdNext;
            p0 = line->GetPointId(0);
            p1 = line->GetPointId(1);

            if(p0 == pId)
                pIdNext= p1;
            else
                pIdNext= p0;
            //pathLength++;
            //vtkIdList: pathPointIds->InsertNextId(pIdNext);//save walked path; not needed here
            if(VERBOSE) std::cerr << std::setw (pathLength+1) << +pathLength << ": " << "pIdNext " << +pIdNext << std::endl;

            if(pIdNext == pIdEnd){//end-vertex reached?
                if(VERBOSE) std::cerr << std::setw (pathLength+1) << +pathLength << "pIdEnd reached! pathLength: " << +pathLength << std::endl;
                return(pathLength);
                }
            else if(pathLength < maxPathLength){//maxPathLength reached?
                unsigned char length= pathExists(pIdNext, pIdEnd, mesh, pathLength+1, maxPathLength);
                if(length)
                    return(length);
                }
            // else return(0); //do not return untill all branches are checked!
            }
        else
            std::cerr << "celltype: " << vtkCellTypes::GetClassNameFromTypeId(mesh->GetCell(cellIdList->GetId(i))->GetCellType()) << std::endl;
        }
    return(0);
    }

int main(int argc, char *argv[]){

    if(argc != 4){
        std::cerr << "Usage: " << argv[0];
        std::cerr << " inputImage";
        std::cerr << " outputGraph";
        std::cerr << " maxPathLength";
        std::cerr << std::endl;
        return EXIT_FAILURE;
        }

    if(!(strcasestr(argv[2],".vtp"))){
        std::cerr << "The output should end with .vtp" << std::endl;
        return EXIT_FAILURE;
        }

    const unsigned char maxPathLength= atoi(argv[3]); //e.g. 2 to only avoid triangular loops (and double lines)

    typedef itk::ImageFileReader<InputImageType>  ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    std::string fileName = argv[1];
    reader->SetFileName(fileName);
    try {
        reader->Update();
        }
    catch ( itk::ExceptionObject &err){
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
        }

    InputImageType::Pointer image= reader->GetOutput();
    InputImageType::RegionType region= image->GetLargestPossibleRegion();
    InputImageType::SizeType size= region.GetSize();

    double bounds[6];
    for(int i= 0; i < Dimension; i++){
        bounds[i]= region.GetIndex()[i];
        bounds[i+1]= region.GetSize()[i];
        }

    typedef itk::NeighborhoodIterator<InputImageType> NeighborhoodIteratorType;
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1); //26-connectivity
    NeighborhoodIteratorType it(radius, image, region);

    InputPixelType fg= 255;
    InputPixelType doneValue= 1;

    vtkSmartPointer<vtkPolyData> polydata= vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points= vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines= vtkSmartPointer<vtkCellArray>::New();

    polydata->SetPoints(points);
    polydata->SetLines(lines);

    vtkSmartPointer<vtkMergePoints> mergePoints= vtkSmartPointer<vtkMergePoints>::New();
    mergePoints->SetDataSet(polydata);
    mergePoints->SetDivisions(10,10,10);
    mergePoints->InitPointInsertion(points, bounds);

    vtkIdType point0Id, pointId;
    double point0[3], point[3];

    std::cout << "Starting iteration!" << std::endl;
    InputImageType::RegionType::SizeValueType counter= 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it){

        InputPixelType cpv= it.GetCenterPixel();
        if (cpv == fg){
            if(VERBOSE) std::cout << "fg point found at: " << it.GetIndex();
            i2dA(it.GetIndex(),point0);//it.GetIndex() == it.GetIndex(it.GetCenterNeighborhoodIndex())
            if(VERBOSE) std::cout << ": " << point0[0] << "; " << point0[1] << "; " << point0[2] << std::endl;

            if(!mergePoints->InsertUniquePoint(point0, point0Id))
                if(VERBOSE) std::cout << "Primary point already existed: " << point0[0] << "; " << point0[1] << "; " << point0[2] <<  std::endl;

            it.SetCenterPixel(doneValue); //same as it.SetPixel(it.GetCenterNeighborhoodIndex(), doneValue); //to keept track where we have worked already
            for (unsigned i= 0; i < it.Size(); i++){ //not skipping center pixel at it.GetCenterNeighborhoodIndex() (==13 for 3x3x3)!!! looping see e.g.: http://www.itk.org/Doxygen47/html/Examples_2Iterators_2NeighborhoodIterators6_8cxx-example.html
		//if(i == it.GetCenterNeighborhoodIndex()) continue; //skips center pixel, not needed with it.SetCenterPixel(doneValue)
                if ((it.GetPixel(i) == fg)){//skips center pixel with it.SetCenterPixel(doneValue)

                    if(VERBOSE) std::cout << "GetIndex" << std::endl;
                    i2dA(it.GetIndex(i), point);
                    if(VERBOSE) std::cout << ": " << point[0] << "; " << point[1] << "; " << point[2] << std::endl;

                    if(!mergePoints->InsertUniquePoint(point, pointId)){//only existing points can cause loops
                        if(VERBOSE) std::cout << "Point already existed: " << point[0] << "; " << point[1] << "; " << point[2] << std::endl;//should never happen when only half the radius is filled
                        // polydata->SetPoints(points);
                        // polydata->SetLines(lines);
                        polydata->DeleteCells();//must be called before a successive call to BuildLinks
                        polydata->BuildLinks();//has to be called before pathExists
                        if(VERBOSE) std::cout << "BuildLinks" << std::endl;
                        unsigned char pathLength= pathExists(point0Id, pointId, polydata, 0, maxPathLength);
                        //std::cout << "Line would create a loop of length: " << +pathLength << std::endl;
                        if(pathLength){
                            //if(VERBOSE)
                            std::cout << "A path of length " << +pathLength << " already exists!" << std::endl; //does not necessarily create a loop of length pathLength+1 !!!!
                            continue;
                            }
                        }

                    vtkSmartPointer<vtkLine> line= vtkSmartPointer<vtkLine>::New();

                    line->GetPointIds()->SetId(0,point0Id);
                    line->GetPointIds()->SetId(1,pointId);
                    lines->InsertNextCell(line);


                    }
                }
            }

        counter++;
        fprintf(stdout, "\rprogress: %6.2f%%", counter*100.0/region.GetNumberOfPixels());

        }
    std::cout << "Iteration done." << std::endl;

    ////should Squeeze + SetPoints + Delete only be done in a filter as in vtkCleanPolyData.cxx?
    // points->Squeeze();//compress dynamic array as no further extension follows
    // polydata->SetPoints(points);//already done above, needed here again?
    // points->Delete();//decrement reference count

    // lines->Squeeze();//compress dynamic array as no further extension follows
    // polydata->SetLines(lines);//already done above, needed here again?
    // lines->Delete();//decrement reference count
    // std::cout << "Squeezing done." << std::endl;

    vtkSmartPointer<vtkXMLPolyDataWriter> Pwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    Pwriter->SetFileName(argv[2]);
    Pwriter->SetInputData(polydata);//will segfault if Delete was called on points and lines before
    Pwriter->Update();

    std::cerr << "# verts: " << polydata->GetNumberOfPoints() << std::endl;
    std::cerr << "# lines: " << polydata->GetNumberOfCells() << std::endl;
    std::cout << "EPC(26)= " << polydata->GetNumberOfPoints() - polydata->GetNumberOfCells() << std::endl;

    return EXIT_SUCCESS;
    }
