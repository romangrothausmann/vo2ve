/////program to convert a 3D voxel-skeleton into a vtkPolyData (better to save than vtkGraph with points)
/////using Examples/Iterators/NeighborhoodIterators5.cxx and Examples/Iterators/NeighborhoodIterators6.cxx (for jumping)
//01: iterating over whole image (i.e. not jumping adjacent fg-pixels), creating coincident points (i.e. lines are not connected to each other)
//02: points then merged by vtkCleanPolyData
//03: clean up


#include <itkImageFileReader.h>
#include <itkNeighborhoodIterator.h>

#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#define VERBOSE 0

const unsigned int Dimension = 3;

typedef  unsigned char InputPixelType;

typedef itk::Image<InputPixelType,  Dimension>   InputImageType;


int main(int argc, char *argv[]){

    if(argc != 3){
        std::cerr << "Usage: " << argv[0];
        std::cerr << " inputImage";
        std::cerr << " outputGraph";
        std::cerr << std::endl;
        return EXIT_FAILURE;
        }

    if(!(strcasestr(argv[2],".vtp"))){
        std::cerr << "The output should end with .vtp" << std::endl;
        return EXIT_FAILURE;
        }

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

    typedef itk::NeighborhoodIterator<InputImageType> NeighborhoodIteratorType;
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType it(radius, image, image->GetRequestedRegion());

    typename InputImageType::IndexType point0ii, pointii;
    double point0[3], point[3];
    bool firstPoint= true;
    InputPixelType fg= 255;
    InputPixelType doneValue= 1;
    vtkSmartPointer<vtkPoints> points= vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines= vtkSmartPointer<vtkCellArray>::New();

    std::cout << "Starting iteration!" << std::endl;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it){

        InputPixelType cpv= it.GetCenterPixel();
        if (cpv == fg){
            if(VERBOSE) std::cout << "fg point found at: " << it.GetIndex() << ": " << it.GetCenterNeighborhoodIndex();

            point0ii= it.GetIndex(it.GetCenterNeighborhoodIndex());
            point0[0]= double(point0ii[0]);
            point0[1]= double(point0ii[1]);
            point0[2]= double(point0ii[2]);
            if(VERBOSE) std::cout << ": " << point0ii << std::endl;
            vtkIdType point0Id= points->InsertNextPoint(point0);

            it.SetPixel(it.GetCenterNeighborhoodIndex(), doneValue); //to keept track where we have worked already
            for (unsigned i = 1; i < it.Size(); i++){ //skipping centre pixel???
                if ((it.GetPixel(i) == fg) && (image->GetPixel(it.GetIndex(i)) != doneValue)){
                    pointii= it.GetIndex(i);
                    if(VERBOSE) std::cout << "Last point: " << point << std::endl;

                    point[0]= double(pointii[0]);
                    point[1]= double(pointii[1]);
                    point[2]= double(pointii[2]);
                    vtkIdType pointId= points->InsertNextPoint(point);

                    vtkSmartPointer<vtkLine> line= vtkSmartPointer<vtkLine>::New();

                    line->GetPointIds()->SetId(0,point0Id);
                    line->GetPointIds()->SetId(1,pointId);
                    lines->InsertNextCell(line);
                    }
                }
            }
        }

    vtkSmartPointer<vtkPolyData> polydata= vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetLines(lines);

    vtkSmartPointer<vtkCleanPolyData> cleanFilter= vtkSmartPointer<vtkCleanPolyData>::New();
    cleanFilter->SetInputData(polydata);
    cleanFilter->PointMergingOn();
    cleanFilter->Update();

    vtkSmartPointer<vtkXMLPolyDataWriter> Pwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    Pwriter->SetFileName(argv[2]);
    Pwriter->SetInputConnection(cleanFilter->GetOutputPort());
    Pwriter->Update();

    return EXIT_SUCCESS;
    }
