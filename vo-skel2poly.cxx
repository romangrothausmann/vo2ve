/////program to convert a 3D voxel-skeleton into a vtkPolyData (better to save than vtkGraph with points)
/////using Examples/Iterators/NeighborhoodIterators5.cxx and Examples/Iterators/NeighborhoodIterators6.cxx (for jumping)
//_01: iterating over whole image (i.e. not jumping adjacent fg-pixels), creating coincident points (i.e. lines are not connected to each other)


//#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNeighborhoodIterator.h>
//#include <itkConstNeighborhoodIteratorWithOnlyIndex.h> //needs itk >= 4.3.x

#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkMergePoints.h>
#include <vtkCellArray.h> //#include <vtkLines.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>

#include <vtkXMLPolyDataWriter.h>//for vtp-files 

#define VERBOSE 0

const unsigned int Dimension = 3;

typedef  unsigned char InputPixelType;
typedef  unsigned char OutputPixelType;

typedef itk::Image<InputPixelType,  Dimension>   InputImageType;
typedef itk::Image<OutputPixelType,  Dimension>   OutputImageType;

  
int main(int argc, char *argv[])
    {
    //InputImageType::Pointer image = InputImageType::New();

    if( argc != 3 )
        {
        std::cerr << "Usage: " << argv[0];
        std::cerr << " inputImage";
        std::cerr << " outputGraph";
        std::cerr << std::endl;  
        return EXIT_FAILURE;
        }

    if(!(strcasestr(argv[2],".vtp"))) {
        std::cerr << "The output should end with .vtp" << std::endl; 
        return -1;
        }

    typedef itk::ImageFileReader<InputImageType>  ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    std::string fileName = argv[1];
    reader->SetFileName(fileName);
    try
        {
        reader->Update();
        }
    catch ( itk::ExceptionObject &err)
        {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return -1;
        }

    InputImageType::Pointer image= reader->GetOutput();


    typedef itk::NeighborhoodIterator< InputImageType > NeighborhoodIteratorType;
    //typedef itk::ConstNeighborhoodIteratorWithOnlyIndex< InputImageType > NeighborhoodIteratorType;


    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType it(radius, image, image->GetRequestedRegion());


    //// specifying the root is not needed 
    // ImageType::IndexType root;
    // root[0] = 27;
    // root[1] = 29;
    // root[2] = 37;
    //// or just
    // ImageType::IndexType root = {{ 0, 0, 0 }};
    // it.SetLocation(root);


    typename InputImageType::IndexType point0ii, pointii;
    double point0[3], point[3];
    bool firstPoint= true;
    InputPixelType fg= 255;
    OutputPixelType doneValue= 1;
    vtkSmartPointer<vtkPoints> points= vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines= vtkSmartPointer<vtkCellArray>::New();

    std::cout << "Starting iteration!" << std::endl;
    //for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out)
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        {

        //std::cout << "Iteration at: " << it.GetIndex() << ": " << it.GetCenterNeighborhoodIndex() << std::endl;

        InputPixelType cpv= it.GetCenterPixel(); //ConstNeighborhoodIterator
        //InputPixelType cpv= it.GetCenterValue(); //ConstNeighborhoodIteratorWithOnlyIndex
        if (cpv == fg){
            if(VERBOSE) std::cout << "fg point found at: " << it.GetIndex() << ": " << it.GetCenterNeighborhoodIndex(); //it.GetIndex() == it.GetIndex(it.GetCenterNeighborhoodIndex())
            //point0ii= image->ComputeIndex(it.GetOffset(0));
            point0ii= it.GetIndex(it.GetCenterNeighborhoodIndex()); //it.GetIndex() == it.GetIndex(it.GetCenterNeighborhoodIndex())
            //point0ii= it.GetCenterNeighborhoodIndex();//ConstNeighborhoodIteratorWithOnlyIndex
            point0[0]= double(point0ii[0]);
            point0[1]= double(point0ii[1]);
            point0[2]= double(point0ii[2]);
            if(VERBOSE) std::cout << ": " << point0ii << std::endl;
            vtkIdType point0Id= points->InsertNextPoint(point0);
            
            it.SetPixel(it.GetCenterNeighborhoodIndex(), doneValue); //to keept track where we have worked already //not available in ConstNeighborhoodIteratorWithOnlyIndex
            //doneMap->SetPixel(it.GetOffset(0), doneValue);
            //for (unsigned i = it.Begin(); i != it.End(); i++)//not skipping centre pixel???
            for (unsigned i = 1; i < it.Size(); i++){ //skipping centre pixel???
                //if ((it.GetPixel(i) == fg) && (doneMap->GetPixel(it.GetOffset(i)) != doneValue)){
                if ((it.GetPixel(i) == fg) && (image->GetPixel(it.GetIndex(i)) != doneValue)){
                    //pointii= image->ComputeIndex(it.GetOffset(i));
                    pointii= it.GetIndex(i);
                    if(VERBOSE) std::cout << "Last point: " << point << std::endl;

                    //points->InsertNextPoint(double(point[0]), double(point[1]), double(point[2]));
                    point[0]= double(pointii[0]);
                    point[1]= double(pointii[1]);
                    point[2]= double(pointii[2]);
                    vtkIdType pointId= points->InsertNextPoint(point);

                    vtkSmartPointer<vtkLine> line= vtkSmartPointer<vtkLine>::New();
                    
                    line->GetPointIds()->SetId(0,point0Id);
                    line->GetPointIds()->SetId(1,pointId); //or GetNumberOf
                    lines->InsertNextCell(line);
                    }
                }
            }
        }

    vtkSmartPointer<vtkPolyData> polydata= vtkSmartPointer<vtkPolyData>::New();
    //polydata->SetPoints(reader->GetOutput()->GetPoints());
    polydata->SetPoints(points);
    polydata->SetLines(lines);

    vtkSmartPointer<vtkCleanPolyData> cleanFilter= vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    cleanFilter->SetInput(polydata);
#else
    cleanFilter->SetInputData(polydata);
#endif
    cleanFilter->PointMergingOn();
    cleanFilter->Update();


    //std::cout << "Input polydata contains " << reader->GetOutput()->GetNumberOfLines() << " lines and the loop-path " << polydata->GetNumberOfLines() << " lines."<< std::endl;

    vtkSmartPointer<vtkXMLPolyDataWriter> Pwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
 
    Pwriter->SetFileName(argv[2]);
    //Pwriter->SetInput(polydata);
    Pwriter->SetInputConnection(cleanFilter->GetOutputPort());
    Pwriter->Update();

    return EXIT_SUCCESS;
    }
