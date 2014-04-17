#ifndef _BME_utility_h
#define _BME_utility_h

#pragma once
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkGenericCell.h>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_matlab_filewrite.h>
//#include <vnl_matlab_fileread.h>

#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkMesh.h>
#include <itkAutomaticTopologyMeshSource.h>
#include <itkMeshSource.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include "itkCastImageFilter.h"
#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkAddImageFilter.h>
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include <itkBinaryThresholdImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkImageDuplicator.h>

#include <vtkLinearExtrusionFilter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkStripper.h>
#include <vtkImageReslice.h>
#include <vtkCutter.h>
#include <vtkImageAccumulate.h>
#include <vtkTransform.h>
#include <vtkMetaImageReader.h>
#include "vtkMetaImageWriter.h"
#include <vtkPolyDataReader.h>
#include <vtkCellArray.h>


template <class T,int Dim> class BME_utility{
private:
public:
	int dim;
	enum cavityNR
	{
		LV,
		RV,
		LA,
		RA,
		Ao,
	}cavity;
	typedef float PixelType;
	typedef itk::Image<PixelType, 3>			ImageType;

	typedef unsigned char						CharPixelType;
	typedef itk::Image<CharPixelType, 3>		CharImageType3D;
	typedef itk::Image<CharPixelType, 2>		CharImageType2D;
	typedef CharImageType3D::Pointer			CharImagePointerType3D;
	typedef CharImageType2D::Pointer			CharImagePointerType2D;

	typedef signed int							signedCharPixelType;
	typedef itk::Image<signedCharPixelType, 2>	INT8ImageType2D;
	typedef itk::Image<signedCharPixelType, 3>	INT8ImageType3D;
	typedef INT8ImageType3D::Pointer		signedCharImagePointerType3D;
	typedef INT8ImageType2D::Pointer		signedCharImagePointerType2D;


	typedef T									BME_PixelType;
	typedef itk::Image<float,Dim>				BME_FLOATImageType;
	typedef itk::Image<signedCharPixelType,Dim>	BME_INT8ImageType;
	typedef itk::Image<T,Dim>					BME_ImageType;

	typedef itk::Image<T,2>						BME_ImageType2D;
	typedef itk::Image<T,3>						BME_ImageType3D;

	typedef itk::CastImageFilter< BME_ImageType, BME_ImageType3D> BME_CastType3D;
	typedef itk::AndImageFilter< BME_ImageType, BME_ImageType, BME_ImageType > BME_AndImageFiltertype;
	typedef itk::OrImageFilter<  BME_ImageType, BME_ImageType, BME_ImageType > BME_OrImageFiltertype;
	typedef itk::AddImageFilter< BME_ImageType, BME_ImageType, BME_ImageType > BME_AddImageFiltertype;

	typedef itk::ExtractImageFilter<BME_ImageType, BME_ImageType2D> BME_ImageFilterNDto2D; //Get Slice out of 3D stack
	typedef itk::ImageDuplicator<BME_ImageType> BME_DuplicatorType;  
	typedef itk::ImageDuplicator<BME_ImageType2D> BME_DuplicatorType2D;  
	typedef itk::ImageDuplicator<BME_ImageType3D> BME_DuplicatorType3D;  
	typedef itk::ImageDuplicator<BME_INT8ImageType> BME_signedDuplicatorType;  

	typedef itk::BinaryThresholdImageFilter<BME_INT8ImageType,BME_ImageType> BME_BinaryThreshold;
	
	typedef itk::SignedDanielssonDistanceMapImageFilter<BME_ImageType,BME_ImageType>  BME_DistanceMapImageFilter;

	typedef itk::ImageToVTKImageFilter<BME_ImageType> BME_ImageConverterITKToVTK;
	typedef itk::VTKImageToImageFilter<BME_ImageType>BME_ImageConverterVTKToITK;

public:
	BME_utility();
	~BME_utility();

	bool LoadVtkObject( const std::string inFileNamePath, short PlaneNR,short cavityNR,vtkSmartPointer<vtkPolyData> mesh);
	void LoadTransform( const std::string inFileNamePath, vnl_matrix_fixed<double, 4, 4>  landTransPlanesToWorld ,vnl_matrix_fixed<double, 4, 4>  landTransModelToWorld);
	void LoadPlaneTransform( const std::string inFileNamePath,short PlaneNR, vnl_matrix_fixed<double, 4, 4> &vnllandTransmatrix,double *origin , double *normal);
	
	void Get2DImageFrom2DPolydata( const vtkSmartPointer<vtkPolyData> mesh,vtkSmartPointer <vtkImageData> MaskImage, vnl_matrix_fixed<double, 4, 4> &vnllandTransmatrix,double *origin , double *normal, BME_DuplicatorType *outImage, BME_DuplicatorType *outMask);
	void Get2DImageFrom3DPolydata( const vtkSmartPointer<vtkPolyData> mesh,vtkSmartPointer <vtkImageData> MaskImage, vnl_matrix_fixed<double, 4, 4> &vnllandTransmatrix,double *origin , double *normal, BME_DuplicatorType *outImage, BME_DuplicatorType *outMask);

	vtkSmartPointer<vtkPolyData> computePolydata( const vtkSmartPointer<vtkPolyData> mesh, double *origin , double *normal);
	vtkSmartPointer<vtkPolyData> ConvertImageToPolydata(BME_ImageType2D * inImage,double *origin );
	
	void TransformPolydataTo2D( vtkSmartPointer<vtkPolyData> input, vtkSmartPointer<vtkPolyData> output,vnl_matrix_fixed<double, 4, 4> &vnlTransmatrix);
	
	void SortPolydata(vtkSmartPointer<vtkPolyData> Inpolydata , vtkSmartPointer<vtkPoints> OutpolydataPoints);
	
	vtkSmartPointer<vtkImageStencil>  ConvertPolydataToStencil(double *spacing, double *bounds,vtkSmartPointer<vtkPolyData> polyData);

	void ConvertPolydataToITKStencil(double *spacing, double *bounds,vtkSmartPointer<vtkPolyData> polyData,BME_DuplicatorType2D *outImage2D,BME_DuplicatorType3D *outImage3D);
	double GetVoxelCount(BME_ImageType *Image);
	
	void GetIntersectionImage(BME_ImageType *inImage1,BME_ImageType *inImage2,BME_DuplicatorType *outImage);
	void GetUnionImage(BME_ImageType *inImage1,BME_ImageType *inImage2,BME_DuplicatorType *outImage);
	void GetMaskedImage(BME_ImageType *inImage1,BME_ImageType *inImage2,BME_DuplicatorType *outImage);

	void GetThresholdImage(BME_INT8ImageType *inImage,BME_DuplicatorType *outImage);
	double * GetPlaneImage(vnl_matrix_fixed<double, 4, 4> vnllandTransmatrix,double * origin,vtkSmartPointer <vtkImageData> InputImage,BME_DuplicatorType *outImage);

	void GetDistanceMap(BME_ImageType *inImage,BME_signedDuplicatorType *outImage);

	vnl_vector<T> vnlShapeMatrixToVector( const vnl_matrix<T> matrix );
	void WriteMHD(const std::string inFileNamePath,BME_ImageType *inImage);
#if 0 //Not OK Yet 
	vtkSmartPointer<vtkImageData> ImageConverterITKToVTK2D( const BME_ImageType2D InputImage );
	vtkSmartPointer<vtkImageData> ImageConverterITKToVTK3D( const BME_ImageType3D InputImage );
	BME_ImageType2D ImageConverterVTKToITK2D( vtkSmartPointer <vtkImageData> InputImage );
	BME_ImageType3D ImageConverterVTKToITK3D(  vtkSmartPointer <vtkImageData> InputImage );
#endif
	vnl_matrix<T> vnlShapeVectorToMatrix( const vnl_vector<T> vector );
	vnl_matrix<T> vtkPointsToVnlPointMatrix( const vtkSmartPointer<vtkPoints> inPts );
	vtkSmartPointer<vtkPoints> vnlPointMatrixToVtkPoints( const vnl_matrix<T> inPts );
	vtkSmartPointer<vtkPolyData> vnlPointMatrixToVtkMesh( const vnl_matrix<T> S, vtkSmartPointer<vtkPolyData> templateMesh );
	vnl_matrix<T> vtkMeshToVnlPointMatrix( vtkSmartPointer<vtkPolyData> poly );
	vnl_vector<T> computeOverlap( const vtkSmartPointer<vtkPolyData> mesh, const CharImagePointerType3D image );
	vnl_vector<T> computeVolume( const vtkSmartPointer<vtkPolyData> mesh, const CharImagePointerType3D image );
	void vtkMeshToItkImage( const vtkSmartPointer<vtkPolyData> polyData, CharImagePointerType3D  & outImage );
	void sortVnlVector( const vnl_vector<T> inVector, vnl_vector<T> & sortedVector, vnl_vector<int> & idxVector );
	vtkSmartPointer<vtkPolyData> getSpecificMesh( const int idx, const vtkSmartPointer<vtkPolyData> inMesh );
	vnl_matrix<T> transformPoints( vnl_matrix_fixed<T, 4, 4> B, vnl_matrix<T> points );

	void matlabSaveVnlVector( const std::string outFileName, vnl_vector<T> outVector, const std::string outVariableName );
	void matlabloadVnlVector( const std::string inFileName, vnl_vector<T> inVector, const std::string inVariableName );
	void matlabSaveVnlMatrix( const std::string outFileName, vnl_matrix<T> outMatrix, const std::string outVariableName );
	void matlabloadVnlMatrix( const std::string inFileName, vnl_matrix<T> inMatrix, const std::string inVariableName );

};

template <class T,int Dim> void BME_utility<T,Dim>:: Get2DImageFrom2DPolydata( const vtkSmartPointer<vtkPolyData> mesh,vtkSmartPointer <vtkImageData> MaskImage, vnl_matrix_fixed<double, 4, 4> &vnlPlaneTransmatrix,double *origin , double *normal, BME_DuplicatorType *outImage, BME_DuplicatorType *outMask)

{
	double *imageBounds;
	double imageSpacing[4];
	//imageSpacing[0]=imageSpacing[1]=imageSpacing[2]=1; // output of ImageReslice , set to 1
	//imageSpacing[3]=0;

	BME_utility<CharPixelType,2>::BME_DuplicatorType2D::Pointer StencilImage2D = BME_utility<CharPixelType,2>::BME_DuplicatorType2D::New();;
	BME_utility<CharPixelType,3>::BME_DuplicatorType3D::Pointer StencilImage3D = BME_utility<CharPixelType,3>::BME_DuplicatorType3D::New();;
	
	BME_utility<CharPixelType,2>::BME_DuplicatorType::Pointer MaskImage2D = BME_utility<CharPixelType,2>::BME_DuplicatorType::New();;
	
		
	imageBounds=GetPlaneImage(vnlPlaneTransmatrix, origin,MaskImage,MaskImage2D);
  
  BME_ImageType2D::SpacingType Spacing = MaskImage2D->GetOutput()->GetSpacing();
  imageSpacing[0]=Spacing[0];
  imageSpacing[1]=Spacing[1];
	
  vtkSmartPointer<vtkPolyData> meshSortedPolyData2D= vtkSmartPointer<vtkPolyData>::New();	
	
	TransformPolydataTo2D(mesh, meshSortedPolyData2D,vnlPlaneTransmatrix);
	
	BME_utility<CharPixelType,2>::BME_ImageType::Pointer Image2DmeshSorted = BME_utility<CharPixelType,2>::BME_ImageType::New();;
	
	ConvertPolydataToITKStencil(imageSpacing,imageBounds,meshSortedPolyData2D,StencilImage2D,StencilImage3D);

	outImage->SetInputImage(StencilImage2D->GetOutput());
	outImage->Update();
	outMask->SetInputImage(MaskImage2D->GetOutput());
	outMask->Update();

}


template <class T,int Dim> void BME_utility<T,Dim>:: Get2DImageFrom3DPolydata( const vtkSmartPointer<vtkPolyData> mesh,vtkSmartPointer <vtkImageData> MaskImage, vnl_matrix_fixed<double, 4, 4> &vnlPlaneTransmatrix,double *origin , double *normal, BME_DuplicatorType *outImage, BME_DuplicatorType *outMask)

{
	double *imageBounds;
	double imageSpacing[4];
	//imageSpacing[0]=imageSpacing[1]=imageSpacing[2]=1; // output of ImageReslice , set to 1
	//imageSpacing[3]=0;
	vtkSmartPointer<vtkPolyData> meshReaded=  vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> meshSorted=  vtkSmartPointer<vtkPolyData>::New();

	meshReaded = computePolydata(mesh,origin , normal); //Get 2D Plane points out of 3D Mesh
	vtkSmartPointer<vtkPoints> sortedPoints2D = vtkSmartPointer<vtkPoints>::New();

	short nrpoints=meshReaded->GetNumberOfPoints() ; 
	sortedPoints2D->SetNumberOfPoints( nrpoints );
	SortPolydata(meshReaded ,sortedPoints2D ); // Sort points
	meshSorted->SetPoints(sortedPoints2D);

	BME_utility<CharPixelType,2>::BME_DuplicatorType::Pointer MaskImage2D = BME_utility<CharPixelType,2>::BME_DuplicatorType::New();;
	
	imageBounds=GetPlaneImage(vnlPlaneTransmatrix, origin,MaskImage,MaskImage2D);
  BME_ImageType2D::SpacingType Spacing = MaskImage2D->GetOutput()->GetSpacing();
  imageSpacing[0]=Spacing[0];
  imageSpacing[1]=Spacing[1];

	vtkSmartPointer<vtkPolyData> meshSortedPolyData2D= vtkSmartPointer<vtkPolyData>::New();	
	
	TransformPolydataTo2D(meshSorted, meshSortedPolyData2D,vnlPlaneTransmatrix);
	
	BME_utility<CharPixelType,2>::BME_ImageType::Pointer Image2DmeshSorted = BME_utility<CharPixelType,2>::BME_ImageType::New();;
	BME_utility<CharPixelType,2>::BME_DuplicatorType2D::Pointer StencilImage2D = BME_utility<CharPixelType,2>::BME_DuplicatorType2D::New();;
	BME_utility<CharPixelType,3>::BME_DuplicatorType3D::Pointer StencilImage3D = BME_utility<CharPixelType,3>::BME_DuplicatorType3D::New();;
	
	ConvertPolydataToITKStencil(imageSpacing,imageBounds,meshSortedPolyData2D,StencilImage2D,StencilImage3D);

	
	typedef itk::ResampleImageFilter<CharImageType2D, CharImageType2D> ResampleFilterType;
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInput( StencilImage2D->GetOutput() );
	resampler->SetReferenceImage( MaskImage2D->GetOutput());
	resampler->UseReferenceImageOn();
	try
	{
		resampler->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Reading mask image failed." << std::endl;
		std::cerr << excp << std::endl;
	}


	outImage->SetInputImage(resampler->GetOutput());
	outImage->Update();
	outMask->SetInputImage(MaskImage2D->GetOutput());
	outMask->Update();

}


template <class T,int Dim> bool BME_utility<T,Dim>::  LoadVtkObject( const std::string inFileNamePath, short PlaneNR,short cavityNR,vtkSmartPointer<vtkPolyData> mesh)
{
	
	char FileNameText[256];
	char* dummypointtype;
	vtkSmartPointer<vtkPolyDataReader> polyReader = vtkSmartPointer<vtkPolyDataReader>::New(); 
	std::vector<char*> RegionParameterName;
	RegionParameterName.assign(5,dummypointtype);
	for (int index=0;index< 5;index++)
		RegionParameterName [index]= new char[64];
	lstrcpy(RegionParameterName[LV],"LVcoords");
	lstrcpy(RegionParameterName[RV],"RVcoords");
	lstrcpy(RegionParameterName[LA],"LAcoords");
	lstrcpy(RegionParameterName[RA],"RAcoords");
	lstrcpy(RegionParameterName[Ao],"Aocoords");
	
	if(PlaneNR<10)
		sprintf(FileNameText,"%s/0%d%s.vtk",inFileNamePath.c_str(),PlaneNR,RegionParameterName[cavityNR]);
	else
		sprintf(FileNameText,"%s/%d%s.vtk",inFileNamePath.c_str(),PlaneNR,RegionParameterName[cavityNR]);
	std::fstream inFile;		

	std::fstream file(FileNameText, ios::binary | ios::in);
	if (!file.is_open())
	{
		return false;
	}
	else
	{
		file.close();

		polyReader->SetFileName( FileNameText );
		polyReader->Update();
		mesh->DeepCopy(polyReader->GetOutput());
		return true;
	}
};

template <class T,int Dim> void BME_utility<T,Dim>::  LoadTransform( const std::string inFileNamePath, vnl_matrix_fixed<double, 4, 4> landTransPlanesToWorld ,vnl_matrix_fixed<double, 4, 4>  landTransModelToWorld)
{

	char FileNameText[256];

	sprintf(FileNameText,"%s/landTransPlanesToWorld.txt",inFileNamePath.c_str());

	std::ifstream landTransPlanes (FileNameText, ios::in );
	landTransPlanesToWorld.read_ascii( landTransPlanes );
	landTransPlanes.close();


	sprintf(FileNameText,"%s/landTransModelToWorld.txt",inFileNamePath.c_str());
	std::ifstream landTransModel (FileNameText, ios::in );
	landTransModelToWorld.read_ascii( landTransModel );
	landTransModel.close();

}

template <class T,int Dim> void BME_utility<T,Dim>:: LoadPlaneTransform( const std::string inFileNamePath,short PlaneNR, vnl_matrix_fixed<double, 4, 4> &vnllandTransPlane,double *origin , double *normal)
{
	
	char FileNameText[256];

	

	if(PlaneNR<10)
		sprintf(FileNameText,"%s/0%dPlaneTransform.txt",inFileNamePath.c_str(),PlaneNR);
	else
		sprintf(FileNameText,"%s/%dPlaneTransform.txt",inFileNamePath.c_str(),PlaneNR);

	std::ifstream landTransPlane (FileNameText, ios::in );
	vnllandTransPlane.read_ascii( landTransPlane );
	landTransPlane.close();

	if(PlaneNR<10)
		sprintf(FileNameText,"%s/0%dPlaneNormal.txt",inFileNamePath.c_str(),PlaneNR);
	else
		sprintf(FileNameText,"%s/%dPlaneNormal.txt",inFileNamePath.c_str(),PlaneNR);

	std::ifstream NorFile (FileNameText, ios::out );


	NorFile >> normal[0] >>  normal[1] >> normal[2] ;

	if(PlaneNR<10)
		sprintf(FileNameText,"%s/0%dPlaneOrigin.txt",inFileNamePath.c_str(),PlaneNR);
	else
		sprintf(FileNameText,"%s/%dPlaneOrigin.txt",inFileNamePath.c_str(),PlaneNR);

	std::ifstream OrgFile (FileNameText, ios::in );

	OrgFile >> origin[0]  >>  origin[1]  >> origin[2]  ;

  // not used
	sprintf(FileNameText,"%s/landTransModelToWorldScale.txt",inFileNamePath.c_str());
	std::ifstream ScaleFile (FileNameText, ios::in );

	int SelectedFrame;
	double Scale[3];
	double Orientation[3];
	double  Position[3];

	ScaleFile >> Scale[0] ;
	ScaleFile >> Orientation[0] >>  Orientation[1]  >> Orientation[2]  ;
	ScaleFile >> Position[0] >>  Position[1] >> Position[2] ;
	ScaleFile >> SelectedFrame ;	


}

template <class T,int Dim> void BME_utility<T,Dim>::SortPolydata(vtkSmartPointer<vtkPolyData> Inpolydata , vtkSmartPointer<vtkPoints> pointsSorted)
{
	vtkSmartPointer<vtkPoints> points	= Inpolydata->GetPoints();
	vtkSmartPointer<vtkPoints> sortedPoints2D = vtkSmartPointer<vtkPoints>::New();
	int numPts=0;
	if(points)
	{			// get ordered points

		numPts = points->GetNumberOfPoints();
		//numPts/=2;
		vtkSmartPointer<vtkIdList> cellIDs = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> pointIDs = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> sortedPointIDs = vtkSmartPointer<vtkIdList>::New();
		sortedPointIDs->SetNumberOfIds( numPts );
		sortedPoints2D->SetNumberOfPoints( numPts );

		vtkIdType currentPtsID, nextPtsID, currentCellID, nextCellID;
		currentPtsID = 0;
		currentCellID = 0;

		for ( int idx = 0; idx < numPts; ++idx )
		{
	
			// cell stuff
			Inpolydata->GetPointCells( currentPtsID, cellIDs);
			if ( cellIDs->GetId(0) != currentCellID )
			{
				nextCellID = cellIDs->GetId(0);
			}
			else if (cellIDs->GetId(1) != currentCellID )
			{
				nextCellID = cellIDs->GetId(1);
			}
			else
			{
				std::cout << "Problem!!" << std::endl;
			}
		
			// points stuff
			Inpolydata->GetCellPoints(nextCellID, pointIDs); 
			if ( pointIDs->GetId(0) != currentPtsID )
			{
				nextPtsID = pointIDs->GetId(0);
			}
			else if (pointIDs->GetId(1) != currentPtsID )
			{
				nextPtsID = pointIDs->GetId(1);
			}
			else
			{
				std::cout << "Problem!!" << std::endl;
			}
			sortedPointIDs->SetId( idx, currentPtsID );
			currentPtsID = nextPtsID;
			currentCellID = nextCellID;
		}

		pointsSorted->SetNumberOfPoints( numPts );
		points->GetPoints( sortedPointIDs,pointsSorted );
	}
}

	
template <class T,int Dim> vtkSmartPointer<vtkPolyData> BME_utility<T,Dim>::ConvertImageToPolydata(BME_ImageType2D * inImage ,double *origin )
{
	//typedef float PixelType;
	//typedef itk::Image<PixelType, 2>		FImageType;
	//
	//typedef itk::CastImageFilter< BME_ImageType2D, FImageType > BME_CastType;
	//BME_CastType::Pointer doubleCaster = BME_CastType::New();
	//doubleCaster->SetInput(inImage);
	//doubleCaster->Update();

	typedef unsigned char				PixelType;
	typedef itk::Image<PixelType, 2>	ImageType;
	typedef itk::ImageFileReader<ImageType> ImageReaderType;

	ImageReaderType::Pointer reader = ImageReaderType::New();
	
//	reader->SetFileName( "Distance.mhd" );
	
	reader->SetFileName( "F:\\Data\\Patients\\TAVI\\COUVEE_20130214_001\\MHD_Data\\IM_0030\\ManualOutline\\Gerard\\test_00LaxLADistanceMapImag0.mhd" );
	ImageType::Pointer image = reader->GetOutput();

  /************************************************************************/
  /* create contour                                                       */
  /************************************************************************/
  typedef itk::ImageToVTKImageFilter<ImageType> ImageConverterITKToVTK2D;
  ImageConverterITKToVTK2D::Pointer imageConverterITKToVTK2D = ImageConverterITKToVTK2D::New();


  typedef itk::ContourExtractor2DImageFilter <ImageType>   ContourExtractor2DImageFilterType;
  ContourExtractor2DImageFilterType::Pointer contourExtractor2DImageFilter = ContourExtractor2DImageFilterType::New();
  contourExtractor2DImageFilter->SetInput( reader->GetOutput() );
  contourExtractor2DImageFilter->SetContourValue(0);



  origin[0]= reader->GetImageIO()->GetOrigin(0);
  origin[1]= reader->GetImageIO()->GetOrigin(1);
  std::cout << "There are " << contourExtractor2DImageFilter->GetNumberOfOutputs() << " contour(s)" << std::endl;

   vtkSmartPointer<vtkPolyData> vertex = vtkSmartPointer<vtkPolyData>::New();


  int numPts = contourExtractor2DImageFilter->GetOutput(0)->GetVertexList()->Size();

#if 0

	typedef itk::ImageFileReader<CharImageType2D> ImageReaderType;

	ImageReaderType::Pointer reader = ImageReaderType::New();
	
//	reader->SetFileName( "Distance.mhd" );
	
	reader->SetFileName( "F:\\Data\\Patients\\TAVI\\COUVEE_20130214_001\\MHD_Data\\IM_0030\\ManualOutline\\Gerard\\test_00LaxLADistanceMapImag0.mhd" );

  typedef itk::ContourExtractor2DImageFilter <CharImageType2D>   ContourExtractor2DImageFilterType;
  


  ContourExtractor2DImageFilterType::Pointer contourExtractor2DImageFilter = ContourExtractor2DImageFilterType::New();
  contourExtractor2DImageFilter->SetInput( reader->GetOutput());
  contourExtractor2DImageFilter->SetContourValue(0);

  //double origin[3];
  //origin[0]= reader->GetImageIO()->GetOrigin(0);
  //origin[1]= reader->GetImageIO()->GetOrigin(1);
  std::cout << "There are " << contourExtractor2DImageFilter->GetNumberOfOutputs() << " contour(s)" << std::endl;

   vtkSmartPointer<vtkPolyData> vertex = vtkSmartPointer<vtkPolyData>::New();


  int numPts = contourExtractor2DImageFilter->GetOutput(0)->GetVertexList()->Size();
  vtkSmartPointer<vtkPoints> curvepoints;
  curvepoints = vtkSmartPointer<vtkPoints>::New();
  vertex->SetPoints(curvepoints);

  vtkUnsignedCharArray *curvecolor=vtkUnsignedCharArray::New();
  curvecolor->SetNumberOfComponents(3);
  curvecolor->SetNumberOfTuples(numPts);
  vtkCellArray   *Curvelines =	vtkCellArray::New();
  double ptr[3];
//  double *ps;
  //double color[3];
  int idx=0;
  curvepoints->SetNumberOfPoints(numPts);
  int i = 0;
 // for(unsigned int i = 0; i < contourExtractor2DImageFilter->GetNumberOfOutputs(); ++i)
  {
    std::cout << "Contour " << i << ": " << std::endl;
    ContourExtractor2DImageFilterType::VertexListType::ConstIterator vertexIterator = 
      contourExtractor2DImageFilter->GetOutput(i)->GetVertexList()->Begin();
    while(vertexIterator != contourExtractor2DImageFilter->GetOutput(i)->GetVertexList()->End())
    {
      std::cout << vertexIterator->Value() << std::endl;
	  ptr[0]=vertexIterator->Value()[0]+origin[0];
	  ptr[1]=vertexIterator->Value()[1]+origin[1];
	  ptr[2]=0;
	  ++vertexIterator;
 
	 		
		  curvepoints->InsertPoint(idx,ptr);
		  curvecolor->SetTuple3(idx ,255,138,126);
		  //else
		  //	curvecolor->SetTuple3(i ,color2[0], color2[1],color2[2]);
		  if ( idx > 0 )
		  {
			  vtkIdType	cell[2]	= {idx-1,idx};
			  Curvelines->InsertNextCell((vtkIdType)2,(vtkIdType*)	cell);
		  }
		  idx++;
	  }
	  vertex->SetPoints(curvepoints);
	  vertex->SetLines(Curvelines);
	  vertex->GetPointData()->SetScalars(curvecolor);

	}
#endif
	return vertex;
	
}


template <class T,int Dim> vtkSmartPointer<vtkPolyData> BME_utility<T,Dim>::computePolydata( const vtkSmartPointer<vtkPolyData> mesh,double *origin , double *normal)
{

	vtkSmartPointer<vtkPlane> polyCutPlane =	vtkSmartPointer<vtkPlane>::New();

	polyCutPlane->SetOrigin( origin );
	polyCutPlane->SetNormal( normal );
	polyCutPlane->Modified();

	vtkSmartPointer<vtkCutter> cutter =	vtkSmartPointer<vtkCutter>::New();
	cutter->SetCutFunction( polyCutPlane );
	cutter->SetInput( mesh);
	cutter->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> splitter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	splitter->SetInputConnection( cutter->GetOutputPort() );
	splitter->SetExtractionModeToLargestRegion();// changed by ahaak 17-Nov-2013
  /*splitter->SetExtractionModeToSpecifiedRegions();
	splitter->AddSpecifiedRegion( 0 ); changed by ahaak 17-Nov-2013 */
	splitter->Modified();
	splitter->Update();

	int numCavities = splitter->GetNumberOfExtractedRegions();

	// do first cavity outside of the loop
	// splitter->AddSpecifiedRegion( 0 ); changed by ahaak 17-Nov-2013
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputConnection( splitter->GetOutputPort() );
	cleaner->PointMergingOff();
	cleaner->Update();

	return cleaner->GetOutput();
}

	
template <class T,int Dim> void  BME_utility<T,Dim>::TransformPolydataTo2D( vtkSmartPointer<vtkPolyData> input, vtkSmartPointer<vtkPolyData> output,vnl_matrix_fixed<double, 4, 4> &vnlTransmatrix)
{
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
	for( int xx=0;xx<4;xx++)
		for( int yy=0;yy<4;yy++)
			matrix->SetElement(xx,yy,vnlTransmatrix[xx][yy]);
	transform->SetMatrix( matrix );
	transform->Update();
	vtkSmartPointer<vtkPoints>		Points2D		= vtkSmartPointer<vtkPoints>::New() ;
	double * point2D;
	double nrpoints =input->GetNumberOfPoints() ;
	Points2D->SetNumberOfPoints( nrpoints);
	vtkCellArray   *Curvelines =	vtkCellArray::New();
	for ( int idx = 0; idx <(nrpoints); ++idx )
	{
		point2D=transform->GetInverse()->TransformDoublePoint(input->GetPoint(idx));
		point2D[2]=0;
		Points2D->SetPoint( idx, point2D );
		if ( idx > 0 )
		{
		  vtkIdType	cell[2]	= {idx-1,idx};
		  Curvelines->InsertNextCell((vtkIdType)2,(vtkIdType*)	cell);
		}
	}
	
	output->SetPoints(Points2D);
	output->SetLines(Curvelines);
}

template <class T,int Dim> vtkSmartPointer<vtkImageStencil> BME_utility<T,Dim>::ConvertPolydataToStencil(double *spacing, double *bounds,vtkSmartPointer<vtkPolyData> polyData)
{
	vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
	stripper->SetInput(polyData); // valid polyData
	stripper->Update();
	// that's our polyData
	vtkSmartPointer<vtkPolyData> strippedPolyData = stripper->GetOutput();

	// prepare the binary image's voxel grid
	vtkSmartPointer<vtkImageData> whiteImage =
		vtkSmartPointer<vtkImageData>::New();
	
	//strippedPolyData->GetBounds(bounds);
	//spacing[0] = 0.5;
	//spacing[1] = 0.5;
	//spacing[2] = 0.5;

	//BME_ImageReslice->GetOutput()->GetBounds(bounds);
	//BME_ImageReslice->GetOutput()->GetSpacing(spacing);
	spacing[2] = 1.0;

	whiteImage->SetSpacing(spacing);

	// compute dimensions
	int dim[3];
	for (int i = 0; i < 3; i++)
	{
		dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) /
			spacing[i])) + 1;
		if (dim[i] < 1)
			dim[i] = 1;
	}
	dim[2]=2;
	whiteImage->SetDimensions(dim);
	whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);
	double origin[3];
	// NOTE: I am not sure whether or not we had to add some offset!
	origin[0] = bounds[0];// + spacing[0] / 2;
	origin[1] = bounds[2];// + spacing[1] / 2;
	origin[2] = bounds[4];// + spacing[2] / 2;
	whiteImage->SetOrigin(origin);
#if VTK_MAJOR_VERSION <= 5
	whiteImage->SetScalarTypeToUnsignedChar();
	whiteImage->AllocateScalars();
#else
	whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
	// fill the image with foreground voxels:
	unsigned char inval = 255;
	unsigned char outval = 0;
	vtkIdType count = whiteImage->GetNumberOfPoints();
	for (vtkIdType i = 0; i < count; ++i)
	{
		whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
	}


  // sweep polygonal data (this is the important thing with contours!)
  vtkSmartPointer<vtkLinearExtrusionFilter> extruder =
    vtkSmartPointer<vtkLinearExtrusionFilter>::New();
  vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();

  int numPts=polyData->GetPoints()->GetNumberOfPoints() ;
  if(numPts)
  {
#if VTK_MAJOR_VERSION <= 5
    extruder->SetInput(polyData);
#else
    extruder->SetInputData(polyData);
#endif

    extruder->SetScaleFactor(1.0);
    extruder->SetExtrusionTypeToNormalExtrusion();
    extruder->SetVector(0, 0, 1);

    extruder->Update();
    // polygonal data --> image stencil:
    pol2stenc->SetTolerance(0); // important if extruder->SetVector(0, 0, 1) !!!
    pol2stenc->SetInputConnection(extruder->GetOutputPort());
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();
  }
  else
    for (vtkIdType i = 0; i < count; ++i)
    {
      whiteImage->GetPointData()->GetScalars()->SetTuple1(i, outval);
    }



    // cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imgstenc =
      vtkSmartPointer<vtkImageStencil>::New();


    origin[0] = bounds[0];// + spacing[0] / 2;
    origin[1] = bounds[2];// + spacing[1] / 2;
    origin[2] = bounds[4];// + spacing[2] / 2;

    whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
    imgstenc->SetInput(whiteImage);
    if(numPts)
      imgstenc->SetStencil(pol2stenc->GetOutput());
#else
    imgstenc->SetInputData(whiteImage);
    if(numPts)
      imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(outval);
    imgstenc->Update();


	vtkSmartPointer<vtkImageAccumulate> imageAccumulate =	vtkSmartPointer<vtkImageAccumulate>::New();
#if VTK_MAJOR_VERSION <= 5
	imageAccumulate->SetStencil(pol2stenc->GetOutput());
	imageAccumulate->SetInputConnection(whiteImage->GetProducerPort());
#else
	imageAccumulate->SetStencilData(pol2stenc->GetOutput());
	imageAccumulate->SetInputData(whiteImage);
#endif
	imageAccumulate->Update();
	//std::cout << "Voxel count: " << imageAccumulate->GetVoxelCount() << std::endl;
	
	//itk::Image<T,3>::Pointer outimage= ImageConverterVTKToITK3D(imgstenc->GetOutput());

	return imgstenc;

};

template <class T,int Dim> void BME_utility<T,Dim>::ConvertPolydataToITKStencil(double *spacing, double *bounds,vtkSmartPointer<vtkPolyData> polyData,BME_DuplicatorType2D *outImage2D,BME_DuplicatorType3D *outImage3D)
{
	vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
	stripper->SetInput(polyData); // valid polyData
	stripper->Update();
	// that's our polyData
	vtkSmartPointer<vtkPolyData> strippedPolyData = stripper->GetOutput();

	// prepare the binary image's voxel grid
	vtkSmartPointer<vtkImageData> whiteImage =
		vtkSmartPointer<vtkImageData>::New();
	
	//strippedPolyData->GetBounds(bounds);
	//spacing[0] = 0.5;
	//spacing[1] = 0.5;
	//spacing[2] = 0.5;

	//BME_ImageReslice->GetOutput()->GetBounds(bounds);
	//BME_ImageReslice->GetOutput()->GetSpacing(spacing);
	spacing[2] = 1.0;

	whiteImage->SetSpacing(spacing);

	// compute dimensions
	int dim[3];
	for (int i = 0; i < 3; i++)
	{
		dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) /
			spacing[i])) + 1;
		if (dim[i] < 1)
			dim[i] = 1;
	}
	dim[2]=2;
	whiteImage->SetDimensions(dim);
	whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);
	double origin[3];
	// NOTE: I am not sure whether or not we had to add some offset!
	origin[0] = bounds[0];// + spacing[0] / 2;
	origin[1] = bounds[2];// + spacing[1] / 2;
	origin[2] = bounds[4];// + spacing[2] / 2;
	whiteImage->SetOrigin(origin);
#if VTK_MAJOR_VERSION <= 5
	whiteImage->SetScalarTypeToUnsignedChar();
	whiteImage->AllocateScalars();
#else
	whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
	// fill the image with foreground voxels:
	unsigned char inval = 1;
	unsigned char outval = 0;
	vtkIdType count = whiteImage->GetNumberOfPoints();
	for (vtkIdType i = 0; i < count; ++i)
	{
		whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
	}


  // sweep polygonal data (this is the important thing with contours!)
  vtkSmartPointer<vtkLinearExtrusionFilter> extruder =
    vtkSmartPointer<vtkLinearExtrusionFilter>::New();
  vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();

  int numPts=polyData->GetPoints()->GetNumberOfPoints() ;
  if(numPts)
  {
#if VTK_MAJOR_VERSION <= 5
    extruder->SetInput(polyData);
#else
    extruder->SetInputData(polyData);
#endif

    extruder->SetScaleFactor(1.0);
    extruder->SetExtrusionTypeToNormalExtrusion();
    extruder->SetVector(0, 0, 1);

    extruder->Update();
    // polygonal data --> image stencil:
    pol2stenc->SetTolerance(0); // important if extruder->SetVector(0, 0, 1) !!!
    pol2stenc->SetInputConnection(extruder->GetOutputPort());
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();
  }
  else
    for (vtkIdType i = 0; i < count; ++i)
    {
      whiteImage->GetPointData()->GetScalars()->SetTuple1(i, outval);
    }



    // cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imgstenc =
      vtkSmartPointer<vtkImageStencil>::New();


    origin[0] = bounds[0];// + spacing[0] / 2;
    origin[1] = bounds[2];// + spacing[1] / 2;
    origin[2] = bounds[4];// + spacing[2] / 2;

    whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
    imgstenc->SetInput(whiteImage);
    if(numPts)
      imgstenc->SetStencil(pol2stenc->GetOutput());
#else
    imgstenc->SetInputData(whiteImage);
    if(numPts)
      imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(outval);
    imgstenc->Update();


	vtkSmartPointer<vtkImageAccumulate> imageAccumulate =	vtkSmartPointer<vtkImageAccumulate>::New();
#if VTK_MAJOR_VERSION <= 5
	imageAccumulate->SetStencil(pol2stenc->GetOutput());
	imageAccumulate->SetInputConnection(whiteImage->GetProducerPort());
#else
	imageAccumulate->SetStencilData(pol2stenc->GetOutput());
	imageAccumulate->SetInputData(whiteImage);
#endif
	imageAccumulate->Update();
	// std::cout << "Voxel count: " << imageAccumulate->GetVoxelCount() << std::endl;
	
	//itk::Image<T,3>::Pointer outimage= ImageConverterVTKToITK3D(imgstenc->GetOutput());

	

	typedef itk::VTKImageToImageFilter<BME_utility<CharPixelType,2>::BME_ImageType3D> BME_ImageConverterVTKToITK3D;

	BME_utility<CharPixelType,3>::BME_ImageConverterVTKToITK::Pointer imageConverterImgstencVTKToITK3D = BME_utility<CharPixelType,3>::BME_ImageConverterVTKToITK::New();
	imageConverterImgstencVTKToITK3D->SetInput(imgstenc->GetOutput());
	imageConverterImgstencVTKToITK3D->Update();


	BME_utility<signedCharPixelType,3>::BME_signedDuplicatorType::Pointer DistanceMapUser13D = BME_utility<signedCharPixelType,3>::BME_signedDuplicatorType::New();
	

	//GetDistanceMap(imageConverterImgstencVTKToITK3D->GetOutput(),DistanceMapUser13D);
	
//	BME_utility<CharPixelType,3>::WriteMHD("F:/Data/Patients/Tavi/COUVEE_20130214_001/MHD_Data/IM_0030/DiceCalculations/User2/userDis1.mhd",DistanceMapUser13D->GetOutput());

	typedef itk::ExtractImageFilter<CharImageType3D, CharImageType2D> ImageFilter3Dto2D; //Get Slice out of 3D stack
			
	CharImageType3D::RegionType iRegion3D = imageConverterImgstencVTKToITK3D->GetOutput()->GetLargestPossibleRegion();
	CharImageType3D::SizeType   iSize3D = iRegion3D.GetSize();
	CharImageType3D::RegionType extractRegion = iRegion3D;
	CharImageType3D::SizeType   extractSize = iSize3D;
	CharImageType3D::IndexType  extractIndex = extractRegion.GetIndex();
	
	ImageFilter3Dto2D::Pointer imageSlicers3Dto2D = ImageFilter3Dto2D::New();
	extractSize[2] = 0;
	extractIndex[0] = extractIndex[1] = 0;
	extractIndex[2]  = 1;
	extractRegion.SetSize(extractSize);
	extractRegion.SetIndex(extractIndex);
				

	imageSlicers3Dto2D->SetExtractionRegion(extractRegion);
	imageSlicers3Dto2D->SetDirectionCollapseToIdentity();
	imageSlicers3Dto2D->SetInput(imageConverterImgstencVTKToITK3D->GetOutput());
	imageSlicers3Dto2D->Update();
	//outImage=imageSlicers3Dto2D->GetOutput();
	
	
		//outImage=imageConverterImgstencVTKToITK3D->GetOutput();

	outImage2D->SetInputImage(imageSlicers3Dto2D->GetOutput());
	outImage2D->Update();
	outImage3D->SetInputImage(imageConverterImgstencVTKToITK3D->GetOutput());
	outImage3D->Update();

};


template <class T,int Dim> double BME_utility<T,Dim>::GetVoxelCount(itk::Image<T,Dim> * Image)
{
	BME_ImageConverterITKToVTK::Pointer imageConverterITKToVTK = BME_ImageConverterITKToVTK::New();
	imageConverterITKToVTK->SetInput(Image);
	imageConverterITKToVTK->Update();

	vtkSmartPointer<vtkImageAccumulate> imageAccumulate =	vtkSmartPointer<vtkImageAccumulate>::New();
#if VTK_MAJOR_VERSION <= 5
	imageAccumulate->SetInput(imageConverterITKToVTK->GetOutput());
#else
	imageAccumulate->SetInputData(imageConverterITKToVTK->GetOutput());
#endif
	imageAccumulate->IgnoreZeroOn  (  ) ;
	imageAccumulate->Update();
	std::cout << "Voxel count:2d " << imageAccumulate->GetVoxelCount() << std::endl;
	return  imageAccumulate->GetVoxelCount() ;

};

template <class T,int Dim> double * BME_utility<T,Dim>::GetPlaneImage(vnl_matrix_fixed<double, 4, 4> vnlPlaneTransmatrix,double * origin,vtkSmartPointer <vtkImageData> InputImage,BME_DuplicatorType *outImage)
{
	vtkSmartPointer<vtkMatrix4x4> Planematrix = vtkSmartPointer<vtkMatrix4x4>::New();
	double *ImageBounds = new double [6];
	double ImageSpacing[6];
	for( int xx=0;xx<4;xx++)
		for( int yy=0;yy<4;yy++)
			Planematrix->SetElement(xx,yy,vnlPlaneTransmatrix[xx][yy]);

	vtkImageReslice * BMET_ImageReslice=vtkImageReslice::New();
	BMET_ImageReslice->SetResliceAxesOrigin(0,0,0);
	BMET_ImageReslice->SetOutputDimensionality(2);
	// BMET_ImageReslice->SetOutputSpacing(1,1,1);
  BMET_ImageReslice->TransformInputSamplingOn();
	BMET_ImageReslice->InterpolateOn();
	BMET_ImageReslice->AutoCropOutputOn();
	BMET_ImageReslice->Modified();
	BMET_ImageReslice->SetInterpolationModeToLinear();

#if VTK_MAJOR_VERSION <= 5
	BMET_ImageReslice->SetInput(InputImage);
#else
	BMET_ImageReslice->SetInputData(imageReader->GetOutput());
#endif

	BMET_ImageReslice->GetResliceAxes()->DeepCopy( Planematrix); 
	BMET_ImageReslice->SetResliceAxesOrigin(origin);  
	BMET_ImageReslice->Modified();
	BMET_ImageReslice->Update();
	int ext[6];
	BMET_ImageReslice ->GetOutput()->GetExtent( ext );
	int offsetx = 28;
	int offsety = 28;
	ext[0] = -offsetx;
	ext[1] += offsetx;
	ext[2] = -offsety;
	ext[3] += offsety;
	ext[0] = (ext[0]*4)/4;
	ext[1] = (ext[1]*4)/4;
	ext[2] = (ext[2]*4)/4;
	ext[3] = (ext[3]*4)/4;

	//BMET_ImageReslice->SetOutputExtent( ext );;
	BMET_ImageReslice->Modified();
	BMET_ImageReslice->Update();
	BMET_ImageReslice->GetOutput()->GetBounds(ImageBounds);
	BMET_ImageReslice->GetOutput()->GetSpacing(ImageSpacing);
  
	BME_utility<CharPixelType,2>::BME_ImageConverterVTKToITK::Pointer imageConverterImgstencVTKToITK2D = BME_utility<CharPixelType,2>::BME_ImageConverterVTKToITK::New();
	imageConverterImgstencVTKToITK2D->SetInput(BMET_ImageReslice->GetOutput());
	imageConverterImgstencVTKToITK2D->Update();

	outImage->SetInputImage(imageConverterImgstencVTKToITK2D->GetOutput());
	outImage->Update();

	return (double *) ImageBounds;
};


template <class T,int Dim> void  BME_utility<T,Dim>:: GetUnionImage(BME_ImageType *inImage1,BME_ImageType *inImage2,BME_DuplicatorType *outImage)
{
	BME_utility<CharPixelType,Dim>::BME_OrImageFiltertype::Pointer Or_ImageFilter = BME_utility<CharPixelType,Dim>::BME_OrImageFiltertype::New();
	Or_ImageFilter->SetInput( 0, inImage1); 
	Or_ImageFilter->SetInput( 1, inImage2 ); 
	Or_ImageFilter->Update(); //Masking
	outImage->SetInputImage(Or_ImageFilter->GetOutput());
	outImage->Update();
};

template <class T,int Dim> void  BME_utility<T,Dim>:: GetIntersectionImage(BME_ImageType *inImage1,BME_ImageType *inImage2,BME_DuplicatorType *outImage)
{
	BME_utility<CharPixelType,Dim>::BME_AndImageFiltertype::Pointer And_ImageFilter = BME_utility<CharPixelType,Dim>::BME_AndImageFiltertype::New();
	And_ImageFilter->SetInput( 0, inImage1); 
	And_ImageFilter->SetInput( 1, inImage2 ); 
	And_ImageFilter->Update(); //Masking
	outImage->SetInputImage(And_ImageFilter->GetOutput());
	outImage->Update();
};


template <class T,int Dim> void  BME_utility<T,Dim>:: GetMaskedImage(BME_ImageType *inImage1,BME_ImageType *inImage2,BME_DuplicatorType *outImage)
{
	BME_utility<CharPixelType,Dim>::BME_AndImageFiltertype::Pointer And_ImageFilter = BME_utility<CharPixelType,Dim>::BME_AndImageFiltertype::New();
	And_ImageFilter->SetInput( 0, inImage1); 
	And_ImageFilter->SetInput( 1, inImage2 ); 
	And_ImageFilter->Update(); //Masking
	outImage->SetInputImage(And_ImageFilter->GetOutput());
	outImage->Update();
};

template <class T,int Dim> void  BME_utility<T,Dim>::GetThresholdImage(BME_INT8ImageType *inImage,BME_DuplicatorType *outImage)
{

	//typedef itk::CastImageFilter< BME_ImageType, BME_INT8ImageType > BME_CastType;
	//BME_CastType::Pointer SignedCaster = BME_CastType::New();
	//SignedCaster->SetInput(ITKStencilUser1);
	//SignedCaster->Update();

	BME_BinaryThreshold::Pointer  binaryThreshold = BME_BinaryThreshold ::New();
	const short lthreshold=SHRT_MIN;
	short uthreshold=1;
	binaryThreshold->SetInput(inImage); 
	binaryThreshold->SetLowerThreshold (lthreshold);
	binaryThreshold->SetUpperThreshold (uthreshold) ;
	binaryThreshold->SetInsideValue(1);
	binaryThreshold->SetOutsideValue(0);
	binaryThreshold->Update();
	outImage->SetInputImage(binaryThreshold->GetOutput());
    outImage->Update();
	
};

template <class T,int Dim> void BME_utility<T,Dim>::GetDistanceMap(BME_ImageType *inImage,BME_signedDuplicatorType *outImage)
{
	typedef itk::SignedDanielssonDistanceMapImageFilter<BME_ImageType,BME_INT8ImageType>  BME_DistanceMapImageFilter;
	BME_DistanceMapImageFilter::Pointer DistanceMapImageFilter = BME_DistanceMapImageFilter::New();
	DistanceMapImageFilter->SetInput(inImage);
	DistanceMapImageFilter->SetInsideIsPositive(false);
	DistanceMapImageFilter->SetUseImageSpacing(true);
	DistanceMapImageFilter->SquaredDistanceOff ();
	DistanceMapImageFilter->Update();
	outImage->SetInputImage(DistanceMapImageFilter->GetOutput());
    outImage->Update();
}



template <class T,int Dim> void BME_utility<T,Dim>::WriteMHD(const std::string inFileNamePath,BME_ImageType *inImage)
{
	BME_ImageConverterITKToVTK::Pointer imageConverterITKToVTK = BME_ImageConverterITKToVTK::New();
	imageConverterITKToVTK->SetInput(inImage);
	imageConverterITKToVTK->Update();
	vtkSmartPointer<vtkMetaImageWriter> imageWriter =vtkSmartPointer<vtkMetaImageWriter>::New();
	imageWriter->SetFileName(inFileNamePath.c_str());
	imageWriter->SetInput(imageConverterITKToVTK->GetOutput());
	imageWriter->Write();


	//typedef itk::ImageFileWriter<BME_ImageType> CharImageWriterType;
	//CharImageWriterType::Pointer writer = CharImageWriterType::New();
	//writer->SetFileName(inFileNamePath.c_str());
	//writer->SetInput(inImage);
	//writer->Write();	
}
#if 0 //Not OK Yet 

template <class T,int Dim> vtkSmartPointer<vtkImageData> BME_utility<T,Dim>::ImageConverterITKToVTK3D( const itk::Image<T,3> InputImage )
{
	
	BME_ImageConverterITKToVTK3D::Pointer 	BME_ConverterITKToVTK ;	
	BME_ConverterITKToVTK=ImageConverterITKToVTK3D::New();
	BME_ConverterITKToVTK->SetInput(InputImage);
	BME_ConverterITKToVTK->Update();
	
	return BME_ConverterITKToVTK->GetOutput();
}

template <class T,int Dim> vtkSmartPointer<vtkImageData> BME_utility<T,Dim>::ImageConverterITKToVTK2D( const itk::Image<T,2> InputImage )
{
	BME_ImageConverterITKToVTK2D::Pointer 	BME_ConverterITKToVTK ;	
	BME_ConverterITKToVTK=ImageConverterITKToVTK3D::New();
	BME_ConverterITKToVTK->SetInput(InputImage);
	BME_ConverterITKToVTK->Update();
	return BME_ConverterITKToVTK->GetOutput();
}


template <class T,int Dim> itk::Image<T,3> BME_utility<T,Dim>::ImageConverterVTKToITK3D(   vtkSmartPointer <vtkImageData> InputImage )
{

	itk::VTKImageToImageFilter<itk::Image<T,3>>::Pointer BME_ConverterVTKToITK ;	
	BME_ConverterVTKToITK=ImageConverterITKToVTK3D::New();
	BME_ConverterVTKToITK->SetInput(InputImage);
	BME_ConverterVTKToITK->Update();
	return BME_ConverterVTKToITK->GetOutput();
}

template <class T,int Dim> itk::Image<T,2> BME_utility<T,Dim>::ImageConverterVTKToITK2D(   vtkSmartPointer <vtkImageData> InputImage )
{
	BME_ImageConverterVTKToITK2D::Pointer 	BME_ConverterVTKToITK ;	
	BME_ConverterVTKToITK=ImageConverterITKToVTK3D::New();
	BME_ConverterVTKToITK->SetInput(InputImage);
	BME_ConverterVTKToITK->Update();

	return BME_ConverterVTKToITK->GetOutput();
}
#endif

template <class T,int Dim> BME_utility<T,Dim>::BME_utility(void)
{
	dim = 3;
}
template <class T,int Dim> BME_utility<T,Dim>::~BME_utility(void)
{

}

template <class T,int Dim> vnl_matrix<T> BME_utility<T,Dim>::vnlShapeVectorToMatrix( const vnl_vector<T> vector )
{
	int numPoints = vector.size() / 3;

	vnl_matrix<T> matrix;
	matrix.set_size( numPoints, 3 );

	int idx = 0;

	for ( int idxRow = 0; idxRow < numPoints; ++idxRow )
	{
		for ( int idxColumn = 0; idxColumn < 3; ++idxColumn )
		{
			matrix[idxRow][idxColumn] = vector[idx];
			++idx;
		}
	}
	return matrix;
}
template <class T,int Dim> vnl_vector<T> BME_utility<T,Dim>::vnlShapeMatrixToVector( const vnl_matrix<T> matrix )
{
	vnl_vector<T> vector;
	vector.set_size( matrix.rows() * matrix.columns() );

	int idx = 0;

	for ( int idxRow = 0; idxRow < matrix.rows(); ++idxRow )
	{
		for ( int idxColumn = 0; idxColumn < matrix.columns(); ++idxColumn )
		{
			vector[idx] = matrix[idxRow][idxColumn];
			++idx;
		}
	}
	return vector;
}
template <class T,int Dim> vnl_matrix<T>  BME_utility<T,Dim>::vtkPointsToVnlPointMatrix( const vtkSmartPointer<vtkPoints> inPts )
{
	vnl_matrix<T> outPts;
	outPts.set_size( inPts->GetNumberOfPoints(), dim );

	double tempPts[3];
	for ( int idxRow = 0; idxRow < inPts->GetNumberOfPoints(); ++idxRow )
	{
		inPts->GetPoint( idxRow, tempPts );
		outPts.set_row( idxRow, tempPts );
	}
	return outPts;
}
template <class T,int Dim> vtkSmartPointer<vtkPoints>  BME_utility<T,Dim>::vnlPointMatrixToVtkPoints( const vnl_matrix<T> inPts )
{
	vtkSmartPointer<vtkPoints> outPts = vtkSmartPointer<vtkPoints>::New();
	outPts->SetNumberOfPoints( inPts.rows() );

	double tempPts[3];
	vnl_vector<T> tempV;

	for ( int idxRow = 0; idxRow < inPts.rows(); ++idxRow )
	{
		tempV = inPts.get_row( idxRow );
		for ( int idx = 0; idx < dim; ++idx )
		{
			tempPts[idx] = tempV[idx];
		}
		outPts->SetPoint( idxRow, tempPts );
	}
	return outPts;
}
template <class T,int Dim> vtkSmartPointer<vtkPolyData> BME_utility<T,Dim>::vnlPointMatrixToVtkMesh( const vnl_matrix<T> S, vtkSmartPointer<vtkPolyData> templateMesh )
{
	vtkSmartPointer<vtkPolyData> outMesh = vtkSmartPointer<vtkPolyData>::New();
	outMesh->DeepCopy( templateMesh );

	if ( outMesh->GetNumberOfPoints() == S.size() / dim )
	{
		// read points into mesh
		vtkSmartPointer<vtkPoints> vtkPts = vnlPointMatrixToVtkPoints( S );
		outMesh->SetPoints( vtkPts );
	} 
	else
	{
		std::cout << "Miss match of number of shape and template mesh points." << std::endl;
	}
	return outMesh;
}
template <class T,int Dim> vnl_matrix<T> BME_utility<T,Dim>::vtkMeshToVnlPointMatrix( vtkSmartPointer<vtkPolyData> poly )
{
	vtkSmartPointer<vtkPoints> points = poly->GetPoints();
	int numPoints = points->GetNumberOfPoints();

	vnl_matrix<T> outPts;

	if ( numPoints > 0 )
	{
		outPts = vtkPointsToVnlPointMatrix( points );
	}
	return outPts;
}
template <class T,int Dim> vnl_vector<T> BME_utility<T,Dim>::computeOverlap( const vtkSmartPointer<vtkPolyData> mesh, const CharImagePointerType3D image )
{
	vnl_vector<T> returnValues;
	vtkSmartPointer<vtkPolyDataConnectivityFilter> splitter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	splitter->SetInput( mesh );
	splitter->SetExtractionModeToSpecifiedRegions();
	splitter->Update();
	int numCavities = splitter->GetNumberOfExtractedRegions();
	returnValues.set_size( numCavities );

	// do first cavity outside of the loop
	splitter->AddSpecifiedRegion( 0 );
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputConnection( splitter->GetOutputPort() );
	cleaner->PointMergingOff();
	cleaner->Update();

	CharImageType3D::Pointer meshImage = CharImageType3D::New();
	meshImage->SetRegions( image->GetLargestPossibleRegion() );
	meshImage->SetSpacing( image->GetSpacing() );
	meshImage->SetOrigin( image->GetOrigin() );
	meshImage->Allocate();
	meshImage->FillBuffer( 0 );

	vtkMeshToItkImage( cleaner->GetOutput(), meshImage );

	typedef itk::ImageRegionConstIterator<CharImageType3D> CharImageIteratorType;
	CharImageIteratorType meshIt( meshImage, meshImage->GetLargestPossibleRegion() );
	CharImageIteratorType maskIt( image, image->GetLargestPossibleRegion() );
	meshIt.GoToBegin();
	maskIt.GoToBegin();
	double overlapVolume = 0.0;

	while( !meshIt.IsAtEnd() )
	{
		overlapVolume += meshIt.Get() * maskIt.Get();
		++meshIt;
		++maskIt;
	}
	returnValues[0] = overlapVolume;
	//// debug
	//typedef itk::ImageFileWriter<CharImageType3D> CharImageWriterType;
	//CharImageWriterType::Pointer writer = CharImageWriterType::New();
	//writer->SetInput( meshImage );
	//writer->SetFileName( "meshImage0.mhd" );
	//writer->Update();

	//vtkSmartPointer<vtkPolyDataWriter> polyWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	//polyWriter->SetInput( cleaner->GetOutput() );
	//polyWriter->SetFileName( "mesh0.vtk" );
	//polyWriter->Update();
	//std::string fileName;
	//std::stringstream convert;
	//// end debug

	
	for( int idx = 1; idx < numCavities; ++idx )
	{
		// do other cavity inside of the loop
		splitter->DeleteSpecifiedRegion( idx - 1 );
		splitter->AddSpecifiedRegion( idx );
		cleaner->Update();

		meshImage->FillBuffer( 0 );
		
		vtkMeshToItkImage( cleaner->GetOutput(), meshImage );

		CharImageIteratorType meshIt( meshImage, meshImage->GetLargestPossibleRegion() );
		CharImageIteratorType maskIt( image, image->GetLargestPossibleRegion() );
		meshIt.GoToBegin();
		maskIt.GoToBegin();
		double overlapVolume = 0.0;

		while( !meshIt.IsAtEnd() )
		{
			overlapVolume += meshIt.Get() * maskIt.Get();
			++meshIt;
			++maskIt;
		}
		returnValues[idx] = overlapVolume;

		//// debug
		//std::string fileName;
		//std::stringstream convert;
		//writer->SetInput( meshImage );
		//
		//convert << "meshImage" << idx << ".mhd";
		//fileName = convert.str();
		//convert.flush();

		//writer->SetFileName( fileName.c_str() );
		//writer->Update();

		//// polyWriter->SetInput( cleaner->GetOutput() );
		//std::string fileName2;
		//std::stringstream convert2;
		//convert2 << "mesh" << idx << ".mhd";
		//fileName2 = convert2.str();
		//polyWriter->SetFileName( fileName2.c_str() );
		//polyWriter->Update();
		//// end debug
	}

	return returnValues;
}
template <class T,int Dim> vnl_vector<T> BME_utility<T,Dim>::computeVolume( const vtkSmartPointer<vtkPolyData> mesh, const CharImagePointerType3D image )
{
	vnl_vector<T> returnValues;
	vtkSmartPointer<vtkPolyDataConnectivityFilter> splitter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	splitter->SetInput( mesh );
	splitter->SetExtractionModeToSpecifiedRegions();
	splitter->Update();
	int numCavities = splitter->GetNumberOfExtractedRegions();
	returnValues.set_size( numCavities );
			

	// do first cavity outside of the loop
	splitter->AddSpecifiedRegion( 0 );
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputConnection( splitter->GetOutputPort() );
	cleaner->PointMergingOff();
	cleaner->Update();

	CharImageType3D::Pointer meshImage = CharImageType3D::New();
	meshImage->SetSpacing( image->GetSpacing() );
	
	// find size for mesh image
	CharImageType3D::RegionType region;
	CharImageType3D::SizeType size;
	CharImageType3D::IndexType index;
	CharImageType3D::PointType origin;
	CharImageType3D::SpacingType spacing = image->GetSpacing(); 
	index.Fill( 0 );
	T temp[6];
	mesh->GetBounds(temp);
	std::cout << "bounds: x1" << temp[0] << "\t x2" <<  temp[1] 
	<< "\t y1" <<  temp[2] << "\t y2" <<  temp[3] 
	<< "\t z1" <<  temp[4] << "\t z2" <<  temp[5]	<<  std::endl;
	
	size[0] = ceil( (temp[1]-temp[0])/spacing[0] ) + 2;
	size[1] = ceil( (temp[3]-temp[2])/spacing[1] ) + 2;
	size[2] = ceil( (temp[5]-temp[4])/spacing[2] ) + 2;

	region.SetIndex( index );
	region.SetSize( size );
	meshImage->SetRegions( region );

	origin[0] = temp[0] - spacing[0];
	origin[1] = temp[2] - spacing[0];
	origin[2] = temp[4] - spacing[0];
	
	meshImage->SetOrigin( origin );
	meshImage->Allocate();
	meshImage->FillBuffer( 0 );

	vtkMeshToItkImage( cleaner->GetOutput(), meshImage );
	
	// debug
	typedef itk::ImageFileWriter<CharImageType3D> CharImageWriterType;
	CharImageWriterType::Pointer writer = CharImageWriterType::New();
	writer->SetInput( meshImage );
	writer->SetFileName( "meshImage0.mhd" );
	writer->Update();

	vtkSmartPointer<vtkPolyDataWriter> polyWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	polyWriter->SetInput( cleaner->GetOutput() );
	polyWriter->SetFileName( "mesh0.vtk" );
	polyWriter->Update();
	std::string fileName;
	std::stringstream convert;
	// end debug
	
	typedef itk::ImageRegionConstIterator<CharImageType3D> CharImageIteratorType;
	CharImageIteratorType meshIt( meshImage, meshImage->GetLargestPossibleRegion() );
	meshIt.GoToBegin();
	double volume = 0.0;

	while( !meshIt.IsAtEnd() )
	{
		volume += meshIt.Get();
		++meshIt;
	}
	returnValues[0] = volume;
	


	for( int idx = 1; idx < numCavities; ++idx )
	{
		// do other cavity inside of the loop
		splitter->DeleteSpecifiedRegion( idx - 1 );
		splitter->AddSpecifiedRegion( idx );
		cleaner->Update();

		meshImage->FillBuffer( 0 );

		vtkMeshToItkImage( cleaner->GetOutput(), meshImage );

		CharImageIteratorType meshIt( meshImage, meshImage->GetLargestPossibleRegion() );
		meshIt.GoToBegin();
		double volume = 0.0;

		while( !meshIt.IsAtEnd() )
		{
			volume += meshIt.Get();
			++meshIt;
		}
		returnValues[idx] = volume;

		// debug
		std::string fileName;
		std::stringstream convert;
		writer->SetInput( meshImage );
		
		convert << "meshImage" << idx << ".mhd";
		fileName = convert.str();
		
		writer->SetFileName( fileName.c_str() );
		writer->Update();

		// polyWriter->SetInput( cleaner->GetOutput() );
		std::string fileName2;
		std::stringstream convert2;
		convert2 << "mesh" << idx << ".vtk";
		fileName2 = convert2.str();
		polyWriter->SetFileName( fileName2.c_str() );
		polyWriter->Update();
		// end debug
	}

	return returnValues;
}
template <class T,int Dim> void BME_utility<T,Dim>::vtkMeshToItkImage( const vtkSmartPointer<vtkPolyData> polyData, CharImagePointerType3D  & outImage )
{
	/************************************************************************/
	/* Create model image                                                    */
	/************************************************************************/
	typedef itk::Mesh< float, 3 >														MeshType;
	typedef itk::AutomaticTopologyMeshSource< MeshType >		MeshSourceType;
	typedef MeshSourceType::Pointer													MeshSourcePointer;

	MeshSourceType::Pointer itkMeshSource = MeshSourceType::New();
	vnl_matrix<T> inShape = vtkMeshToVnlPointMatrix( polyData );

	// add points to itkMeshSource
	MeshSourceType::PointType tempPt;

	MeshSourceType::IdentifierArrayType idArray;
	idArray.set_size( inShape.rows() );


	for ( int idxPt = 0; idxPt < inShape.rows(); ++idxPt )
	{
		for ( int idxColumn = 0; idxColumn < 3; ++idxColumn )
		{
			tempPt[idxColumn] = inShape[idxPt][idxColumn];
		}
		idArray[idxPt] = itkMeshSource->AddPoint( tempPt );
	}

	// add triangle cells to itkMeshSource
	vtkSmartPointer<vtkGenericCell> tempCell = vtkSmartPointer<vtkGenericCell>::New();

	for ( int idxCell = 0; idxCell < polyData->GetNumberOfCells(); ++idxCell )
	{
		polyData->GetCell(idxCell, tempCell);
		itkMeshSource->AddTriangle( tempCell->GetPointIds()->GetId(0), tempCell->GetPointIds()->GetId(1), tempCell->GetPointIds()->GetId(2) );
		//std::cout << tempCell->GetPointIds()->GetId(0) << "\t" << tempCell->GetPointIds()->GetId(1) << "\t" << tempCell->GetPointIds()->GetId(2) << std::endl;
	}
	MeshType::Pointer itkMesh = itkMeshSource->GetOutput();
	itkMeshSource->Update();

	// Mesh to image stuff
	typedef itk::TriangleMeshToBinaryImageFilter<MeshType,CharImageType3D> FilterType;
	FilterType::Pointer imageFilter = FilterType::New();
	imageFilter->SetInput( itkMesh );
	imageFilter->SetSpacing( outImage->GetSpacing() );
	imageFilter->SetSize( outImage->GetLargestPossibleRegion().GetSize() );
	imageFilter->SetOrigin( outImage->GetOrigin() );

	imageFilter->SetInsideValue( 255 );
	imageFilter->SetOutsideValue( 0 );
	imageFilter->Update();

	
	outImage = imageFilter->GetOutput();
}
template <class T,int Dim> vnl_matrix<T> BME_utility<T,Dim>::transformPoints( vnl_matrix_fixed<T, 4, 4> B, vnl_matrix<T> points )
{
	// create transformed points
	vnl_matrix<T> transPoints;
	transPoints.set_size( points.rows(), points.columns() );
	transPoints.fill( 0.0 );

	vnl_vector_fixed<T, 4> point;
	vnl_vector_fixed<T, 4> transPoint;


	for ( int idxPoint = 0; idxPoint < points.rows(); ++idxPoint )
	{
		// creat affine point
		for ( int idx = 0; idx < 3; ++idx )
		{
			point[idx] = points[idxPoint][idx];
		}
		point[3] = 1.0;
		// transform point using the 4x4 affine transform B
		transPoint = B * point;
		// read transformed point into output matrix
		for ( int idx = 0; idx < 3; ++idx )
		{
			transPoints[idxPoint][idx] = transPoint[idx];
		}
	}
	return transPoints;
}
template <class T,int Dim> void BME_utility<T,Dim>::sortVnlVector( const vnl_vector<T> inVector, vnl_vector<T> & sortedVector, vnl_vector<int> & idxVector )
{
	sortedVector.set_size( inVector.size() );
	idxVector.set_size( inVector.size() );
	vnl_vector<T> tempVector( inVector );

	for( int idx = 0; idx < inVector.size(); ++idx )
	{
		T maxValue = tempVector.max_value();
		sortedVector[idx] = maxValue;
		for ( int idx2 = 0; idx2 < tempVector.size(); ++idx2 )
		{
			if ( tempVector[idx2] >= maxValue*0.9999 && tempVector[idx2] <= maxValue*1.0001 )
			{
				idxVector[idx] = idx2;
				tempVector[idxVector[idx]] = -1;
				break;
			}
		}
	}
}
template <class T,int Dim> vtkSmartPointer<vtkPolyData> BME_utility<T,Dim>::getSpecificMesh( const int idx, const vtkSmartPointer<vtkPolyData> inMesh )
{
	vtkSmartPointer<vtkPolyDataConnectivityFilter> splitter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	splitter->SetInput( inMesh );
	splitter->SetExtractionModeToSpecifiedRegions();
	splitter->Update();
	splitter->AddSpecifiedRegion( idx );
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputConnection( splitter->GetOutputPort() );
	cleaner->PointMergingOff();
	vtkSmartPointer<vtkPolyData> returnMesh = cleaner->GetOutput();
	cleaner->Update();

	return returnMesh;
}
template <class T,int Dim> void BME_utility<T,Dim>::matlabSaveVnlVector( const std::string outFileName, vnl_vector<T> outVector, const std::string outVariableName )
{
	vnl_matlab_filewrite writer( outFileName.c_str() );
	writer.write( outVector, outVariableName.c_str() );
}
template <class T,int Dim> void BME_utility<T,Dim>::matlabSaveVnlMatrix( const std::string outFileName, vnl_matrix<T> outMatrix, const std::string outVariableName )
{
	vnl_matlab_filewrite writer( outFileName.c_str() );
	vnl_matrix<double> temp(outMatrix);
	writer.write( temp, outVariableName.c_str() );
}

template <class T,int Dim> void BME_utility<T,Dim>::matlabloadVnlVector( const std::string inFileName, vnl_vector<T> inVector, const std::string inVariableName )
{
	//vnl_matlab_fileread reader( inFileName.c_str() );
	//reader.read( inVector, inVariableName.c_str() );
}


template <class T,int Dim> void BME_utility<T,Dim>::matlabloadVnlMatrix( const std::string inFileName, vnl_matrix<T> inMatrix, const std::string inVariableName )
{
	//vnl_matlab_fileread reader( inFileName.c_str() );
	//vnl_matrix<double> temp(inMatrix);
	//reader.read( temp, inVariableName.c_str() );
}
#endif