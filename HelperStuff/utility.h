#ifndef _utility_h
#define _utility_h

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
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkImageReslice.h>
#include <vtkGenericCell.h>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_matlab_filewrite.h>

#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkMesh.h>
#include <itkAutomaticTopologyMeshSource.h>
#include <itkMeshSource.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkContourExtractor2DImageFilter.h>

#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <itkVTKImageToImageFilter.h>
#include <itkResampleImageFilter.h>
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#define PI 3.14159

template <class T> class utility{
private:
	int dim;
	typedef float FloatPixelType;
		
	typedef unsigned char					CharPixelType;	
	typedef itk::Image<CharPixelType, 3>	CharImageType;
	typedef itk::Image<CharPixelType, 2>	CharImageType2D;
	typedef CharImageType::Pointer			CharImagePointerType;
	typedef CharImageType2D::Pointer		CharImagePointerType2D;
	typedef itk::ImageDuplicator<CharImageType2D> DuplicatorType;  

  int m_inVal;
  int m_outVal;

public:
	utility();
	~utility();

  void SetInValue( const int val ){ m_inVal = val; }
  void SetOutValue( const int val ){ m_outVal = val; }

	vnl_vector<T> vnlShapeMatrixToVector( const vnl_matrix<T> matrix );
	vnl_matrix<T> vnlShapeVectorToMatrix( const vnl_vector<T> vector );
	vnl_matrix<T> vtkPointsToVnlPointMatrix( const vtkSmartPointer<vtkPoints> inPts );
	vtkSmartPointer<vtkPoints> vnlPointMatrixToVtkPoints( const vnl_matrix<T> inPts );
	vtkSmartPointer<vtkPolyData> vnlPointMatrixToVtkMesh( const vnl_matrix<T> S, vtkSmartPointer<vtkPolyData> templateMesh );
	vnl_matrix<T> vtkMeshToVnlPointMatrix( vtkSmartPointer<vtkPolyData> poly );
	vnl_vector<T> computeOverlap( const vtkSmartPointer<vtkPolyData> mesh, const CharImagePointerType image );
	vnl_vector<T> computeVolume( const vtkSmartPointer<vtkPolyData> mesh, const CharImagePointerType image );
	void vtkMeshToItkImage( const vtkSmartPointer<vtkPolyData> polyData, CharImagePointerType  & outImage );
	void sortVnlVector( const vnl_vector<T> inVector, vnl_vector<T> & sortedVector, vnl_vector<int> & idxVector );
	vtkSmartPointer<vtkPolyData> getSpecificMesh( const int idx, const vtkSmartPointer<vtkPolyData> inMesh );
	vnl_matrix<T> transformPoints( vnl_matrix_fixed<T, 4, 4> B, vnl_matrix<T> points );

	void matlabSaveVnlVector( const std::string outFileName, vnl_vector<T> outVector, const std::string outVariableName );
	void matlabSaveVnlMatrix( const std::string outFileName, vnl_matrix<T> outMatrix, const std::string outVariableName );

	vnl_vector<T> elementWiseMult( const vnl_vector<T> a, const vnl_vector<T> b );

	void loadAutoInitTransform( const std::string fileName, T & initScale, vnl_matrix_fixed<T, 3, 3> & initRot, vnl_vector_fixed<T, 3> & initTranslation );
  void loadScaleAnglesTranslation( const std::string fileName, T & scale, vnl_vector_fixed<T, 3> &  angles, vnl_vector_fixed<T, 3> & translation );
	vnl_vector_fixed<T, 3> getAngles( const vnl_matrix_fixed<T, 3, 3> Rot );

	void SortPolydata(vtkSmartPointer<vtkPolyData> Inpolydata , vtkSmartPointer<vtkPoints> OutpolydataPoints);
	double * GetPlaneImage(vnl_matrix_fixed<double, 4, 4> vnllandTransmatrix,double * origin,vtkSmartPointer <vtkImageData> InputImage, DuplicatorType *outImage);

	void Convert2DBitmapTo2DPoly(CharImagePointerType2D  InImg , vtkSmartPointer<vtkPolyData> polyData );
	bool Convert2DPolyTo2DBitmap(vtkSmartPointer<vtkPolyData> polyData, CharImagePointerType2D  & outImg);
	bool Convert3DPolyTo2DPoly(vtkSmartPointer<vtkPolyData> polyDataIn, vnl_matrix_fixed<double, 4, 4> &vnlPlaneTransmatrix, vtkSmartPointer<vtkPolyData> &polyDataOut );
	void Convert3DImageTo2DImage(const CharImagePointerType  InImg ,vnl_matrix_fixed<double, 4, 4> &vnlPlaneTransmatrix, CharImagePointerType2D  &outImg );
	void Convert3DPolyTo3DImage(vtkSmartPointer<vtkPolyData> polyDataIn, CharImagePointerType  &outImg );
 	void TransformPolydataTo2D( vtkSmartPointer<vtkPolyData> input, vtkSmartPointer<vtkPolyData> output,vnl_matrix_fixed<double, 4, 4> &vnlTransmatrix);
  vtkSmartPointer<vtkPolyData> transformVtkPolyData( const vtkSmartPointer<vtkPolyData>, const vtkSmartPointer<vtkTransform> );
  vtkSmartPointer<vtkMatrix4x4> convertVNLMatrixFixedToVtkMatrix( const vnl_matrix_fixed<T, 4, 4> vnlMatrix );
  vnl_matrix_fixed<T, 4, 4> convertVtkMatrixToVNLMatrixFixed( const vtkSmartPointer<vtkMatrix4x4> vtkMatrix );
  
 
  
};
template <class T> utility<T>::utility(void)
{
	dim = 3;
  m_inVal = 1;
  m_outVal = 0;
}
template <class T> utility<T>::~utility(void)
{

}

template <class T> vnl_matrix<T> utility<T>::vnlShapeVectorToMatrix( const vnl_vector<T> vector )
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
template <class T> vnl_vector<T> utility<T>::vnlShapeMatrixToVector( const vnl_matrix<T> matrix )
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
template <class T> vnl_matrix<T>  utility<T>::vtkPointsToVnlPointMatrix( const vtkSmartPointer<vtkPoints> inPts )
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
template <class T> vtkSmartPointer<vtkPoints>  utility<T>::vnlPointMatrixToVtkPoints( const vnl_matrix<T> inPts )
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
template <class T> vtkSmartPointer<vtkPolyData> utility<T>::vnlPointMatrixToVtkMesh( const vnl_matrix<T> S, vtkSmartPointer<vtkPolyData> templateMesh )
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
template <class T> vnl_matrix<T> utility<T>::vtkMeshToVnlPointMatrix( vtkSmartPointer<vtkPolyData> poly )
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
template <class T> vnl_vector<T> utility<T>::computeOverlap( const vtkSmartPointer<vtkPolyData> mesh, const CharImagePointerType image )
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

	CharImageType::Pointer meshImage = CharImageType::New();
	meshImage->SetRegions( image->GetLargestPossibleRegion() );
	meshImage->SetSpacing( image->GetSpacing() );
	meshImage->SetOrigin( image->GetOrigin() );
	meshImage->Allocate();
	meshImage->FillBuffer( 0 );

	vtkMeshToItkImage( cleaner->GetOutput(), meshImage );

	typedef itk::ImageRegionConstIterator<CharImageType> CharImageIteratorType;
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
	//typedef itk::ImageFileWriter<CharImageType> CharImageWriterType;
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
template <class T> vnl_vector<T> utility<T>::computeVolume( const vtkSmartPointer<vtkPolyData> mesh, const CharImagePointerType image )
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

	CharImageType::Pointer meshImage = CharImageType::New();
	meshImage->SetSpacing( image->GetSpacing() );
	
	// find size for mesh image
	CharImageType::RegionType region;
	CharImageType::SizeType size;
	CharImageType::IndexType index;
	CharImageType::PointType origin;
	CharImageType::SpacingType spacing = image->GetSpacing(); 
	index.Fill( 0 );
	T temp[6];
	mesh->GetBounds(temp);
	/*std::cout << "bounds: x1" << temp[0] << "\t x2" <<  temp[1] 
	<< "\t y1" <<  temp[2] << "\t y2" <<  temp[3] 
	<< "\t z1" <<  temp[4] << "\t z2" <<  temp[5]	<<  std::endl;*/
	
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
	
	//// debug
	//typedef itk::ImageFileWriter<CharImageType> CharImageWriterType;
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
	
	typedef itk::ImageRegionConstIterator<CharImageType> CharImageIteratorType;
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

		//// debug
		//std::string fileName;
		//std::stringstream convert;
		//writer->SetInput( meshImage );
		//
		//convert << "meshImage" << idx << ".mhd";
		//fileName = convert.str();
		//
		//writer->SetFileName( fileName.c_str() );
		//writer->Update();

		//// polyWriter->SetInput( cleaner->GetOutput() );
		//std::string fileName2;
		//std::stringstream convert2;
		//convert2 << "mesh" << idx << ".vtk";
		//fileName2 = convert2.str();
		//polyWriter->SetFileName( fileName2.c_str() );
		//polyWriter->Update();
		//// end debug
	}

	return returnValues;
}
template <class T> void utility<T>::vtkMeshToItkImage( const vtkSmartPointer<vtkPolyData> polyData, CharImagePointerType  & outImage )
{
  // create vtk image which contains the binary volume
  vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    
  double bounds[6];
  polyData->GetBounds(bounds);
  double spacing[3]; // desired volume spacing
  spacing[0] = 0.5;
  spacing[1] = 0.5;
  spacing[2] = 0.5;
  whiteImage->SetSpacing(spacing);

  // compute dimensions
  int dim[3];
  for (int i = 0; i < 3; i++)
  {
    dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
  }
  whiteImage->SetDimensions(dim);
  whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

  double origin[3];
  origin[0] = bounds[0] + spacing[0] / 2;
  origin[1] = bounds[2] + spacing[1] / 2;
  origin[2] = bounds[4] + spacing[2] / 2;
  whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
  whiteImage->SetScalarTypeToUnsignedChar();
  whiteImage->AllocateScalars();
#else
  whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif

  // fill the image with foreground voxels:
  unsigned char inval = m_inVal;
  unsigned char outval = m_outVal;
  vtkIdType count = whiteImage->GetNumberOfPoints();
  for (vtkIdType i = 0; i < count; ++i)
  {
    whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
  }

  // polygonal data --> image stencil:
  vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
  pol2stenc->SetInput( polyData );
#else
  pol2stenc->SetInputData( polyData );
#endif
  pol2stenc->SetOutputOrigin(origin);
  pol2stenc->SetOutputSpacing(spacing);
  pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  pol2stenc->Update();

  // cut the corresponding white image and set the background:
  vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
  imgstenc->SetInput(whiteImage);
  imgstenc->SetStencil(pol2stenc->GetOutput());
#else
  imgstenc->SetInputData(whiteImage);
  imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
  imgstenc->ReverseStencilOff();
  imgstenc->SetBackgroundValue(outval);
  imgstenc->Update();

  // convert to ITK image
  typedef itk::VTKImageToImageFilter<CharImageType> VtkToItkImageConverterType;
  VtkToItkImageConverterType::Pointer vtkToItkConverter = VtkToItkImageConverterType::New();
  vtkToItkConverter->SetInput( imgstenc->GetOutput() );
  try
  {
    vtkToItkConverter->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Vtk to ITK image conversion failed." << std::endl;
    std::cerr << excp << std::endl;
  }

  // resample binary image
  typedef itk::ResampleImageFilter<CharImageType, CharImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( vtkToItkConverter->GetOutput() );
  resampler->SetReferenceImage( outImage );
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
  outImage = resampler->GetOutput();
}
template <class T> vnl_matrix<T> utility<T>::transformPoints( vnl_matrix_fixed<T, 4, 4> B, vnl_matrix<T> points )
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
template <class T> void utility<T>::sortVnlVector( const vnl_vector<T> inVector, vnl_vector<T> & sortedVector, vnl_vector<int> & idxVector )
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
template <class T> vtkSmartPointer<vtkPolyData> utility<T>::getSpecificMesh( const int idx, const vtkSmartPointer<vtkPolyData> inMesh )
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
template <class T> void utility<T>::matlabSaveVnlVector( const std::string outFileName, vnl_vector<T> outVector, const std::string outVariableName )
{
	vnl_matlab_filewrite writer( outFileName.c_str() );
	vnl_vector<double> temp(outVector.size());
	for ( int idx = 0; idx < temp.size(); ++idx )
	{
		temp[idx] = outVector[idx];
	}
	writer.write( temp, outVariableName.c_str() );
}
template <class T> void utility<T>::matlabSaveVnlMatrix( const std::string outFileName, vnl_matrix<T> outMatrix, const std::string outVariableName )
{
	vnl_matlab_filewrite writer( outFileName.c_str() );
	vnl_matrix<double> temp(outMatrix);
	writer.write( temp, outVariableName.c_str() );
}
template <class T> vnl_vector<T> utility<T>::elementWiseMult( const vnl_vector<T> a, const vnl_vector<T> b )
{
	vnl_vector<T> outVector( a.size() );
	for ( int idx = 0; idx < a.size(); ++idx )
	{
		outVector[idx] = a[idx] * b[idx];
	}
	return outVector;
}
template <class T> void utility<T>::loadAutoInitTransform( const std::string fileName, T & initScale, vnl_matrix_fixed<T, 3, 3> & initRot, vnl_vector_fixed<T, 3> & initTranslation )
{
  std::ifstream transRead ( fileName.c_str() , ifstream::in );
  if ( transRead.good() )
  {
    initScale = 1.0;
    transRead >> initRot[0][0]; transRead >> initRot[0][1]; transRead >> initRot[0][2];
    transRead >> initRot[1][0]; transRead >> initRot[1][1]; transRead >> initRot[1][2];
    transRead >> initRot[2][0]; transRead >> initRot[2][1]; transRead >> initRot[2][2];
    transRead >> initTranslation[0]; transRead >> initTranslation[1]; transRead >> initTranslation[2];
    transRead.close();
  } 
  else
  {
    std::cout << "Init transform could not be read." << std::endl;
  }
}
template <class T> void utility<T>::loadScaleAnglesTranslation( const std::string fileName, T & scale, vnl_vector_fixed<T, 3> &  angles, vnl_vector_fixed<T, 3> & translation )
{
  std::ifstream transRead ( fileName.c_str() , ifstream::in );
  if ( transRead.good() )
  {
    transRead >> scale;
    transRead >> angles[0]; transRead >> angles[1]; transRead >> angles[2];
    transRead >> translation[0]; transRead >> translation[1]; transRead >> translation[2];
    transRead.close();
  } 
  else
  {
    std::cout << "Scale, angles, and translation could not be read." << std::endl;
  }
}

template <class T> vnl_vector_fixed<T, 3> utility<T>::getAngles( const vnl_matrix_fixed<T, 3, 3> Rot )
{
  itk::Euler3DTransform<double>::Pointer itkRot = itk::Euler3DTransform<double>::New();
  itk::Euler3DTransform<double>::MatrixType itkMatrix;
  for (int idxC = 0; idxC<3; ++idxC)
  {
    for ( int idxR = 0; idxR < 3; ++idxR )
    {
      itkMatrix[idxR][idxC] = Rot[idxR][idxC];
    }
  }

  itkRot->SetMatrix( itkMatrix );
  vnl_vector_fixed<T, 3> angles;
  angles[0] = itkRot->GetAngleX()*180.0/PI; angles[1] = itkRot->GetAngleY()*180.0/PI; angles[2] = itkRot->GetAngleZ()*180.0/PI;
  return angles;
}

template <class T> void utility<T>::SortPolydata(vtkSmartPointer<vtkPolyData> Inpolydata , vtkSmartPointer<vtkPoints> pointsSorted)
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

template <class T> double * utility<T>::GetPlaneImage(vnl_matrix_fixed<double, 4, 4> vnlPlaneTransmatrix,double * origin,vtkSmartPointer <vtkImageData> InputImage,DuplicatorType *outImage)
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

	//BMET_ImageReslice->SetOutputSpacing(1,1,1);
	//BMET_ImageReslice->SetOutputSpacing(InputImage->GetSpacing());

	//BMET_ImageReslice->InterpolateOff ();
	BMET_ImageReslice->AutoCropOutputOn ();
	BMET_ImageReslice->Modified();
	BMET_ImageReslice->SetInterpolationModeToLinear();
	BMET_ImageReslice->TransformInputSamplingOn();

#if VTK_MAJOR_VERSION <= 5
	BMET_ImageReslice->SetInput(InputImage);
#else
	BMET_ImageReslice->SetInputData(imageReader->GetOutput());
#endif

	BMET_ImageReslice->GetResliceAxes()->DeepCopy( Planematrix); 
	BMET_ImageReslice->SetResliceAxesOrigin(origin);  
	BMET_ImageReslice->Modified();
	BMET_ImageReslice->Update();
	BMET_ImageReslice->GetOutput()->GetBounds(ImageBounds);
	BMET_ImageReslice->GetOutput()->GetSpacing(ImageSpacing);
	ImageSpacing[3]=0;
	typedef itk::VTKImageToImageFilter<CharImageType2D>ImageConverterVTKToITK;
	ImageConverterVTKToITK::Pointer imageConverterResliceVTKToITK2D = ImageConverterVTKToITK::New();
	imageConverterResliceVTKToITK2D->SetInput(BMET_ImageReslice->GetOutput());
	imageConverterResliceVTKToITK2D->Update();
	outImage->SetInputImage(imageConverterResliceVTKToITK2D->GetOutput());
	outImage->Update();

	return (double *) ImageBounds;
};

template <class T> void utility<T>::Convert2DBitmapTo2DPoly(CharImagePointerType2D  inImage , vtkSmartPointer<vtkPolyData> polyData )
{

	typedef unsigned char				PixelType;
	typedef itk::Image<PixelType, 2>	ImageType;
	typedef itk::ImageFileReader<ImageType> ImageReaderType;

	ImageReaderType::Pointer reader = ImageReaderType::New();

	ImageType::PointType origin=inImage->GetOrigin(); ;
	ImageType::SpacingType spacing = inImage->GetSpacing(); 

	/************************************************************************/
	/* create contour                                                       */
	/************************************************************************/
	typedef itk::ImageToVTKImageFilter<ImageType> ImageConverterITKToVTK2D;
	ImageConverterITKToVTK2D::Pointer imageConverterITKToVTK2D = ImageConverterITKToVTK2D::New();
	vtkSmartPointer<vtkPolyData> vertex = vtkSmartPointer<vtkPolyData>::New();

	typedef itk::ContourExtractor2DImageFilter <ImageType>   ContourExtractor2DImageFilterType;
	ContourExtractor2DImageFilterType::Pointer contourExtractor2DImageFilter = ContourExtractor2DImageFilterType::New();
	contourExtractor2DImageFilter->SetInput( inImage );
	contourExtractor2DImageFilter->SetContourValue(0);
	contourExtractor2DImageFilter->Update();
	vtkSmartPointer<vtkPoints> curvepoints;
	curvepoints = vtkSmartPointer<vtkPoints>::New();
	vertex->SetPoints(curvepoints);

	int i = 0;
	int maxnrpoints=0;
	int maxnrpointsIndex=0;
	int  numPts=0;
	for(unsigned int i = 0; i < contourExtractor2DImageFilter->GetNumberOfOutputs(); ++i)
	{
		numPts = contourExtractor2DImageFilter->GetOutput(i)->GetVertexList()->Size();
		if(maxnrpoints<numPts)
		{
			//std::cout << "Max Contour " << i << ": "<< numPts << std::endl;
			maxnrpoints=numPts;
			maxnrpointsIndex=i;
		}
		//else
		//std::cout << "Contour " << i << ": "<< numPts << std::endl;
	}
	vtkUnsignedCharArray *curvecolor=vtkUnsignedCharArray::New();
	curvecolor->SetNumberOfComponents(3);
	curvecolor->SetNumberOfTuples(maxnrpoints);
	vtkCellArray   *Curvelines =	vtkCellArray::New();
	double ptr[3];
	double ptr1[3];
	int idx=0;
	curvepoints->SetNumberOfPoints(maxnrpoints);


	// for(unsigned int i = 0; i < contourExtractor2DImageFilter->GetNumberOfOutputs(); ++i)
	if(maxnrpoints)
	{
		i= maxnrpointsIndex;
		//std::cout << "Contour " << i << ": " << std::endl;

		ContourExtractor2DImageFilterType::VertexListType::ConstIterator vertexIterator = 
			contourExtractor2DImageFilter->GetOutput(i)->GetVertexList()->Begin();
		while(vertexIterator != contourExtractor2DImageFilter->GetOutput(i)->GetVertexList()->End())
		{
			// std::cout << vertexIterator->Value() << std::endl;
			ptr1[0]=vertexIterator->Value()[0];
			ptr1[1]=vertexIterator->Value()[1];
			ptr[0]=vertexIterator->Value()[0]*spacing[0]+origin[0];
			ptr[1]=vertexIterator->Value()[1]*spacing[1]+origin[1];
			ptr[2]=0;

			//std::cout << ptr[0]  << ": " << ptr[1]<< std::endl;
			++vertexIterator;

			curvepoints->InsertPoint(idx,ptr);
			curvecolor->SetTuple3(idx ,255,138,126);
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
	polyData->DeepCopy( vertex);
}

template <class T> bool utility<T>::Convert2DPolyTo2DBitmap(vtkSmartPointer<vtkPolyData> polyData,CharImagePointerType2D  & outImg)
{

	CharImageType2D::RegionType iRegion2D = outImg->GetLargestPossibleRegion();
	CharImageType2D::SizeType   iSize2D = iRegion2D.GetSize();

	// create vtk image which contains the binary volume
	vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    
	double bounds[6];
	polyData->GetBounds(bounds);
	double spacing[3]; // desired volume spacing
	spacing[0] = 0.5;
	spacing[1] = 0.5;
	spacing[2] = 0.5;
	whiteImage->SetSpacing(spacing);

	// compute dimensions
	int dim[3];
	for (int i = 0; i < 3; i++)
	{
		dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
	}
	int numPts=polyData->GetPoints()->GetNumberOfPoints() ;
	if(!numPts)
	{
		dim[0]=iSize2D[0];
		dim[1]=iSize2D[1];
	}

	dim[2]=2;
	whiteImage->SetDimensions(dim);
	whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

	double origin[3];
	origin[0] = bounds[0] + spacing[0] / 2;
	origin[1] = bounds[2] + spacing[1] / 2;
	origin[2] = bounds[4] + spacing[2] / 2;
	whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
	whiteImage->SetScalarTypeToUnsignedChar();
	whiteImage->AllocateScalars();
#else
	whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif

	// fill the image with foreground voxels:
	unsigned char inval = m_inVal = 1;
	unsigned char outval = m_outVal =0;
	vtkIdType count = whiteImage->GetNumberOfPoints();
	for (vtkIdType i = 0; i < count; ++i)
	{
		whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
	}
	vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();

	if(numPts)
	{

		// polygonal data --> image stencil:
#if VTK_MAJOR_VERSION <= 5
		pol2stenc->SetInput( polyData );
#else
		pol2stenc->SetInputData( polyData );
#endif
		pol2stenc->SetOutputOrigin(origin);
		pol2stenc->SetOutputSpacing(spacing);
		pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
		pol2stenc->Update();
	}
	// cut the corresponding white image and set the background:
	vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
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


	// convert to ITK image
	typedef itk::VTKImageToImageFilter<CharImageType> VtkToItkImageConverterType;
	VtkToItkImageConverterType::Pointer vtkToItkConverter = VtkToItkImageConverterType::New();
	vtkToItkConverter->SetInput( imgstenc->GetOutput() );
	try
	{
		vtkToItkConverter->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Vtk to ITK image conversion failed." << std::endl;
		std::cerr << excp << std::endl;
	}

	typedef itk::ExtractImageFilter<CharImageType, CharImageType2D> ImageFilter3Dto2D; //Get Slice out of 3D stack

	CharImageType::RegionType iRegion3D = vtkToItkConverter->GetOutput()->GetLargestPossibleRegion();
	CharImageType::SizeType   iSize3D = iRegion3D.GetSize();
	CharImageType::RegionType extractRegion = iRegion3D;
	CharImageType::SizeType   extractSize = iSize3D;
	CharImageType::IndexType  extractIndex = extractRegion.GetIndex();

	ImageFilter3Dto2D::Pointer imageSlicers3Dto2D = ImageFilter3Dto2D::New();
	extractSize[2] = 0;
	extractIndex[0] = extractIndex[1] = 0;
	extractIndex[2]  = 0;
	extractRegion.SetSize(extractSize);
	extractRegion.SetIndex(extractIndex);


	imageSlicers3Dto2D->SetExtractionRegion(extractRegion);
	imageSlicers3Dto2D->SetDirectionCollapseToIdentity();
	imageSlicers3Dto2D->SetInput(vtkToItkConverter->GetOutput());
	imageSlicers3Dto2D->Update();
	//vtkSmartPointer<vtkMetaImageWriter> imageWriter =vtkSmartPointer<vtkMetaImageWriter>::New();
	//imageWriter->SetFileName("F:/Data/userDis1.mhd");
	//imageWriter->SetInput(imgstenc->GetOutput());
	//imageWriter->Write();
	//WriteMHD("F:/Data/userDis12.mhd",imageSlicers3Dto2D->GetOutput());


	// resample binary image
	typedef itk::ResampleImageFilter<CharImageType2D, CharImageType2D> ResampleFilterType;
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInput( imageSlicers3Dto2D->GetOutput() );
	resampler->SetReferenceImage( outImg );
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
	outImg = resampler->GetOutput();
	if(numPts>2)
		return true;
	else 
		return false;
}

template <class T> bool utility<T>::Convert3DPolyTo2DPoly(vtkSmartPointer<vtkPolyData> polyDataIn,vnl_matrix_fixed<double, 4, 4> &vnlPlaneTransmatrix, vtkSmartPointer<vtkPolyData> &polyDataOut )
{
	double origin[3];
	double normal[3];
	for (int k=0;k<3;k++)
	{
		origin[k]=vnlPlaneTransmatrix[k][3];
		normal[k]=vnlPlaneTransmatrix[k][2];
	}

	vtkSmartPointer<vtkPlane> polyCutPlane =	vtkSmartPointer<vtkPlane>::New();

	polyCutPlane->SetOrigin( origin );
	polyCutPlane->SetNormal( normal );
	polyCutPlane->Modified();

	vtkSmartPointer<vtkCutter> cutter =	vtkSmartPointer<vtkCutter>::New();
	cutter->SetCutFunction( polyCutPlane );
	cutter->SetInput( polyDataIn);
	cutter->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> splitter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	splitter->SetInputConnection( cutter->GetOutputPort() );
	splitter->SetExtractionModeToSpecifiedRegions();
	splitter->AddSpecifiedRegion( 0 );
	splitter->Modified();
	splitter->Update();


	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputConnection( splitter->GetOutputPort() );
	cleaner->PointMergingOff();
	cleaner->Update();


	int numCavities = splitter->GetNumberOfExtractedRegions();
	short MaxPoints=0;
	short selectedRegion=0;
	vtkSmartPointer<vtkPolyData> intersect = vtkSmartPointer<vtkPolyData>::New();
	intersect->DeepCopy(cleaner->GetOutput());

	vtkSmartPointer<vtkPoints> points	= intersect->GetPoints();
	if(points)	MaxPoints = points->GetNumberOfPoints();

	if(numCavities>1)
	{

		for( int idx = 1; idx < numCavities; ++idx )
		{
			// do other cavity inside of the loop


			vtkSmartPointer<vtkPoints> points1	= intersect->GetPoints();
			short nrPoints1=0;
			if(points1)
			{
				nrPoints1=points1->GetNumberOfPoints();
			}


			vtkSmartPointer<vtkPolyDataConnectivityFilter> BME_Spliter2 = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
			BME_Spliter2->SetInputConnection( cutter->GetOutputPort() );
			BME_Spliter2->SetExtractionModeToSpecifiedRegions();
			BME_Spliter2->AddSpecifiedRegion( idx );
			BME_Spliter2->Update();
			vtkSmartPointer<vtkCleanPolyData>	BME_Cleaner2	= vtkSmartPointer<vtkCleanPolyData>::New() ;
			BME_Cleaner2->SetInputConnection( BME_Spliter2->GetOutputPort() );
			BME_Cleaner2->Modified();
			BME_Cleaner2->Update();
			vtkSmartPointer<vtkPolyData> intersect2 = BME_Cleaner2->GetOutput();
			vtkSmartPointer<vtkPoints> points2	= intersect2->GetPoints();
			short nrPoints2=0;
			if(points2)
			{
				nrPoints2=points2->GetNumberOfPoints();
			}
			if(nrPoints1<nrPoints2)
				intersect->DeepCopy(intersect2);

			BME_Cleaner2->RemoveAllInputs();
			BME_Cleaner2->RemoveAllObservers();
			BME_Spliter2->RemoveAllInputs();
			BME_Spliter2->RemoveAllObservers();

		}
	}

	vtkSmartPointer<vtkPolyData> polyDataOut2D = vtkSmartPointer<vtkPolyData>::New();	
	vtkSmartPointer<vtkPoints> sortedPoints2D = vtkSmartPointer<vtkPoints>::New();
	sortedPoints2D->SetNumberOfPoints( intersect->GetNumberOfPoints() ); // Sort points
	SortPolydata(intersect ,sortedPoints2D ); // Sort points
	intersect->SetPoints(sortedPoints2D);

	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
	for( int xx=0;xx<4;xx++)
		for( int yy=0;yy<4;yy++)
			matrix->SetElement(xx,yy,vnlPlaneTransmatrix[xx][yy]);
	transform->SetMatrix( matrix );
	transform->Update();
	
	vtkSmartPointer<vtkPoints>		Points2D		= vtkSmartPointer<vtkPoints>::New() ;
	double * point2D;
	double nrpoints =intersect->GetNumberOfPoints() ;
	Points2D->SetNumberOfPoints( nrpoints);
	vtkCellArray   *Curvelines =	vtkCellArray::New();
	for ( int idx = 0; idx <(nrpoints); ++idx )
	{
		point2D=transform->GetInverse()->TransformDoublePoint(intersect->GetPoint(idx));
		point2D[2]=0;
		Points2D->SetPoint( idx, point2D );
		if ( idx > 0 )
		{
			vtkIdType	cell[2]	= {idx-1,idx};
			Curvelines->InsertNextCell((vtkIdType)2,(vtkIdType*)	cell);
		}
	}

	polyDataOut2D->SetPoints(Points2D);
	polyDataOut2D->SetLines(Curvelines);

	polyDataOut->DeepCopy( polyDataOut2D);
	if(nrpoints>2)
		return true;
	else 
		return false;

}

template <class T> void utility<T>::Convert3DImageTo2DImage(const CharImagePointerType  InImg ,vnl_matrix_fixed<double, 4, 4> &vnlPlaneTransmatrix,  CharImagePointerType2D  &outImg )
{

	double origin[3];
	double normal[3];
	for (int k=0;k<3;k++)
	{
		origin[k]=vnlPlaneTransmatrix[k][3];
		normal[k]=vnlPlaneTransmatrix[k][2];
	}

	
	typedef itk::ImageToVTKImageFilter<CharImageType> ImageConverterITKToVTK;
	typedef itk::ImageDuplicator<CharImageType2D> DuplicatorType;  
	
	ImageConverterITKToVTK::Pointer imageConverterITKToVTK3D = ImageConverterITKToVTK::New();
	imageConverterITKToVTK3D->SetInput(InImg);
	imageConverterITKToVTK3D->Update();
	DuplicatorType::Pointer Image2D = DuplicatorType::New();;

	GetPlaneImage(vnlPlaneTransmatrix, origin,imageConverterITKToVTK3D->GetOutput(),Image2D);
	// resample binary image
	typedef itk::ResampleImageFilter<CharImageType2D, CharImageType2D> ResampleFilterType;
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInput( Image2D->GetOutput() );
	resampler->SetReferenceImage( outImg );
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
	outImg = resampler->GetOutput();

}

template <class T> void utility<T>::Convert3DPolyTo3DImage(vtkSmartPointer<vtkPolyData> polyData, CharImagePointerType  &outImage )
{

	// create vtk image which contains the binary volume
	vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    
	double bounds[6];
	polyData->GetBounds(bounds);
	double spacing[3]; // desired volume spacing
	spacing[0] = 0.5;
	spacing[1] = 0.5;
	spacing[2] = 0.5;
	whiteImage->SetSpacing(spacing);

	// compute dimensions
	int dim[3];
	for (int i = 0; i < 3; i++)
	{
		dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
	}
	whiteImage->SetDimensions(dim);
	whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

	double origin[3];
	origin[0] = bounds[0] + spacing[0] / 2;
	origin[1] = bounds[2] + spacing[1] / 2;
	origin[2] = bounds[4] + spacing[2] / 2;
	whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
	whiteImage->SetScalarTypeToUnsignedChar();
	whiteImage->AllocateScalars();
#else
	whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif

	// fill the image with foreground voxels:
	unsigned char inval = m_inVal;
	unsigned char outval = m_outVal;
	vtkIdType count = whiteImage->GetNumberOfPoints();
	for (vtkIdType i = 0; i < count; ++i)
	{
		whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
	}

	// polygonal data --> image stencil:
	vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
	pol2stenc->SetInput( polyData );
#else
	pol2stenc->SetInputData( polyData );
#endif
	pol2stenc->SetOutputOrigin(origin);
	pol2stenc->SetOutputSpacing(spacing);
	pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
	pol2stenc->Update();

	// cut the corresponding white image and set the background:
	vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
	imgstenc->SetInput(whiteImage);
	imgstenc->SetStencil(pol2stenc->GetOutput());
#else
	imgstenc->SetInputData(whiteImage);
	imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
	imgstenc->ReverseStencilOff();
	imgstenc->SetBackgroundValue(outval);
	imgstenc->Update();

	// convert to ITK image
	typedef itk::VTKImageToImageFilter<CharImageType> VtkToItkImageConverterType;
	VtkToItkImageConverterType::Pointer vtkToItkConverter = VtkToItkImageConverterType::New();
	vtkToItkConverter->SetInput( imgstenc->GetOutput() );
	try
	{
		vtkToItkConverter->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Vtk to ITK image conversion failed." << std::endl;
		std::cerr << excp << std::endl;
	}

  // check if outImage has already been allocated
  CharImageType::SizeType outSize = outImage->GetLargestPossibleRegion().GetSize();
  int accuSize = 0;
  for ( int itr = 0; itr < 3; itr++ )
  {
    accuSize += outSize[itr];
  }

  typedef itk::ResampleImageFilter<CharImageType, CharImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( vtkToItkConverter->GetOutput() );


	if ( accuSize > 0 )
  {
	  // resample binary image if it has already been allocated
		resampler->SetReferenceImage( outImage );  
  } 
  else
  {
    resampler->SetReferenceImage( vtkToItkConverter->GetOutput() );
  }

  resampler->UseReferenceImageOn();
  try
  {
    resampler->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "utility::getSpecificMesh:resampling image failed." << std::endl;
    std::cerr << excp << std::endl;
  }
  outImage = resampler->GetOutput();
}

template <class T> void utility<T>::TransformPolydataTo2D( vtkSmartPointer<vtkPolyData> input, vtkSmartPointer<vtkPolyData> output,vnl_matrix_fixed<double, 4, 4> &vnlTransmatrix)
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

template <class T> vtkSmartPointer<vtkPolyData> utility<T>::transformVtkPolyData( const vtkSmartPointer<vtkPolyData> mesh, const vtkSmartPointer<vtkTransform> transform )
{
  vtkSmartPointer<vtkTransformPolyDataFilter> tranfsformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  tranfsformer->SetTransform( transform );
  tranfsformer->SetInput( mesh );
  tranfsformer->Update();

  vtkSmartPointer<vtkPolyData> outMesh = vtkSmartPointer<vtkPolyData>::New();
  outMesh = tranfsformer->GetOutput();
  return outMesh;
}

template <class T> vtkSmartPointer<vtkMatrix4x4> utility<T>::convertVNLMatrixFixedToVtkMatrix( const vnl_matrix_fixed<T, 4, 4> vnlMatrix )
{
  vtkSmartPointer<vtkMatrix4x4> vtkMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  for ( int rowIt = 0; rowIt < 4; ++rowIt )
  {
    for ( int colIt = 0; colIt < 4; ++colIt )
    {
      vtkMatrix->SetElement( rowIt, colIt, vnlMatrix[rowIt][colIt] ); 
    }
  }

  return vtkMatrix;
}

template <class T> vnl_matrix_fixed<T, 4, 4> utility<T>::convertVtkMatrixToVNLMatrixFixed( const vtkSmartPointer<vtkMatrix4x4> vtkMatrix )
{
  vnl_matrix_fixed<T, 4, 4> vnlMatrix;
  for ( int rowIt = 0; rowIt < 4; ++rowIt )
  {
    for ( int colIt = 0; colIt < 4; ++colIt )
    {
      vnlMatrix[rowIt][colIt] = vtkMatrix->GetElement( rowIt, colIt );
    }
  }
  return vnlMatrix;
}

#endif