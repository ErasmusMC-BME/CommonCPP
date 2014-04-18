############################################################################
# <summary> CMakeLists.txt for OpenCV and CUDA. </summary>
# <date>    2012-02-07          </date>
# <author>  Quan Tran Minh. edit by Johannes Kast, Michael Sarahan </author>
# <email>   quantm@unist.ac.kr  kast.jo@googlemail.com msarahan@gmail.com</email>
############################################################################

# collect source files

file( GLOB  HW5_hdr *.hpp *.h )
file( GLOB  HW5_cu  *.cu)

SET (HW5_files main.cu student.cu reference_calc.cpp)


SET_SOURCE_FILES_PROPERTIES( ${HW5_hdr} PROPERTIES HEADER_FILE_ONLY TRUE )
SOURCE_GROUP( "headers" FILES ${HW5_hdr} )
SOURCE_GROUP( "Cuda sources" FILES ${HW5_cu} )
SOURCE_GROUP( "sources" FILES ${HW5_files} )

SET(CUDA_NVCC_FLAGS ${MY_NVCC_OPTIONS})

CUDA_ADD_EXECUTABLE(HW5 ${HW5_files} ${HW5_hdr} ${HW5_cu})




