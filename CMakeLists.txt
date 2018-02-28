cmake_minimum_required (VERSION 2.6)
project (BB3D)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS -std=c++0x)

set (OCE_INCLUDE_DIR "" CACHE PATH "Path to Open CASCADE include directory.")
set (OCE_LIB_DIR "" CACHE PATH "Path to Open CASCADE library directory.")
include_directories(include)
include_directories(${OCE_INCLUDE_DIR})
link_directories(${OCE_LIB_DIR})
file(GLOB SOURCES "src/*.cpp")
add_executable(make_blade ${SOURCES})
target_link_libraries(make_blade TKBinL TKBin TKBinTObj TKBinXCAF TKBool TKBO TKBRep TKCAF TKCDF TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKIVtk TKLCAF TKMath TKMesh TKMeshVS TKOffset TKOpenGl TKPrim TKService TKShHealing TKStdL TKStd TKSTEP209 TKSTEPAttr TKSTEPBase TKSTEP TKSTL TKTObj TKTopAlgo TKV3d TKVCAF TKVRML TKXCAF TKXDEIGES TKXDESTEP TKXMesh TKXmlL TKXml TKXmlTObj TKXmlXCAF TKXSBase)
install(TARGETS make_blade DESTINATION bin)
#target_link_libraries(TKernel "${OCE_LIB}")
#target_link_libraries(make_blade TKBinL TKBin TKBinTObj TKBinXCAF TKBool TKBO TKBRep TKCAF TKCDF TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKIVtk TKLCAF TKMath TKMesh TKMeshVS TKOffset TKOpenGl TKPrim TKService TKShHealing TKStdL TKStd TKSTEP209 TKSTEPAttr TKSTEPBase TKSTEP TKSTL TKTObj TKTopAlgo TKV3d TKVCAF TKVRML TKXCAF TKXDEIGES TKXDESTEP TKXMesh TKXmlL TKXml TKXmlTObj TKXmlXCAF TKXSBase)


#include_directories(${OCE_INCLUDE_DIR})
#include_directories(include)
#link_directories(${OCE_LIB_DIR})
#
#
#add_executable(BB3D ${SOURCES})
#target_link_libraries(BB3D ${OCE_LIBS})
