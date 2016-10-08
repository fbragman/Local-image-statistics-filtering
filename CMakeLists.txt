cmake_minimum_required (VERSION 2.6)

project(trunk)

add_library(Image Image.cpp)
add_library(Stats Stats.cpp)

FIND_PATH(nifticlib nifticlib)
INCLUDE_DIRECTORIES(${nifticlib}/include)
LINK_DIRECTORIES(${nifticlib}/lib)

add_executable(VolumeQuery VolumeQuery.cpp)
add_executable(VolumeMedian VolumeMedian.cpp)

target_link_libraries(VolumeQuery Image Stats niftiio znz z)
target_link_libraries(VolumeMedian Image Stats niftiio znz z)







