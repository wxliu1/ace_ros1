cmake_minimum_required(VERSION 3.15)
project(test_sync_cb)

set(CMAKE_BUILD_TYPE "Release") #lwx add

# for gdb debug 2023-3-27.
# SET(CMAKE_BUILD_TYPE "Debug")
# SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
# the end.

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

find_package(catkin REQUIRED COMPONENTS 
  roscpp
  #tf2_ros
  #tf2_eigen
  #cv_bridge
  #nav_msgs
  sensor_msgs
  #message_generation
  #message_runtime
  std_msgs
  roslib
  message_filters
)

catkin_package(
  #INCLUDE_DIRS include
#  LIBRARIES binocular_sensor
  CATKIN_DEPENDS roscpp std_msgs rosmsg
#  DEPENDS system_lib
)

include_directories(
  # lib
  # ${catkin_LIB_DIRS}
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  "./"
  "/usr/local/include"
)
MESSAGE(STATUS "aaaa This is sophus library dir " ${Sophus_LIBRARIES})

add_executable(test_sync_cb
  ${CMAKE_CURRENT_SOURCE_DIR}/src/test_sync_cb.cpp
)

target_link_libraries(test_sync_cb 
    ${catkin_LIBRARIES}
)