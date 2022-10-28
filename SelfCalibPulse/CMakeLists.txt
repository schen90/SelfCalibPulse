
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(SelfCalib)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
#option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
#if(WITH_GEANT4_UIVIS)
#  find_package(Geant4 REQUIRED ui_all vis_all)
#else()
#  find_package(Geant4 REQUIRED)
#endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
#include(${Geant4_USE_FILE})
include(FindROOT.cmake)
#find_package(ROOT)
set(ROOT_REQUIRED_COMPONENTS Core Hist Tree Rint Matrix MathCore Thread)
find_package(ROOT REQUIRED ${ROOT_REQUIRED_COMPONENTS})

include_directories(${PROJECT_SOURCE_DIR}/include
                    ${ROOT_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/src/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(SelfCalib src/SelfCalib.cc ${sources} ${headers})
target_link_libraries(SelfCalib ${ROOT_LIBRARIES})

SET ( CMAKE_CXX_FLAGS "-std=c++11 -mavx -pthread" CACHE STRING "compile flags" FORCE)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SelfCalib DESTINATION bin)


