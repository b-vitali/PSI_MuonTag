
#----------------------------------------------------------------------------
# Setup project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(Entrance)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or cmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else()
    find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup inclue directory for this project#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#ROOT_GENERATE_DICTIONARY(G__Event Event.h LINKDEF EventLinkDef.h)

include_directories(
    ${PROJECT_SOURCE_DIR}/SiPM/inc
    ${PROJECT_SOURCE_DIR}/SiPM/cc
    ${PROJECT_SOURCE_DIR}/inc
    ${Geant4_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}
)

link_directories(${ROOT_LIBRARY_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in include_directories#
file(GLOB sources 
    ${PROJECT_SOURCE_DIR}/SiPM/src/*.cc
    ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers 
    ${PROJECT_SOURCE_DIR}/SiPM/inc/*.hh
    ${PROJECT_SOURCE_DIR}/inc/*.hh
)

#----------------------------------------------------------------------------
# Add the executable and link it to the Geant4 libraries
#
add_executable(entrance entrance.cc ${sources} ${headers})
target_link_libraries(entrance ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we 
# build element. This is so that we can run the executable directly because it 
# relies on these scripts being the current working sirectory.
#
file(COPY mac DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this standalone
#
add_custom_target(Entrance DEPENDS entrance)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS entrance DESTINATION bin)

