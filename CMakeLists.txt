#set minimum cmake version and project name
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TRexGeant4)

#set compile flags
set( CMAKE_CXX_FLAGS "-Wall -g " )


#The aforementioned find_package command is used to locate and configure Geant4 
#(we'll see how to specify the location of the Geant4 installation to be used later when we run CMake),
# the REQUIRED argument being supplied so that CMake will fail with an error if it cannot find Geant4. 
#The option command specifies a boolean variable which defaults to ON , and which can be set when running CMake via 
#a -D command line argument, or toggled in the CMake GUI interfaces. We wrap the calls to find_package in a conditional 
#block on the option value. This allows us to configure the use of Geant4 UI and Visualization drivers by exampleB1 via 
#the ui_all vis_all "component" arguments to find_package . These components and their usage is described in detail in Section 4.1. 
 option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
 if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
 else()
    find_package(Geant4 REQUIRED)
 endif()

 # To automatically configure the header path, compiler flags and compiler definitions needed for linking to Geant4,
 # we use the include command to load a CMake script supplied by Geant4. The CMake variable named Geant4_USE_FILE
 # is set to the path to this module when Geant4 is located by find_package . We recommend the usage of this file as
 # it will automatically set the compiler definitions and flags recommended for use with Geant4. We use the include_directories
 # command to add the header directory to the compiler's header search path. The CMake variable PROJECT_SOURCE_DIR points to
 # the top level directory of the project and is set by the earlier call to the project command. 
 message(STATUS ${Geant4_USE_FILE})
 include(${Geant4_USE_FILE})


# Load some basic macros which are needed later on for root
include(${CMAKE_SOURCE_DIR}/cmakeScripts/FindROOT.cmake)

 include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})
message(STATUS ${ROOT_INCLUDE_DIR})

 #    Use the globbing functionality of the file command to prepare lists of the B1 source and header files.

# Note however that CMake globbing is only used here as a convenience. The expansion of the glob only happens
# when CMake is run, so if you later add or remove files, the generated build scripts will not know a change
# has taken place. Kitware strongly recommend listing sources explicitly as CMake automatically makes the build
# depend on the CMakeLists.txt file. This means that if you explicitly list the sources in CMakeLists.txt, any 
# changes you make will be automatically picked when you rebuild. This is most useful when you are working on a 
# project with sources under version control and multiple contributors.    
 file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
 file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#  The add_executable command defines the build of an application, outputting an executable named by its first argument,
# with the sources following. Note that we add the headers to the list of sources so that they will appear in IDEs like Xcode.
# After adding the executable, we use the target_link_libraries command to link it with the Geant4 libraries. 
# The Geant4_LIBRARIES variable is set by find_package when Geant4 is located, and is a list of all the 
# libraries needed to link against to use Geant4. 
 add_executable(TRexGeant4 TRexGeant4.cc ${sources} ${headers})
 target_link_libraries(TRexGeant4 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})



#  install(TARGETS exampleB1 DESTINATION bin)
