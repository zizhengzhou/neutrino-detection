#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.27)
project(MKID_reactor_sim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# See the documentation for a guide on how to enable/disable specific components
#
find_package(Geant4 REQUIRED ui_all vis_all)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, use our local headers, and link it to the Geant4 libraries
#
add_executable(mkid_sim main.cc ${sources} ${headers})
target_include_directories(mkid_sim PRIVATE include)
target_link_libraries(mkid_sim PRIVATE ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build mkid_sim. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
# set(EXAMPLEB1_SCRIPTS
#   exampleB1.in
#   exampleB1.out
#   init_vis.mac
#   run1.mac
#   run2.mac
#   vis.mac
#   )

# foreach(_script ${EXAMPLEB1_SCRIPTS})
#   configure_file(
#     ${PROJECT_SOURCE_DIR}/${_script}
#     ${PROJECT_BINARY_DIR}/${_script}
#     COPYONLY
#     )
# endforeach()

#----------------------------------------------------------------------------
# Copy all macro files to the build directory
file(GLOB macros ${PROJECT_SOURCE_DIR}/*.mac)
file(COPY ${macros} DESTINATION ${CMAKE_BINARY_DIR})
