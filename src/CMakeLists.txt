cmake_minimum_required(VERSION 3.8)

# Define the project name.
project(numsim_parallel)

# Specify the name of the executable (${PROJECT_NAME} which is equal to what was set in the project() command).
# Also specify the source files.
set(PROJECT_FILES
    partitioning/partitioning.cpp    

    computation/computation.cpp  
    computation/computation_parallel.cpp  
    
    discretization/0_staggered_grid.cpp
    discretization/1_discretization.cpp
    discretization/2_central_differences.cpp
    discretization/2_donor_cell.cpp

    output_writer/output_writer.cpp
    output_writer/output_writer_text.cpp
    output_writer/output_writer_text_parallel.cpp
    output_writer/output_writer_paraview.cpp
    output_writer/output_writer_paraview_parallel.cpp

    pressure_solver/pressure_solver.cpp
    pressure_solver/gauss_seidel.cpp
    pressure_solver/sor.cpp
    pressure_solver/pressure_solver_parallel.cpp
    pressure_solver/red_black_sor.cpp
    pressure_solver/conjugate_gradient.cpp

    settings/settings.cpp

    storage/array2d.cpp
    storage/field_variable.cpp
)

add_executable(${PROJECT_NAME} main.cpp ${PROJECT_FILES})

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

# Add the project directory to include directories, to be able to include all project header files from anywhere
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

# Search for the external package "VTK"
find_package(VTK)

# Output various cmake variables for demonstration purpose
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")
message("The directory of VTK: VTK_DIR: ${VTK_DIR}")
message("The include directory of VTK: VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})        # add the include directory where the header files are for the compiler
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
endif(VTK_FOUND)

# install numsim executable in build directory
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)

# Add additional compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")