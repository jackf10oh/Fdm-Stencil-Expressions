cmake_minimum_required(VERSION 3.16)
project(Bindings LANGUAGES CXX)

# set cxx standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# allow pybind11 to find python 
set(PYBIND11_FINDPYTHON ON)

# adding Eigenlib to compiler
find_package(Eigen3 REQUIRED)

# Find pybind11
find_package(pybind11 REQUIRED)

add_library(LinOps INTERFACE)
target_link_libraries(LinOps
 INTERFACE Eigen3::Eigen 
)

add_library(FDStencils INTERFACE)
target_link_libraries(FDStencils
 INTERFACE Eigen3::Eigen 
 INTERFACE LinOps
)

# Link the shared library to your pybind11 module
pybind11_add_module(PyFds
  PyFds.cpp  
)

# Link the Python module with your shared library
target_link_libraries(PyFds 
  PRIVATE
  FDStencils
  LinOps  
)

# Apply static linking flags as compile options
target_link_options(PyFds PRIVATE 
  -static-libgcc 
  -static-libstdc++ 
  -lpthread
)


