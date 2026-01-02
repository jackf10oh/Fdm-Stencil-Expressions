cmake_minimum_required(VERSION 3.16)
project(LinOps LANGUAGES CXX)

# set cxx standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# adding Eigenlib to compiler
find_package(Eigen3 REQUIRED)

# define the library as LinOps
add_library(LinOps
  INTERFACE
  Mesh.hpp
  Discretization.hpp
  LinOpTraits.hpp
  LinearOpBase.hpp
  LinOpExpr.hpp
  Operators/IOp.hpp
  Operators/RandLinOp.hpp
)

target_include_directories(LinOps 
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(LinOps
  INTERFACE Eigen3::Eigen 
)

# adding  google test
enable_testing()

#find the gtest library 
find_package(GTest CONFIG REQUIRED) # will fail if not found

if(NOT GTest_FOUND)
  message(WARNING "GTest not found. Skipping targets.")

else()
  add_executable(tests 
    tests/tests.cpp
  )

  target_link_libraries(tests  
    PRIVATE LinOps 
    PRIVATE GTest::gtest 
    PRIVATE GTest::gtest_main 
    # PRIVATE GTest::gmock
    # PRIVATE GTest::gmock_main
    # PRIVATE pthread # required by GTest
  )

  include(GoogleTest)
  gtest_discover_tests(tests)
  
endif()