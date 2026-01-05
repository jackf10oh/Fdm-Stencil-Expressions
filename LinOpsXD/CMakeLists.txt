cmake_minimum_required(VERSION 3.16)
project(LinOpsXD LANGUAGES CXX)

# set cxx standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# adding Eigenlib to compiler
find_package(Eigen3 REQUIRED)

# load library LinOps into current CMakeLists.txt
add_library(LinOps INTERFACE)
target_link_libraries(LinOps
 INTERFACE Eigen3::Eigen 
)

# define the library as LinOpsXD
add_library(LinOpsXD
  INTERFACE
  MeshXD.hpp
  DiscretizationXD.hpp
  LinOpXDTraits.hpp
  LinearOpXDBase.hpp
  # LinOpExpr.hpp
  OperatorsXD/IOpXD.hpp
  OperatorsXD/DirectionalRandOp.hpp
)

target_include_directories(LinOpsXD 
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(LinOpsXD
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
    PRIVATE LinOpsXD 
    PRIVATE GTest::gtest 
    PRIVATE GTest::gtest_main 
    # PRIVATE GTest::gmock
    # PRIVATE GTest::gmock_main
    # PRIVATE pthread # required by GTest
  )

  include(GoogleTest)
  gtest_discover_tests(tests)
  
endif() 

# on windows with mingw 
# cmake .. `
#   -G "MinGW Makefiles" `
#   -DCMAKE_C_COMPILER=C:/msys64/mingw64/bin/gcc.exe `
#   -DCMAKE_CXX_COMPILER=C:/msys64/mingw64/bin/g++.exe `
#   -DEigen3_DIR=C:/eigen-install/share/eigen3/cmake `
#   -DGTest_DIR=C:\gtest-mingw\lib\cmake\GTest
      