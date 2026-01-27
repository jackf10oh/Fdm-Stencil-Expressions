cmake_minimum_required(VERSION 3.16)
project(FdmStencilExprs LANGUAGES CXX)
add_library(FdmStencilExprs INTERFACE)

# set cxx standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# options 
option(USE_OPENMP "Enable OpenMP Parallelization" OFF) 

# adding OpenMP is enable 
if(USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    message(STATUS "OpenMP found successfully!")
    target_link_libraries(FdmStencilExprs INTERFACE OpenMP::OpenMP_CXX)
  endif()
else()
  message(STATUS "Cmake args disabling OpenMP...")
endif()

# adding Eigenlib to compiler
find_package(Eigen3 REQUIRED)

# header only libraries ===============================
# Utilities ------------- 
add_library(Utilities INTERFACE) 
target_include_directories(Utilities 
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Utilities/include 
)
target_link_libraries(Utilities 
  INTERFACE Eigen3::Eigen
)


# LinOps ------------- 
add_library(LinOps INTERFACE)
target_include_directories(LinOps
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/LinOps/include
)
target_link_libraries(LinOps  
  INTERFACE Eigen3::Eigen
)

# LinOpsXD ------------- 
add_library(LinOpsXD INTERFACE)
target_include_directories(LinOpsXD
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/LinOpsXD/include
)
target_link_libraries(LinOpsXD
  INTERFACE LinOps
)

# FDStencils ------------- 
add_library(FDStencils INTERFACE)
target_include_directories(FDStencils
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/FDStencils/include
)
target_link_libraries(FDStencils 
  INTERFACE LinOps
  INTERFACE Utilities
)

# FDStencilsXD ------------- 
add_library(FDStencilsXD INTERFACE)
target_include_directories(FDStencilsXD 
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/FDStencilsXD/include
)
target_link_libraries(FDStencilsXD
  INTERFACE LinOps
  INTERFACE LinOpsXD
  INTERFACE FDStencils
  INTERFACE Utilities
)

# TExprs --------------------------
add_library(TExprs INTERFACE)
target_include_directories(TExprs
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/TExprs/include
)
target_link_libraries(TExprs
  INTERFACE LinOps 
  INTERFACE LinOpsXD
  INTERFACE FDStencils 
  INTERFACE FDStencilsXD
  INTERFACE Utilities
)

# Tie all sub directories into main library ============================================== 
target_link_libraries(FdmStencilExprs
  INTERFACE LinOps 
  INTERFACE LinOpsXD
  INTERFACE FDStencils 
  INTERFACE FDStencilsXD
  INTERFACE TExprs 
  INTERFACE Utilities
)

# adding google test
enable_testing()
find_package(GTest REQUIRED) #find the gtest library  # will fail if not found

# adding google benchmark  
# find_package(benchmark) # Google Benchmark will fail if not found
# find_package(benchmark REQUIRED) # Google Benchmark will fail if not found

# if(NOT (benchmark_FOUND AND GTest_FOUND))
#   message(WARNING "Google Benchmark not found. Skipping benchmark targets")
# else()
#   message(STATUS "Google Benchmark found.") #  building benchmarks
#   # add_executable(benchmarks
#   #   benchmarks/benchmarks.cpp
#   # ) 
#   # target_link_libraries(benchmarks
#   #   PRIVATE LinOps 
#   #   PRIVATE benchmark::benchmark
#   #   PRIVATE GTest::gtest 
#   #   PRIVATE GTest::gtest_main 
#   #   PRIVATE pthread # required by GTest
#   # )

#   # if(USE_OPENMP)
#   #   target_link_libraries(benchmarks
#   #     PRIVATE OpenMP::OpenMP_CXX
#   #   )
#   # endif()

# endif()

if(NOT GTest_FOUND)
  message(WARNING "GTest not found. Skipping targets.")

else()
  message(STATUS "GTest found.") # building tests 

  add_subdirectory(tests)
  # add_executable(tests 
  #   tests/tests.cpp
  # )

  # target_link_libraries(tests  
  #   PRIVATE LinOps 
  #   PRIVATE LinOpsXD 
  #   PRIVATE FDStencils 
  #   PRIVATE FDStencilsXD 
    # PRIVATE GTest::gtest 
    # PRIVATE GTest::gtest_main 
    # PRIVATE pthread # required by GTest
  # )

  # include(GoogleTest)
  # gtest_discover_tests(tests)
  
endif()

# on windows with mingw 
# cmake .. `
#   -G "MinGW Makefiles" `
#   -DCMAKE_C_COMPILER=C:/msys64/mingw64/bin/gcc.exe `
#   -DCMAKE_CXX_COMPILER=C:/msys64/mingw64/bin/g++.exe `
#   -DEigen3_DIR=C:/eigen-install/share/eigen3/cmake `
#   -DGTest_DIR=C:\gtest-mingw\lib\cmake\GTest