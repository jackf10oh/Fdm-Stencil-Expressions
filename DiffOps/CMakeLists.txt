cmake_minimum_required(VERSION 3.16)
project(DiffOps LANGUAGES CXX)

# set cxx standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# options 
option(USE_OPENMP "Enable OpenMP Parallelization" ON) 

# adding OpenMP is enable 
if(USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    message(STATUS "OpenMP found successfully!")
  endif()
else()
  message(STATUS "Cmake args disabling OpenMP...")
endif()


# adding Eigenlib to compiler
find_package(Eigen3 REQUIRED)

# load library LinOps into current CMakeLists.txt
add_library(LinOps INTERFACE)
target_link_libraries(LinOps
 INTERFACE Eigen3::Eigen 
)

# define the library as FdmStencilExressions
add_library(DiffOps INTERFACE)

target_include_directories(DiffOps 
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(DiffOps
  INTERFACE LinOps
  INTERFACE OpenMP::OpenMP_CXX
)

# adding  google test
enable_testing()

#find the gtest library 
find_package(GTest REQUIRED) # will fail if not found
find_package(benchmark) # Google Benchmark will fail if not found
# find_package(benchmark REQUIRED) # Google Benchmark will fail if not found

if(NOT (benchmark_FOUND AND GTest_FOUND))
  message(WARNING "Google Benchmark not found. Skipping benchmark targets")
else()
  message(STATUS "Google Benchmark found. building benchmarks")
  add_executable(benchmarks
    benchmarks/benchmarks.cpp
  ) 
  target_link_libraries(benchmarks
    PRIVATE LinOps 
    PRIVATE benchmark::benchmark
    PRIVATE GTest::gtest 
    PRIVATE GTest::gtest_main 
    PRIVATE pthread # required by GTest
  )
endif()

if(NOT GTest_FOUND)
  message(WARNING "GTest not found. Skipping targets.")

else()
  message(STATUS "GTest found. building tests") 
  add_executable(tests 
    tests/tests.cpp
  )

  target_link_libraries(tests  
    PRIVATE DiffOps 
    PRIVATE GTest::gtest 
    PRIVATE GTest::gtest_main 
    PRIVATE pthread # required by GTest
  )

  include(GoogleTest)
  gtest_discover_tests(tests)
  
endif()