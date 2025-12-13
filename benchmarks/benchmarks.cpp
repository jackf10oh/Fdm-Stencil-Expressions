// benchmarks.cpp
//
//
//
// JAF 12/12/2025

// compilation command 
// g++ -O3 -DNDEBUG     ./benchmarks/benchmarks.cpp  -lbenchmark -lpthread     -o benchmarks_main

#include<benchmark/benchmark.h>

#include "../incl/All.hpp"

static void BENCHMARK_NthDerivOp_set_mesh(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto my_mesh = make_mesh(0.0,double(r), r+1); 
    NthDerivOp D(1); 
    state.ResumeTiming(); 
    D.set_mesh(my_mesh); 
  }
};

static void BENCHMARK_ConvectionDiffusionExpression_set_mesh(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming(); 
    auto r = 10.0; 
    int n_gridpoints = state.range(0);

    // uniform mesh from 0.0 to r with n_gridpoints
    MeshPtr_t my_mesh = make_mesh(0.0,r,n_gridpoints);

    // create an fdm_scheme for convection diffusion equation 
    using D = NthDerivOp; 
    double dt = 0.01;
    auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2) + 0.5*D(1));

    state.ResumeTiming();
    fdm_scheme.set_mesh(my_mesh);
  };
}

static void BENCHMARK_ConvectionDiffusionExpression_explicit_step(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming(); 
    auto r = 10.0; 
    int n_gridpoints = state.range(0);

    // uniform mesh from 0.0 to r with n_gridpoints
    MeshPtr_t my_mesh = make_mesh(0.0,r,n_gridpoints);

    Discretization1D my_vals;
    
    // Bump centered at r/2. Zero at 0.0 and r 
    auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; 
    
    // fill my_vals with func(x)
    my_vals.set_init(my_mesh, func); 
    
    // create an fdm_scheme for convection diffusion equation 
    using D = NthDerivOp; 
    double dt = 0.01;
    auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2) + 0.5*D(1));

    // set the boundary conditions to Dirichlet 0
    auto left = std::make_shared<DirichletBC>(0.0);
    auto right = left; 
    fdm_scheme.lbc_ptr = left; 
    fdm_scheme.rbc_ptr = right; 

    // fill out the stencil matrix
    fdm_scheme.set_mesh(my_mesh);

    state.ResumeTiming();
    benchmark::DoNotOptimize(fdm_scheme.apply(my_vals));

  };
}

BENCHMARK(BENCHMARK_NthDerivOp_set_mesh)->Arg(1e4)->Arg(1e5)->Arg(1e6);
BENCHMARK(BENCHMARK_ConvectionDiffusionExpression_set_mesh)->Arg(1e4)->Arg(1e5)->Arg(1e6);
BENCHMARK(BENCHMARK_ConvectionDiffusionExpression_explicit_step)->Arg(1e4)->Arg(1e5)->Arg(1e6);

BENCHMARK_MAIN();