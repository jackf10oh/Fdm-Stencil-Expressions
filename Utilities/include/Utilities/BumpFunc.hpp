// BumpFunc.hpp
//
// Bump function on [L,R] with center c and height h
//
// JAF 1/13/2026

#ifndef BUMPFUNC_H
#define BUMPFUNC_H 

struct BumpFunc
{
    double L = 0.0;
    double R = 1.0; 
    double c = 0.5;
    double h = 1.0; 
    double focus = 3.0; 
    
    // Member Funcs 
    double operator()(double x) const 
    {
      // if( c<L || x<L ||c>R || x>R) throw std::runtime_error("Error: Bad args to BumpFunc operator().");
      if( c<L || x<L ||c>R || x>R) return 0.0;
      // get a,b in Beta(a,b) distribution  
      double a = focus * ((c-L)/(R-L)); 
      double b = focus * ((R-c)/(R-L)); 
      // scaling coefficient to make pdf(mode) == h
      double scale = h / (std::pow(R*(a/(a+b)), a)*std::pow(R*(b/(a+b)), b)); 
      // value of h*Beta(a,b) at (x-L) / (R-L) 
      return scale * std::pow(x-L,a) * std::pow(R-x+L,b); 
    }
};

#endif // BumpFunc.hpp