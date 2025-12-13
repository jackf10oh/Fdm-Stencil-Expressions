// Fornberg.hpp
//
// header file for stateful fornberg calc that only allocates on construction 
//
// JAF 10/29/2025

#ifndef FORNBERGCALC_H
#define FORNBERGCALC_H

#include<cmath>
#include<vector>
#include<string> 
#include<iostream> 

// traits 
using Array_t = std::vector<double>;
using Iter_t = Array_t::iterator;
using CIter_t = Array_t::const_iterator;
using RIter_t = Array_t::reverse_iterator;
using CRIter_t = Array_t::const_reverse_iterator;

// non owning view into array of memory
class WeightsResult_t{
  private:
    // member data 
    const Iter_t m_begin_it; 
    const Iter_t m_end_it;
  public:
    // Constructors 
    WeightsResult_t()=delete; 
    WeightsResult_t(Iter_t begin, Iter_t end): m_begin_it(begin), m_end_it(end){}; 
    WeightsResult_t(const WeightsResult_t& other)=delete; 
    
    // destructors
    ~WeightsResult_t()=default;
    
    // member funcs 
    std::size_t size() const {return std::distance(m_begin_it,m_end_it);};
    double at(std::size_t i){if(i>=size()) throw std::out_of_range("index >= size");return *(m_begin_it+i);}
    // Iterators 
    const Iter_t& begin() const {return m_begin_it;};  
    const Iter_t& end() const {return m_end_it;};  
    CIter_t cbegin() const {return m_begin_it;};  
    CIter_t cend() const {return m_end_it;};  
    RIter_t rbegin() const { return std::make_reverse_iterator(std::prev(m_end_it));}; 
    RIter_t rend() const { return std::make_reverse_iterator(std::prev(m_begin_it));};
    CRIter_t crbegin() const { return std::make_reverse_iterator(std::prev(m_end_it));}; 
    CRIter_t crend() const { return std::make_reverse_iterator(std::prev(m_end_it));}; 

    // operators
    double operator[](std::size_t i){return *(m_begin_it+i);}
};

class FornCalc
{
  public:
    // member data
    // maximum order of derivative stencil 
    std::size_t m_order; 
    std::size_t m_n_nodes; 
    // single allocation of memory rows*cols big 
    std::vector<double> m_arr; 
  public:
    FornCalc()=delete;
    FornCalc(std::size_t max_nodes, std::size_t max_order_init=1)
      : m_n_nodes(max_nodes), m_order(max_order_init), m_arr(m_n_nodes*(m_order+1))
    {
      // if(m_order+1>max_nodes) throw std::invalid_argument("# Nodes must be > Order deriv in Fornberg algo"); 
      // for(double& val : m_arr) val=0.0; 
    };
    FornCalc(const FornCalc& other)=delete; 
    ~FornCalc()=default; 
    template<typename Input_Iter>
    WeightsResult_t GetWeights(double x_bar, Input_Iter start, Input_Iter end, std::size_t order=1)
    {
      // Matrix of Order+1 rows, N cols
      // row m from Weights is the coeffs for derivative of order m (m= 0, ... , order)
      // ---> we need atleast Order+1 * N entries in m_arr 
      if(m_n_nodes*m_order > m_arr.size()) throw std::runtime_error("nodes*(order+1) exceeds size of stored array!"); 
      // m_arr.resize(m_n_nodes*(order+1));

      // utility lambdas 
      auto entryRef = [this](std::size_t i, std::size_t j)->double&{ return this->m_arr[i*this->m_n_nodes+j];}; 

      // number of nodes 
      m_n_nodes = std::distance(start,end); 
      m_order = order;
      auto nodeRef = [start](std::size_t i)->double& {return *(start+i);};
      
      // zero all stored entries
      for(auto& entry : m_arr) entry=0.0; 

      // Using 1 nodes is just a flat line
      entryRef(0,0) = 1; 

      // c1 holds old c2 for next loop
      double c1=1.0; 
      // c2 will hold an accumulation of (node[n]-node[0]) * ... * (node[n]-node[n-1])
      double c2; 
      // c3 holds the difference (nodes[new]-nodes[old])
      double c3; 

      // for number of nodes n=2, ..., N (first node was zero index)
      for(std::size_t n=1; n<m_n_nodes; n++)
      {
        // std::cout << "N=" << n << std::endl;
        // c1 *= (x_bar-nodes[n-1]);
        // reset c5. it depends on node[n]
        c2=1.0; 
        // loop over all previous nodes 0, ..., n-1
        for(std::size_t old_n=0; old_n<n; old_n++)
        { 
          // accumulation is updated 
          c3 = (nodeRef(n)-nodeRef(old_n)); 
          c2 *= c3;

          // if old_n == n-1 we must use the very last old column 
          // to update the newest nodes column
          // this is because of the formulas. we cannot overwrite them yet!
          if(old_n==n-1)
          {
            // std::cout << "SPECIAL" << std::endl;
            // the newest weight at weight n is calculated
            // its derivates are calculated now too 
            for(std::size_t m=std::min(n,m_order); m>=1; m--)
            {
              // (3.9) ------------------------------------- 
              entryRef(m,n) = (c1/c2) * ( m*entryRef(m-1,n-1) - (nodeRef(n-1)-x_bar)*entryRef(m,n-1) ); 
              // std::cout <<"old_n:"<<old_n+1<<" m:"<<m<<" weight[m][n]:" << Weights[m][n]<<std::endl;
            }
            // (3.7) ------------------------------
            entryRef(0,n) = (c1/c2) * (x_bar - nodeRef(n-1)) * entryRef(0,n-1);
            // std::cout <<"old_n:"<<old_n+1<<" m:0 weight[0][n]:" << Weights[0][n]<<std::endl;
            // std::cout << "END SPECIAL" << std::endl;      
          }

        // all previous weights F(x_bar) are updated according to the rule 
        // F[n,v](x_bar) = ((x_bar-node[n]) / (node[v]-node[n])) * F[n-1,v](x)
        // the derivates are updated as well 
        for(std::size_t m=std::min(n,order); m>=1; m--)
        {
          // (3.8) ----------------------------------------
          entryRef(m,old_n) = (((nodeRef(n)-x_bar)*entryRef(m,old_n)) - (m*entryRef(m-1,old_n))) / c3;
          // std::cout <<"old_n:"<<old_n<<" m:"<<m<<" weight[m][old_n]:" << Weights[m][old_n]<<std::endl;
        }
        // (3.6) -------------------------------------------------
        entryRef(0,old_n) = ((x_bar-nodeRef(n))/(nodeRef(old_n)-nodeRef(n)))*entryRef(0,old_n);
        // std::cout <<"old_n:"<<old_n<<" m:0 weight[0][old_n]:" << Weights[0][old_n] <<std::endl;

        // next old node
        }

        // accumulator is stored for next loop 
        c1 = c2; 
      
        // next new node
        // std::cout << std::endl;
      }

      // Weights now contains LaGrange Interpolant Polynomials for 
      // nodes a0, a1, ..., an evaluated at x_bar
      // return non owning view of last row of Weights;
      return WeightsResult_t(m_arr.begin()+(m_n_nodes*order), m_arr.begin()+(m_n_nodes*(order+1))); 
    };

};

#endif // FornbergCalc.hpp