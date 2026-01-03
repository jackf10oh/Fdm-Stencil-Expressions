// MemView.hpp
//
// Non owning view into section of memory according to stored iterators to begin/end 
//
// JAF 1/3/2026 

#include<cstdint>
#include<cmath>
#include<vector>
#include<string> 
#include<iostream> 

// non owning view into array of memory
// template<typename Iter_t> 
class MemView{
  private:
    // type defs 
    using Array_t = std::vector<double>;
    using Iter_t = Array_t::iterator;
    using CIter_t = Array_t::const_iterator;
    using RIter_t = Array_t::reverse_iterator;
    using CRIter_t = Array_t::const_reverse_iterator;

    // member data 
    const Iter_t m_begin_it; 
    const Iter_t m_end_it;

  public:
    // Constructors 
    MemView()=delete; 
    MemView(Iter_t begin, Iter_t end): m_begin_it(begin), m_end_it(end){}; 
    MemView(const MemView& other) : m_begin_it(other.m_begin_it), m_end_it(other.m_end_it){}; 
    
    // destructors
    ~MemView()=default;
    
    // member funcs 
    std::size_t size() const {return std::distance(m_begin_it,m_end_it);};
    double at(std::size_t i){if(i>=size()) throw std::out_of_range("index >= size");return *(m_begin_it+i);}
    // Iterators 
    const Iter_t& begin() const {return m_begin_it;};  
    const Iter_t& end() const {return m_end_it;};  
    CIter_t cbegin() const {return m_begin_it;};  
    CIter_t cend() const {return m_end_it;};  
    RIter_t rbegin() const { return std::make_reverse_iterator(m_end_it);}; 
    RIter_t rend() const { return std::make_reverse_iterator(m_begin_it);};
    CRIter_t crbegin() const { return std::make_reverse_iterator(m_end_it);}; 
    CRIter_t crend() const { return std::make_reverse_iterator(m_end_it);}; 

    // operators
    double operator[](std::size_t i){return *(m_begin_it+i);}
};