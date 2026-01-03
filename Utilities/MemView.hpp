// MemView.hpp
//
// Non owning view into section of memory according to stored iterators to begin/end 
//
// JAF 1/3/2026 

#ifndef MEMVIEW_H
#define MEMVIEW_H 

#include<cstdint>
#include<cmath>
#include<vector>
#include<string> 
#include<iostream> 

// non owning view into array of memory
template<typename Iter_t> 
class MemView{
  public:
    // typedefs 
    using iterator = Iter_t; 
    using const_iterator = Iter_t; // allow Iter_t itself to define const ness  
    using reverse_iterator = std::reverse_iterator<Iter_t>; 
    using const_reverse_iterator =  std::reverse_iterator<const_iterator>; // delegate to const_iterator 

    using value_type = typename std::iterator_traits<Iter_t>::value_type;
    using difference_type = typename std::iterator_traits<Iter_t>::difference_type;
    using reference = typename std::iterator_traits<Iter_t>::reference;
    using pointer = typename std::iterator_traits<Iter_t>::pointer;

  private:
    // member data 
    const iterator m_begin_it; 
    const iterator m_end_it;

  public:
    // Constructors 
    MemView()=delete; 
    MemView(iterator begin, iterator end): m_begin_it(begin), m_end_it(end){}; 
    MemView(const MemView& other) : m_begin_it(other.m_begin_it), m_end_it(other.m_end_it){}; 
    
    // destructors
    ~MemView()=default;
    
    // member funcs 
    difference_type size() const {return std::distance(m_begin_it,m_end_it);};
    value_type at(difference_type i){if(i>=size()) throw std::out_of_range("index >= size");return *(m_begin_it+i);}
    // Iterators 
    const iterator& begin() const {return m_begin_it;};  
    const iterator& end() const {return m_end_it;};  
    const_iterator cbegin() const {return m_begin_it;};  
    const_iterator cend() const {return m_end_it;};  
    reverse_iterator rbegin() const { return std::make_reverse_iterator(m_end_it);}; 
    reverse_iterator rend() const { return std::make_reverse_iterator(m_begin_it);};
    const_reverse_iterator crbegin() const { return std::make_reverse_iterator(m_end_it);}; 
    const_reverse_iterator crend() const { return std::make_reverse_iterator(m_end_it);}; 

    // operators
    value_type operator[](std::size_t i){return *(m_begin_it+i);}
};

// Utility functions ================================================================================
template<typename Cont>
void print_vec(const Cont& v, std::string comment="")
{
  if(!comment.empty()) std::cout << comment << ": "; 
  auto it = v.begin(); 
  std::cout<< "["; 
  while(it!=std::prev(v.end()))
  {
    std::cout << *it << ", ";
    it++;
  }
  std::cout << *it << "]" << std::endl; 
}

template<typename Iter>
void print_vec(const Iter& start, const Iter& end, std::string comment="")
{
  if(!comment.empty()) std::cout << comment << ": "; 
  auto it = start; 
  std::cout<< "["; 
  while(it!=end-1)
  {
    std::cout << *it << ", ";
    it++;
  }
  std::cout << *it << "]" << std::endl; 
}; 

#endif // MemView.hpp 
