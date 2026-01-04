// StrideIterator
//
// wrapper for iterator that overrides operator++(), etc... to be in multiples of a stride 
// currently requires Iter_t be a random_access iterator 
//
// JAF 1/3/2026 

#ifndef STRIDEITERATOR_H
#define STRIDEITERATOR_H 

#include<iterator>
#include<memory> // std::addressof()
#include "MemView.hpp" 

template<typename Iter_t>
class StrideIterator
{
  public:
    // typedefs (required to be STL compatible) --------------------------------
    using traits = std::iterator_traits<Iter_t>;

    using difference_type = typename traits::difference_type;
    using value_type = typename traits::value_type;
    using pointer = typename traits::pointer;
    using reference = typename traits::reference;
    using iterator_category = typename traits::iterator_category;
     
    // constructors ------------------------------------------------------------
    // no default constructor
    StrideIterator()=delete; 
    // from base iterator 
    StrideIterator(Iter_t it_init, difference_type stride_init=1) : m_base_iterator(it_init), m_stride(stride_init){} 
    // copy 
    StrideIterator(const StrideIterator& other) : m_base_iterator(other.underlying_iter()), m_stride(other.m_stride){}
    // copy different type of iter_t 
    template<typename U>
    StrideIterator(const StrideIterator<U>& other) : m_base_iterator(other.underlying_iter()),  m_stride(other.stride()){}
    // destructors ------------------------------------------------------------------
    ~StrideIterator()=default; 
    // member functions / operators ==========================================================
    Iter_t& underlying_iter() { return m_base_iterator; }
    const Iter_t& underlying_iter() const { return m_base_iterator; } 
    difference_type stride() const {return m_stride;}
    // write / read 
    reference operator*() const {
        return *m_base_iterator;
    }
    // -> to member func 
    pointer operator->() const {return std::addressof(*m_base_iterator);}

    // pre/post increment 
    StrideIterator& operator++(){m_base_iterator+=m_stride; return *this; }
    StrideIterator operator++(int){StrideIterator temp(*this); m_base_iterator += m_stride; return temp; } 

    // pre/post decrement
    StrideIterator& operator--(){m_base_iterator-=m_stride; return *this; }
    StrideIterator operator--(int){StrideIterator temp(*this); m_base_iterator -= m_stride; return temp; } 

    // differencing 
    difference_type operator-(const StrideIterator& other) const {
      return (m_base_iterator - other.m_base_iterator) / m_stride;
    }

    // equality / inequality 
    bool operator==(const StrideIterator& other) const { return m_base_iterator == other.m_base_iterator; } 
    bool operator!=(const StrideIterator& other) const { return m_base_iterator != other.m_base_iterator; } 

    // member data ===========================================================
    Iter_t m_base_iterator;
    difference_type m_stride; 

}; // end StrideIterator class 


// operators for random access iterators --------------------------------------------------
template<typename Iter_t>
StrideIterator<Iter_t> operator+(const StrideIterator<Iter_t>& s_it, typename StrideIterator<Iter_t>::difference_type n){
  return StrideIterator(s_it.underlying_iter() + s_it.stride() * n, s_it.stride());  
}

template<typename Iter_t>
const StrideIterator<Iter_t>& operator+=(StrideIterator<Iter_t>& s_it, typename StrideIterator<Iter_t>::difference_type n){
  s_it = s_it + n;
  return s_it; 
}

template<typename Iter_t>
StrideIterator<Iter_t> operator-(const StrideIterator<Iter_t>& s_it, typename  StrideIterator<Iter_t>::difference_type n){
  return StrideIterator(s_it.underlying_iter() - s_it.stride() * n, s_it.stride());  
}

// entry points --------------------------------------------------------------------- 
template<typename Cont, typename difference_type = typename Cont::difference_type> 
MemView<StrideIterator<typename Cont::iterator>> make_strided_MemView(Cont& c, difference_type stride=1, difference_type offset=0){
  StrideIterator<typename Cont::iterator> begin(c.begin() + offset, stride); 
  difference_type size = ((c.size()-offset) + stride - 1) / stride; 
  auto end = begin + size;  
  return MemView(begin, end); 
} 

#endif // StrideIterator.hpp 
