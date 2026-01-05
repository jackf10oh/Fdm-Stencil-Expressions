// PrintVec.hpp 
//
//
//
// JAF 1/4/2026 

#ifndef PRINTVEC_H
#define PRINTVEC_H

#include<iostream>
#include<Eigen/Core> 

// from container with .begin() .end() 
template<typename Cont>
void print_vec(const Cont& v, std::string comment="", bool new_line=true){
  if(!comment.empty()) std::cout << comment << ": ";
  auto it = v.begin(); 
  auto end = std::prev(v.end()); 
  std::cout << "["; 
  while(it!= end) std::cout << *(it++) << ", ";
  std::cout << *it << "]"; 
  if(new_line) std::cout << std::endl;   
};

// from iterators. 
template<typename Iter>
void print_vec(const Iter& start, const Iter& stop, std::string comment="", bool new_line=true){
  if(!comment.empty()) std::cout << comment << ": ";
  auto it = start; 
  auto end = std::prev(stop); 
  std::cout << "["; 
  while(it!= end) std::cout << *(it++) << ", ";
  std::cout << *it << "]";   
  if(new_line) std::cout << std::endl;   
};

// for an eigen MatrixXd
void print_mat(const Eigen::Ref<const Eigen::MatrixXd>& A, std::string comment="")
{
  if(!comment.empty()) std::cout << "--------------" << comment << "--------------" << std::endl;
  for(auto row : A.rowwise()){
    print_vec(row); 
    std::cout << "," << std::endl; 
  }
  std::cout << std::endl; 
}

#endif // PrintVec.hpp