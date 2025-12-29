// Traits.hpp
//
//
//
// JAF 12/11/2025

#ifndef FDMSTENCIL_TRAITS_H
#define   FDMSTENCIL_TRAITS_H

#include<type_traits>
#include "FdmPlugin.hpp"
#include "CoeffOps/CoeffOpBase.hpp"

// Given a type, see if it is derived from CoeffOpBase's crtp scheme
template<typename T, typename = void> 
struct is_coeffop_crtp : std::false_type{}; 

template<typename T>
struct is_coeffop_crtp<T, std::void_t<typename std::remove_cv_t<std::remove_reference_t<T>>::is_coeff_flag>>: std::true_type{}; 

// given a type, see if it has a .sparseView() method
template<typename T, typename = void>
struct has_sparseview_method : std::false_type{}; 

template<typename T>
struct has_sparseview_method<T, std::void_t<decltype(std::declval<T>().sparseView())>> : std::true_type{}; 

#endif