#pragma once 

#include "CommonTypes.h"

#include <iostream>

using Real3  = CommonTypes::Real3;
using Real   = CommonTypes::Real;

// printing statements 
template <typename T, std::size_t dim>
inline std::ostream& operator<<(std::ostream &out, const std::array<T,dim>& v) {
    for (int i=0;i<dim;i++)
    {
        out << v[i] << " ";
    }
    return out;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
    for (int i=0;i<v.size();i++)
    {
        out << v[i] << " ";
    }

    return out;
}


// std::array
template <typename T, std::size_t dim>
inline std::array<T,dim> operator+(const std::array<T,dim>& v1, const std::array<T,dim>& v2){
    std::array<T,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] + v2[i];
    }

    return ret;
}

template <typename T, std::size_t dim>
inline std::array<T,dim> operator-(const std::array<T,dim>& v1, const std::array<T,dim>& v2){
    std::array<T,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] - v2[i];
    }

    return ret;
}

template <typename T, std::size_t dim>
inline std::array<T,dim> operator*(const std::array<T,dim>& v1, const std::array<T,dim>& v2)
{
    std::array<T,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] * v2[i];
    }
    return ret;
}

template <typename T1, typename T2>
inline std::vector<T1> operator*(const std::vector<T1>& v1, const T2& f)
{
    ASSERT((v1.size() !=0), "Multiply an empty std::vector.");
    std::vector<T1> ret(v1.size());
    for (int i=0;i<v1.size();i++){
        ret[i] = v1[i] * f;
    }
    return ret;
}

template <typename T1, typename T2>
inline std::vector<T1> operator/(const std::vector<T1>& v1, const T2& f)
{
    T2 factor = 1.0/f;
    return v1 * factor;
}

template <typename T1, typename T2>
inline std::vector<T1> operator+(const std::vector<T1>& v1, const T2& f)
{
    ASSERT((v1.size() !=0), "adding an empty std::vector.");
    std::vector<T1> ret(v1.size());
    for (int i=0;i<v1.size();i++){
        ret[i] = v1[i] + f;
    }
    return ret;
}

template <typename T1, typename T2>
inline std::vector<T1> operator+(const std::vector<T1>& v1, std::vector<T2>& v2){
    ASSERT((v1.size() == v2.size()), "The vector operation cannot be done because of size mismatch.");
    std::vector<T1> res(v1.size());
    for (int i=0;i<v1.size();i++){
        res[i] = v1[i] + v2[i];
    }

    return res;
}

template <typename T, std::size_t dim>
inline std::array<T,dim> operator/(const std::array<T,dim>& v1, const std::array<T,dim>& v2)
{
    std::array<T,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] / v2[i];
    }
    return ret;
}

template <typename T1, typename T2, std::size_t dim>
inline std::array<T1,dim> operator+(const std::array<T1,dim>& v1, T2 value)
{
    std::array<T1,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] + value;
    }
    return ret;
}

template <typename T1, typename T2, std::size_t dim>
inline std::array<T1,dim> operator*(const std::array<T1,dim>& v1, T2 value)
{
    std::array<T1,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] * value;
    }

    return ret;
}

template <typename T1, typename T2, std::size_t dim>
inline std::array<T1,dim> operator*(T2 value, const std::array<T1,dim>& v1)
{
    return v1 * value;
}

template <typename T1, typename T2, std::size_t dim>
inline std::array<T1,dim> operator/(const std::array<T1,dim>& v1, T2 value)
{
    std::array<T1,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] / value;
    }

    return ret;
}