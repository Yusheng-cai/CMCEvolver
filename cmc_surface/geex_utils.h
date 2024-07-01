///////////////////////////////////////////////////////////////////////////////
//                                                                            //
// CMC Surface Evolver                                                        //
// http://i.cs.hku.hk/~hpan/code/cmc_source.zip        						  //
//                                                                            //
// CMC Surface Evolver is a demo implementation of the surface modeling method//
// proposed in the paper "Robust Modeling of Constant Mean Curvature Surfaces"//
//                                                                            //
// Author: Hao Pan (hpan@cs.hku.hk)                                           //
// Last Update: 2014-6-13                                                     //
//																			  //
// The software is freely available for non-commercial purposes.              //
//																			  //
////////////////////////////////////////////////////////////////////////////////



//! Code in this file is mostly adapted from the GEEX free software (Graphics Environment for EXperimentations.
//! Copyright (C) 2006 INRIA - Project ALICE). It contains utility functions used in CMC surface computation.



#ifndef _GEEX_UTIL_ADAPTION_
#define _GEEX_UTIL_ADAPTION_

#include <iostream>
#include <cfloat>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>

#include <assert.h>

#define gx_debug_assert assert

namespace Geex {

	//-------------------------- Check number validity ----------------------------------------------------
	template <class T>
	bool is_nan(T x) {
	#ifdef WIN32
		return (_std::isnan(x) != 0) || (_finite(x) == 0) ;
	#else
		return std::isnan(x) || !finite(x);
	#endif
	}    
    
    //--------------------------- VEC3 ------------------------------------------------------------------------------------------

    template <class T>
    struct vec3g {
        typedef vec3g<T> thisclass ;
        
        vec3g() : x(0), y(0), z(0) { }
        vec3g(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {  }
            
        inline T length2() const { return x*x+y*y+z*z ; }
        inline T length() const { return sqrt(x*x+y*y+z*z) ; }
            
        // operators
        inline thisclass& operator+=(const thisclass& v) { x += v.x ; y += v.y ; z += v.z ; return *this ; }
        inline thisclass& operator-=(const thisclass& v) { x -= v.x ; y -= v.y ; z -= v.z ; return *this ; }
        inline thisclass& operator*=(const thisclass& v) { x *= v.x ; y *= v.y ; z *= v.z ; return *this ; }
        inline thisclass& operator/=(const thisclass& v) { x /= v.x ; y /= v.y ; z /= v.z ; return *this ; }
        template <class T2> inline thisclass& operator*=(T2 s) { x *= T(s) ; y *= T(s) ; z *= T(s) ; return *this ; }
        template <class T2> inline thisclass& operator/=(T2 s) { x /= T(s) ; y /= T(s) ; z /= T(s) ; return *this ; }

		//<<<<<

		inline bool operator==(const thisclass& v) { return (x - v.x)*(x - v.x) + (y - v.y)*(y - v.y) + (z - v.z)*(z - v.z) < 1e-50 ? true:false;  }

		//>>>>>
    
        inline thisclass operator+ (const thisclass& v) const {return thisclass(x+v.x, y+v.y, z+v.z); }
        inline thisclass operator- (const thisclass& v) const {return thisclass(x-v.x, y-v.y, z-v.z); }
        template <class T2> inline thisclass operator* (T2 s) const {return thisclass(x*T(s), y*T(s), z*T(s)); }
        template <class T2> inline thisclass operator/ (T2 s) const {return thisclass(x/T(s), y/T(s), z/T(s)); }

        inline thisclass operator- () const {return thisclass(-x, -y, -z);}

        inline T& operator[](int idx) {
            switch(idx) {
            case 0: return x ; break ;
            case 1: return y ; break ;
            case 2: return z ; break ;
            default:  ;
            }
            return x ;
        }

        inline const T& operator[](int idx) const {
            switch(idx) {
            case 0: return x ; break ;
            case 1: return y ; break ;
            case 2: return z ; break ;
            default:  ;
            }
            return x ;
        }
        
        T x ;
        T y ;
        T z ;
    } ;
        
    template <class T> inline T dot(const vec3g<T>& v1, const vec3g<T>& v2) {  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;  }

    template <class T> inline  vec3g<T> cross(const vec3g<T>& v1, const vec3g<T>& v2) {
        return vec3g<T>(
            v1.y*v2.z - v1.z*v2.y,
            v1.z*v2.x - v1.x*v2.z,
            v1.x*v2.y - v1.y*v2.x
        ) ;
    }
        
    template <class T> inline vec3g<T> operator-(const vec3g<T>& v1) { return vec3g<T>(-v1.x, -v1.y, -v1.z) ; }
    template <class T2, class T> inline vec3g<T> operator*(T2 s, const vec3g<T>& v) { 
        return vec3g<T>(T(s)*v.x, T(s)*v.y, T(s)*v.z) ;   
    }

    template <class T> inline vec3g<T> operator+(const vec3g<T>& v1, const vec3g<T>& v2) { 
        return vec3g<T>(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z) ; 
    }

    template <class T> inline vec3g<T> operator-(const vec3g<T>& v1, const vec3g<T>& v2) { 
        return vec3g<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z) ; 
    }

    // Compatibility with GLSL
    template <class T> inline T length(const vec3g<T>& v) { return v.length() ; }
    template <class T> inline T length2(const vec3g<T>& v) { return v.length2() ; }
    template <class T> inline T distance2(const vec3g<T>& v1, const vec3g<T>& v2) { return length2(v2 - v1) ; }
    template <class T> inline T distance(const vec3g<T>& v1, const vec3g<T>& v2) { return length(v2 - v1) ; }
    template <class T> inline vec3g<T> normalize(const vec3g<T>& v) { return (T(1) / length(v)) * v ;   }
    template <class T> inline vec3g<T> mix(const vec3g<T>& v1, const vec3g<T>& v2, T s) { return (T(1) - s) * v1 + s * v2 ; }

	 template <class T> inline std::ostream& operator<<(std::ostream& out, const Geex::vec3g<T>& v) {
        return out << v.x << "  " << v.y << "  " << v.z  ;
    }

    template <class T> inline std::istream& operator>>(std::istream& in, Geex::vec3g<T>& v) {
        return in >> v.x >> v.y >> v.z ;
    }

	typedef vec3g<double> vec3 ;

	//--------------------------- triangle area ---------------------------------------------------------------------------
	template<class T>
	T tri_area(const vec3g<T>& v0, const vec3g<T>& v1, const vec3g<T>& v2){
		return cross(v1-v0,v2-v0).length()*0.5;
	}

	//--------------------------- CVT energy -------------------------------------------------------------------------------

	// Lloyd energy of a triangle which belong to Voronoi region of seed
	template <class T> inline T Lloyd_energy(
		const vec3g<T>& seed, const vec3g<T>& p1, const vec3g<T>& p2, const vec3g<T>& p3,
		vec3g<T>& gradient, T& V
		) {
			vec3g<T> v0 = p1 - seed ;
			vec3g<T> v1 = p2 - seed ;
			vec3g<T> v2 = p3 - seed ;

			V = tri_area(p1, p2, p3);

			gradient = 2.0 * V *(seed - (1.0/3.0) * (p1+p2+p3));		

			vec3g<T> v_12 = v1+v2;
			vec3g<T> v_012 = v0+v_12;
			return V * (dot(v0,v_012)+dot(v1,v_12)+dot(v2,v2)) / 6.;
			//i.e return V * (dot(v0,v0)+dot(v0,v1)+dot(v0,v2)+dot(v1,v1)+dot(v1,v2)+dot(v2,v2)) / 6;
	}

	//------------------------------ Small set ----------------------------------------------

	  /**
     * Similar to std::set, but with fixed maximum size
     * (and no dynamic memory allocation). Used to store
     * vertices equations (represented as plane indices 
     * triplets).
     */
    template <class T, int DDIM> class small_set {
    public:
        typedef small_set<T,DDIM> thisclass ;
        typedef T* iterator ;
        typedef const T* const_iterator ;
        typedef T& reference ;
        typedef T value_type ;
        
        small_set() : end_(data_) { }

        small_set(const thisclass& rhs) {
            copy(rhs) ;
        }
        
        thisclass& operator=(const thisclass& rhs) {
            copy(rhs) ; return *this ;
        }

        template <int DDIM2> small_set(const small_set<T,DDIM2>& rhs) {
            copy(rhs) ;
        }
        
        template <int DDIM2> thisclass& operator=(const small_set<T,DDIM2>& rhs) {
            copy(rhs) ; return *this ;
        }

        unsigned int size() const { return (unsigned int)(end_ - data_) ; }
        unsigned int capacity() const { return (unsigned int)DDIM ; }

        iterator begin() { return data_ ; }
        iterator end() { return end_ ; }
        iterator end_of_storage() { return data_ + DDIM ; }

        iterator* end_ptr() { return &end_ ; }

        const_iterator begin() const { return data_ ; }
        const_iterator end() const { return end_ ; }
        const_iterator end_of_storage() const { return data_ + DDIM ; }

        iterator insert(const T& x) { 
            return insert(x, find_i(x)) ; 
        }

        iterator insert(const T& x, iterator where) {
            if(where == end()) { *where = x ; grow() ; return where ; }
            if(*where == x) { return where ; }
            grow() ;
            if(where == end() - 1) {
				std::cout << ">>>>>>>>Error out of set bound<<<<<<<<" << std::endl; //hpan
                *where = x ; return where ; 
            }
            for(iterator i = end()-1; i != where; i--) {
                gx_debug_assert(i != begin()) ;
                *i = *(i-1) ;
            }
            *where = x ;
#ifdef GX_DEBUG
            for(iterator i = begin(); i != end() - 1; i++) {
                gx_debug_assert(*i < *(i+1)) ;
            }
#endif
            return where ;
        }

		//<<<<<
		void insert(const thisclass& rhs) {
			for (const_iterator i = rhs.begin(); i != rhs.end(); ++i)
			{
				this->insert(*i, find_i(*i));
			}
		}

		inline bool operator<(const thisclass& rhs) const{
			int com_size = (size() < rhs.size()) ? size() : rhs.size();
			for (int i = 0; i < com_size; i++)
			{
				if ((*this)[i] < rhs[i])
				{
					return true;
				}
				else
					if ((*this)[i] > rhs[i])
					{
						return false;
					}
			}
			if (size() < rhs.size())
				return true;
			return false;
		}

		inline bool operator==(const thisclass& rhs) const{
			if (size() != rhs.size())
			{
				return false;
			}
			
			for (int i = 0; i < size(); i++)
			{
				if ((*this)[i] != rhs[i])
				{
					return false;
				}
			}

			return true;
		}

		//>>>>>

        void clear() { end_ = data_ ; }

        iterator find(const T& x) {
            iterator result = find_i(x) ;
            if(*result != x) { result = end() ; }
            return result ;
        }

        const_iterator find(const T& x) const {
            const_iterator result = find_i(x) ;
            if(*result != x) { result = end() ; }
            return result ;
        }

        void push_back(const T& x) {
#ifdef GX_DEBUG
            for(iterator i = begin(); i != end(); i++) {
                gx_debug_assert(*i < x) ;
            }
#endif
            *end_ = x ;
            grow() ;
        }

        void print(std::ostream& out) const {
            out << "[ " ;
            for(const_iterator it=begin(); it!=end(); it++) {
                out << *it << " " ;
            }
            out << "]" ;
        }

        T& operator[](int i) {
            gx_debug_assert(i >= 0) ;
            gx_debug_assert(begin() + i < end()) ;
            return begin()[i] ;
        }

        const T& operator[](int i) const {
            gx_debug_assert(i >= 0) ;
            gx_debug_assert(begin() + i < end()) ;
            return begin()[i] ;
        }

    protected:

        void grow() {
            gx_debug_assert(end() != end_of_storage()) ;
            end_++ ;
        }


        template <int DDIM2> void copy(const small_set<T,DDIM2>& rhs) {
            end_ = data_ ;
            for(typename small_set<T,DDIM2>::const_iterator it=rhs.begin(); it!=rhs.end(); it++) {
                push_back(*it) ;
            }
        }


        // Note: maybe we should start from end() instead of begin()
        // since negative indices are inserted first.
        iterator find_i(const T& x) {
            iterator result = begin() ;
            while(result != end() && *result < x) {
                result++ ;
            }
            return result ;
        }

        const_iterator find_i(const T& x) const {
            const_iterator result = begin() ;
            while(result != end() && *result < x) {
                result++ ;
            }
            return result ;
        }

    protected:
        T data_[DDIM] ;
        iterator end_ ;
    } ;

    template <class T, int DDIM> inline std::ostream& operator<<(
        std::ostream& out,
        const small_set<T, DDIM>& S
    ) {
        S.print(out) ; return out ;
    }


    template <class T, int DDIM1, int DDIM2, int DDIM3> inline void sets_intersect(
        const small_set<T,DDIM1>& S1, const small_set<T,DDIM2>& S2, small_set<T,DDIM3>& I
    ) {
        I.clear() ;
        typename small_set<T,DDIM1>::const_iterator i1 = S1.begin() ;
        typename small_set<T,DDIM2>::const_iterator i2 = S2.begin() ; 
        while(i1 < S1.end() && i2 < S2.end()) {
            if(*i1 < *i2) 
                ++i1 ;
            else if (*i2 < *i1) 
                ++i2 ;
            else {
                I.push_back(*i1) ;
                ++i1 ;
                ++i2 ;
            }
        }
    }


    template <class T, int DDIM> inline std::ostream& operator<<(
        std::ostream& out,
        const std::set<T>& S
    ) {
        out << "[ " ;
        for(typename std::set<T>::const_iterator it = S.begin(); it != S.end(); it++) {
            out << *it << " " ;
        }
        out << "]" ;
        return out ;
    }
    
    template <class T> inline void sets_intersect(
        const std::set<T>& S1, const std::set<T>& S2, std::set<T>& I
    ) {
        I.clear() ;
        std::set_intersection(
            S1.begin(), S1.end(), S2.begin(), S2.end(), std::inserter(I, I.begin())
        ) ;
    }

	//--------------------------------------------------------------------------
}

#endif