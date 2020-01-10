#pragma once

// Must include before Eigen\Dense for additional Quaternion functions
#include "eigen_QuaternionBaseAddons.h"

#ifdef _WIN32
#include <Eigen\Dense>
#else
#include <eigen3/Eigen/Dense>
#endif // _WIN32

#define _USE_MATH_DEFINES
#include <math.h>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>
#include <complex>

#ifndef M_PI
	#define M_PI		3.14159265358979323846
#endif

// Override these defines
#define M_3_PI_2	4.71238898038468985769
#define M_2_PI_		6.28318530717958647692

// Eigen Placeholders
#define NO_VEC2     mg::Vec2::Constant(-1)
#define NO_VEC3     mg::Vec3::Constant(-1)
#define NO_VEC2i    mg::Vec2i::Constant(INT_MIN)
#define NO_VEC3i    mg::Vec3i::Constant(INT_MIN)

#define _i std::complex<double>(0,1)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#if defined(LOG_STDOUT) && defined(LOG_DEBUG)
    #define DLOG(message, ...) printf("[DEBUG] " message, ## __VA_ARGS__)
    #define DLOG_(message, ...) printf("[DEBUG] " message, ## __VA_ARGS__)
#else
    #define DLOG(message, ...)
#endif

#ifdef LOG_STDOUT
#define WLOG(message, ...) printf("[WARN] " message, ## __VA_ARGS__)
#define ELOG(message, ...) printf("[ERROR] " message, ## __VA_ARGS__)
#else
#define WLOG(message, ...)
#define ELOG(message, ...)
#endif

#ifndef NDEBUG
#   define ASSERT(condition, message, ...) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": "; \
            printf(message "\n", ## __VA_ARGS__); \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

namespace mg
{
	template <typename Derived>
	struct Eig_hash { // Pseudo-unique integer hasher
		size_t operator()(const Derived& t) const
		{
			size_t val = 0;
			for (int i = 0; i < _Rows; i++) {
				for (int j = 0; j < _Cols; j++) {
					val = (41 * val) + (unsigned int)t(i, j);
				}
			}
			return val;
		}
	};
    
    template <typename Derived>
    struct Eig_hash2X { // Pseudo-unique integer hasher
        size_t operator()(const Derived& t) const
        {
            return int(41 * t(0)) ^ int(t(1));
        }
    };

    template<typename T, typename K>
    using enum_map = std::unordered_map<T, K, std::hash<int> >;

    template<typename T>
    using enum_set = std::unordered_set<T, std::hash<int> >;

	typedef std::complex<double> cdouble;

	typedef Eigen::Quaternion<double> Quaternion;

	/**
	 * Vec<T,N>
	 **/
	template<typename T, int N>
    using Vec = Eigen::Matrix<T, N, 1>;

	typedef Eigen::Matrix<int, 2, 1> Vec2i;
	typedef Eigen::Matrix<int, 3, 1> Vec3i;

	typedef Eigen::Matrix<double, 2, 1> Vec2;
	typedef Eigen::Matrix<double, 3, 1> Vec3;
	typedef Eigen::Matrix<double, 4, 1> Vec4;
	typedef Eigen::Matrix<double, 5, 1> Vec5;
	typedef Eigen::Matrix<double, -1, 1> VecX;

	typedef Eigen::Matrix<double, 1, 2> Vec2T;
	typedef Eigen::Matrix<double, 1, 3> Vec3T;
	typedef Eigen::Matrix<double, 1, 4> Vec4T;
	typedef Eigen::Matrix<double, 1, -1> VecXT;

	/**
	 * Mat<T,N>
	 **/
	template<typename T, int N, int M>
	using Mat = Eigen::Matrix<T, N, M>;
    
    typedef Eigen::Matrix<double, 2, 2> Mat2;
    typedef Eigen::Matrix<double, 3, 3> Mat3;
    typedef Eigen::Matrix<double, 4, 4> Mat4;
    typedef Eigen::Matrix<double, 2, 3> Matrix23;
    typedef Eigen::Matrix<double, 3, 2> Matrix32;
    typedef Eigen::Matrix<double, 3, 4> Matrix34;
    typedef Eigen::Matrix<double, 4, 3> Matrix43;
    typedef Eigen::Matrix<double, 2, 4> Matrix24;
    typedef Eigen::Matrix<double, -1, -1> MatrixXX;
    
    template<
        typename Derived,
        typename _alloc = Eigen::aligned_allocator<Derived>
    >
        using EigList = std::vector<Derived, _alloc>;

    template<
        typename Derived,
        typename _hash = Eig_hash<Derived>,
        typename _equals = std::equal_to<Derived >,
        typename _alloc = Eigen::aligned_allocator<Derived >
    >
        using EigSet = std::unordered_set<Derived, _hash, _equals, _alloc >;

    template<
        typename T,
        int N,
        int M = 1,
        typename _cmpr = std::less<Eigen::Matrix<T, N, M>>,
        typename _alloc = Eigen::aligned_allocator<Eigen::Matrix<T, N, M> > 
    >
        using EigOrdSet = std::set<Eigen::Matrix<T, N, M>, _cmpr, _alloc >;

    template<
        typename T,
        int N = 2,
        int M = 1,
        typename _hash = Eig_hash2X<T>,
        typename _equals = std::equal_to<Eigen::Matrix<T, N, M> >,
        typename _alloc = Eigen::aligned_allocator<Eigen::Matrix<T, N, M> >
    >
        using EigSet2X = std::unordered_set<Eigen::Matrix<T, N, M>, _hash, _equals, _alloc >;

    template<
        typename KT,
        typename VT,
        typename _hash = std::hash<KT>,
        typename _equals = std::equal_to<KT>,
        typename _alloc = Eigen::aligned_allocator<std::pair<const KT, VT> >
    >
        using EigMap = std::unordered_map<KT, VT, _hash, _equals, _alloc>;

    template<
        typename KT,
        typename VT,
        typename _cmpr = std::less<KT>,
        typename _alloc = Eigen::aligned_allocator<std::pair<const KT, VT> >
    >
        using EigOrdMap = std::map<KT, VT, _cmpr, _alloc>;

	typedef EigList<Vec2> VecList2f;
	typedef EigSet2X<double> VecSet2f;
	typedef EigMap<int, Vec2> VecMap2f;
	typedef EigMap<int, VecSet2f> VecSetMap2f;

	typedef EigList<Vec3> VecList3f;
	typedef EigSet<Vec3> VecSet3f;
	typedef EigMap<int, Vec3> VecMap3f;
	typedef EigMap<int, VecSet3f> VecSetMap3f;
    
    typedef EigList<Vec2i> VecList2i;
    typedef EigSet2X<Vec2i> VecSet2i;
    typedef EigMap<int, Vec2i> VecMap2i;
    typedef EigMap<int, VecList2i> VecSetMap2i;
    
    typedef EigList<Vec3i> VecList3i;
    typedef EigSet<Vec3i> VecSet3i;
    typedef EigMap<int, Vec3i> VecMap3i;
    typedef EigMap<int, VecSet3i> VecSetMap3i;
    
    template<typename _enum, typename VK>
    using enum_VecMap = EigMap<_enum, VK, std::hash<int> >;
    
    template<typename _enum, typename VK>
    using enum_VecSetMap = EigMap<_enum, VK, std::hash<int> >;
    
    template<typename _enum, typename VK>
    using enum_VecSet2i = EigMap<_enum, VK, std::hash<int> >;

    /**
     * Pointer map
     **/
//    template<typename T>
//    struct ptr_deref_hash {
//        inline size_t operator()(const T& pointer) const {
//            return std::hash(*pointer);
//        }
//    };
    
//    template<typename T>
//    struct ptr_deref_equalto {
//        inline bool operator()(const T& lhs, const T& rhs) const {
//            return *lhs == *rhs;
//        }
//    };
    
    template<typename K,typename V>
    using ptr_unordered_map = std::unordered_map<K,V>;

    
//#ifdef _WIN32
//	typedef concurrency::concurrent_unordered_set < Vec2i, std::hash<Vec2i>, std::equal_to<Vec2i>,
//		Eigen::aligned_allocator < Vec2i >> concurrent_VecSet;
//	typedef concurrency::concurrent_unordered_map<int, concurrent_VecSeti> concurrent_VecMap;
//	typedef concurrency::concurrent_unordered_set < Vec2i, std::hash<Vec2i>, std::equal_to<Vec2i>,
//		Eigen::aligned_allocator < Vec2i >> concurrent_VecSeti;
//	typedef concurrency::concurrent_unordered_map<int, concurrent_VecSeti> concurrent_VecMapi;
//#endif

}
