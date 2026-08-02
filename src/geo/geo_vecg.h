// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/basic/vecg.h
// Copied rather than reimplemented: generic vector template, unchanged

#pragma once

#include "geo_basic.h"
#include "geo_determinant.h"
#include <initializer_list>

#include <iostream>
#include <cfloat>
#include <cmath>

/**
 * \file geogram/basic/vecg.h
 * \brief Generic implementation of geometric vectors
 */

namespace floatTetWild {
namespace geo {

    // Only vecng<3, double> is ever instantiated, so the primary template is declared but not
    // defined: the DIM == 3 specialization below is the whole implementation. The free functions
    // between the two are still written against the primary, as they are in geogram.
    template <index_t DIM, class T>
    class vecng;

    // Compatibility with GLSL

    /**
     * \brief Gets the norm of a vector
     * \param[in] v a vector
     * \return the norm of vector \p v
     * \see vecng::length()
     * \relates vecng
     */
    template <index_t DIM, class T>
    inline T length(const vecng<DIM, T>& v) {
        return v.length();
    }

    /**
     * \brief Gets the square norm of a vector
     * \param[in] v a vector
     * \return the square norm of vector \p v
     * \see vecng::length2()
     * \relates vecng
     */
    template <index_t DIM, class T>
    inline T length2(const vecng<DIM, T>& v) {
        return v.length2();
    }

    /**
     * \brief Gets the square distance between 2 vectors
     * \param[in] v1 the first vector
     * \param[in] v2 the second vector
     * \return the square distance between \p v1 and \p v2.
     * \see vecng::distance2()
     * \relates vecng
     */
    template <index_t DIM, class T>
    inline T distance2(
        const vecng<DIM, T>& v1, const vecng<DIM, T>& v2
    ) {
        return v2.distance2(v1);
    }

    /**
     * \brief Gets the distance between 2 vectors
     * \param[in] v1 the first vector
     * \param[in] v2 the second vector
     * \return the distance between \p v1 and \p v2.
     * \see vecng::distance()
     * \relates vecng
     */
    template <index_t DIM, class T>
    inline T distance(
        const vecng<DIM, T>& v1, const vecng<DIM, T>& v2
    ) {
        return v2.distance(v1);
    }

    /**
     * \brief Normalizes a vector
     * \details Returns a normalized vector constructed by dividing
     * coordinates of vector \p v by it norm. If the norm is 0, then the
     * result is undefined.
     * \param[in] v a vector
     * \return the normalized vector
     * \relates vecng
     */
    template <index_t DIM, class T>
    inline vecng<DIM, T> normalize(
        const vecng<DIM, T>& v
    ) {
        T s = length(v);
        if(s > 1e-30) {
            s = T(1) / s;
        }
        return s * v;
    }






    /************************************************************************/

    /**
     * \brief Specialization of class vecng for DIM == 3
     * \see vecng
     */
    template <class T>
    class vecng<3, T> {
    public:
        /** \copydoc vecng::dim */
        static constexpr index_t dim = 3;

        /** \copydoc vecng::vector_type */
        typedef vecng<dim, T> vector_type;

        /** \copydoc vecng::value_type */
        typedef T value_type;

        /** \copydoc vecng::vecng() */
        vecng() :
            x(T(0.0)),
            y(T(0.0)),
            z(T(0.0)) {
        }

        /**
         * \brief Constructs a vector from coordinates
         * \param[in] x_in , y_in , z_in vector coordinates
         */
        vecng(T x_in, T y_in, T z_in) :
            x(x_in),
            y(y_in),
            z(z_in) {
        }

        /** \copydoc vecng::vecng(const vecng<DIM, T2>&) */
        template <class T2>
        explicit vecng(const vecng<dim, T2>& v) :
            x(v.x),
            y(v.y),
            z(v.z) {
        }

        /** \copydoc vecng::vecng(const T2*) */
        template <class T2>
        explicit vecng(const T2* v) :
            x(v[0]),
            y(v[1]),
            z(v[2]) {
        }

        /** \copydoc vecng::vecng(const std::initializer_list<T>) */
        vecng(const std::initializer_list<T>& Vi) {
            index_t i = 0;
            for(auto& it: Vi) {
                geo_debug_assert(i < dim);
                data()[i] = it;
                ++i;
            }
        }

        /** \copydoc vecng::length2() const */
        inline T length2() const {
            return x * x + y * y + z * z;
        }

        /** \copydoc vecng::length() const */
        inline T length() const {
            return sqrt(x * x + y * y + z * z);
        }

        /** \copydoc vecng::distance2(const vector_type&) const */
        inline T distance2(const vector_type& v) const {
            T dx = v.x - x;
            T dy = v.y - y;
            T dz = v.z - z;
            return dx * dx + dy * dy + dz * dz;
        }

        /** \copydoc vecng::distance(const vector_type&) const */
        inline T distance(const vector_type& v) const {
            return sqrt(distance2(v));
        }

        /** \copydoc vecng::operator+=(const vector_type&) */
        inline vector_type& operator+= (const vector_type& v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        /** \copydoc vecng::operator-=(const vector_type&) */
        inline vector_type& operator-= (const vector_type& v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }

        /** \copydoc vecng::operator*=(T2) */
        template <class T2>
        inline vector_type& operator*= (T2 s) {
            x *= T(s);
            y *= T(s);
            z *= T(s);
            return *this;
        }

        /** \copydoc vecng::operator/=(T2) */
        template <class T2>
        inline vector_type& operator/= (T2 s) {
            x /= T(s);
            y /= T(s);
            z /= T(s);
            return *this;
        }

        /** \copydoc vecng::operator+(const vector_type&) const */
        inline vector_type operator+ (const vector_type& v) const {
            return vector_type(x + v.x, y + v.y, z + v.z);
        }

        /** \copydoc vecng::operator-(const vector_type&) const */
        inline vector_type operator- (const vector_type& v) const {
            return vector_type(x - v.x, y - v.y, z - v.z);
        }

        /** \copydoc vecng::operator*(T2) const */
        template <class T2>
        inline vector_type operator* (T2 s) const {
            return vector_type(x * T(s), y * T(s), z * T(s));
        }

        /** \copydoc vecng::operator/(T2) const */
        template <class T2>
        inline vector_type operator/ (T2 s) const {
            return vector_type(x / T(s), y / T(s), z / T(s));
        }

        /** \copydoc vecng::operator-() const */
        inline vector_type operator- () const {
            return vector_type(-x, -y, -z);
        }

        /** \copydoc vecng::dimension() const */
        index_t dimension() const {
            return dim;
        }

        /** \copydoc vecng::data() */
        T* data() {
            return &x;
        }

        /** \copydoc vecng::data() const */
        const T* data() const {
            return &x;
        }

        /** \copydoc vecng::operator[](index_t) */
        inline T& operator[] (index_t i) {
            geo_debug_assert(i < dim);
            return data()[i];
        }

        /** \copydoc vecng::operator[](index_t) const */
        inline const T& operator[] (index_t i) const {
            geo_debug_assert(i < dim);
            return data()[i];
        }

        /** \brief Optimizes coordinate representation */
        void optimize() {
            Numeric::optimize_number_representation(x);
            Numeric::optimize_number_representation(y);
            Numeric::optimize_number_representation(z);
        }

        /** \brief Vector x coordinate */
        T x;
        /** \brief Vector y coordinate */
        T y;
        /** \brief Vector z coordinate */
        T z;
    };

    /**
     * \copydoc vecng::dot(const vecng<DIM,T>&, const vecng<DIM,T>&)\
     * \relates vecng
     */
    template <class T>
    inline T dot(
        const vecng<3, T>& v1, const vecng<3, T>& v2
    ) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    /**
     * \brief Computes the cross product of 2 vectors
     * \param[in] v1 the first vector
     * \param[in] v2 the second vector
     * \return the cross product (\p v1 x \p v2)
     * \relates vecng
     */
    template <class T>
    inline vecng<3, T> cross(
        const vecng<3, T>& v1, const vecng<3, T>& v2
    ) {
        return vecng<3, T>(
            det2x2(v1.y, v2.y, v1.z, v2.z),
            det2x2(v1.z, v2.z, v1.x, v2.x),
            det2x2(v1.x, v2.x, v1.y, v2.y)
        );
    }

    /**
     * \copydoc vecng::operator*(T2, const vecng<DIM,T>&)
     * \relates vecng
     */
    template <class T2, class T>
    inline vecng<3, T> operator* (
        T2 s, const vecng<3, T>& v
    ) {
        return vecng<3, T>(T(s) * v.x, T(s) * v.y, T(s) * v.z);
    }

    /************************************************************************/

    namespace Numeric {


        template<class T>
        inline void optimize_number_representation(
            vecng<3,T>& v
        ) {
            v.optimize();
        }

    }

    /************************************************************************/

} }

