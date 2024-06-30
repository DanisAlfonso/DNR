#ifndef DNR_H
#define DNR_H

//#define _CHECKBOUNDS_ 1
//#define _USESTDVECTOR_ 1
//#define _USENRERRORCLASS_ 1
//#define _TURNONFPES_ 1

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fcntl.h>
#include <cstring>
#include <cctype>
#include <stdexcept>
#include <concepts>
#include <utility>
#include <memory>

namespace DNR {
// macro-like inline functions

    template<typename T>
    requires std::is_arithmetic_v<T>
    inline T SQR(const T a) { return a * a; }

    template<typename T>
    requires std::is_arithmetic_v<T>
    inline const T& MAX(const T& a, const T& b) { return b > a ? b : a; }

    inline float MAX(const double& a, const float& b) { return b > a ? b : float(a); }
    inline float MAX(const float& a, const double& b) { return b > a ? float(b) : a; }

    template<typename T>
    requires std::is_arithmetic_v<T>
    inline const T& MIN(const T& a, const T& b) { return b < a ? b : a; }

    inline float MIN(const double& a, const float& b) { return b < a ? b : float(a); }
    inline float MIN(const float& a, const double& b) { return b < a ? float(b) : a; }

    template<typename T>
    requires std::is_arithmetic_v<T>
    inline T SIGN(const T& a, const T& b) { return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

    inline float SIGN(const float& a, const double& b) { return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }
    inline float SIGN(const double& a, const float& b) { return float(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a)); }

    template<typename T>
    inline void SWAP(T& a, T& b) { std::swap(a, b); }

// exception handling

#ifndef _USEERRORCLASS_
#define throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw std::runtime_error(message);}
#else
struct error {
    const char* message;
    const char* file;
    int line;
    error(const char* m, const char* f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(error(message,__FILE__,__LINE__));
void catch_error(error err) {
    printf("ERROR: %s\n     in file %s at line %d\n",
        err.message, err.file, err.line);
    std::exit(1);
}
#endif

// Vector and Matrix Classes

#ifdef _USESTDVECTOR_
#define vector Vector
#else

    template <class T>
    class Vector {
    private:
        int nn; // size of array. upper index is nn-1
        std::unique_ptr<T[]> v;
    public:
        Vector() : nn(0), v(nullptr) {}
        explicit Vector(int n) : nn(n), v(n > 0 ? std::make_unique<T[]>(n) : nullptr) {}
        Vector(int n, const T& a) : nn(n), v(n > 0 ? std::make_unique<T[]>(n) : nullptr) { std::fill_n(v.get(), n, a); }
        Vector(int n, const T* a) : nn(n), v(n > 0 ? std::make_unique<T[]>(n) : nullptr) { std::copy(a, a + n, v.get()); }
        Vector(const Vector& rhs) : nn(rhs.nn), v(nn > 0 ? std::make_unique<T[]>(nn) : nullptr) { std::copy(rhs.v.get(), rhs.v.get() + nn, v.get()); }
        Vector& operator=(const Vector& rhs) {
            if (this != &rhs) {
                nn = rhs.nn;
                v = nn > 0 ? std::make_unique<T[]>(nn) : nullptr;
                std::copy(rhs.v.get(), rhs.v.get() + nn, v.get());
            }
            return *this;
        }
        typedef T value_type; // make T available externally
        inline T& operator[](const int i) {
#ifdef _CHECKBOUNDS_
            if (i < 0 or i >= nn) {
                throw std::out_of_range("Vector subscript out of bounds");
            }
#endif
            return v[i];
        }
        inline const T& operator[](const int i) const {
#ifdef _CHECKBOUNDS_
            if (i < 0 or i >= nn) {
                throw std::out_of_range("Vector subscript out of bounds");
            }
#endif
            return v[i];
        }
        inline int size() const { return nn; }
        void resize(int newn) {
            nn = newn;
            v = nn > 0 ? std::make_unique<T[]>(nn) : nullptr;
        }
        void assign(int newn, const T& a) {
            nn = newn;
            v = nn > 0 ? std::make_unique<T[]>(nn) : nullptr;
            std::fill_n(v.get(), nn, a);
        }
        ~Vector() = default;
    };

#endif //ifdef _USESTDVECTOR_

    template <typename T>
    class Matrix {
    private:
        int nn;
        int mm;
        std::unique_ptr<T[]> v;
    public:
        Matrix() : nn(0), mm(0), v(nullptr) {}
        Matrix(int n, int m) : nn(n), mm(m), v((n > 0 && m > 0) ? std::make_unique<T[]>(n * m) : nullptr) {}
        Matrix(int n, int m, const T& a) : nn(n), mm(m), v((n > 0 && m > 0) ? std::make_unique<T[]>(n * m) : nullptr) { std::fill_n(v.get(), n * m, a); }
        Matrix(int n, int m, const T* a) : nn(n), mm(m), v((n > 0 && m > 0) ? std::make_unique<T[]>(n * m) : nullptr) { std::copy(a, a + n * m, v.get()); }
        Matrix(const Matrix& rhs) : nn(rhs.nn), mm(rhs.mm), v((nn > 0 && mm > 0) ? std::make_unique<T[]>(nn * mm) : nullptr) { std::copy(rhs.v.get(), rhs.v.get() + nn * mm, v.get()); }
        Matrix& operator=(const Matrix& rhs) {
            if (this != &rhs) {
                nn = rhs.nn;
                mm = rhs.mm;
                v = (nn > 0 && mm > 0) ? std::make_unique<T[]>(nn * mm) : nullptr;
                std::copy(rhs.v.get(), rhs.v.get() + nn * mm, v.get());
            }
            return *this;
        }
        typedef T value_type; // make T available externally
        inline T* operator[](const int i) { return v.get() + i * mm; }
        inline const T* operator[](const int i) const { return v.get() + i * mm; }
        inline int nrows() const { return nn; }
        inline int ncols() const { return mm; }
        void resize(int newn, int newm) {
            nn = newn;
            mm = newm;
            v = (nn > 0 && mm > 0) ? std::make_unique<T[]>(nn * mm) : nullptr;
        }
        void assign(int newn, int newm, const T& a) {
            nn = newn;
            mm = newm;
            v = (nn > 0 && mm > 0) ? std::make_unique<T[]>(nn * mm) : nullptr;
            std::fill_n(v.get(), nn * mm, a);
        }
        ~Matrix() = default;
    };

    template <typename T>
    class Matrix3D {
    private:
        int nn;
        int mm;
        int kk;
        std::unique_ptr<T[]> v;
    public:
        Matrix3D() : nn(0), mm(0), kk(0), v(nullptr) {}
        Matrix3D(int n, int m, int k) : nn(n), mm(m), kk(k), v((n > 0 && m > 0 && k > 0) ? std::make_unique<T[]>(n * m * k) : nullptr) {}
        inline T* operator[](const int i) { return v.get() + i * mm * kk; }
        inline const T* operator[](const int i) const { return v.get() + i * mm * kk; }
        inline int dim1() const { return nn; }
        inline int dim2() const { return mm; }
        inline int dim3() const { return kk; }
        ~Matrix3D() = default;
    };

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

    static const double NaN = std::numeric_limits<double>::quiet_NaN();

    //Uint proto_nan[2]={0xffffffff, 0x7fffffff};
    //double NaN = *( double* )proto_nan;

    //Doub NaN = sqrt(-1.);
// complex
    using Complex = std::complex<double>;

// Vector types

    using VectorInt = Vector<int>;
    using VectorUnsignedInt = Vector<unsigned int>;
    using VectorLongLongInt = Vector<long long int>;
    using VectorUnsignedLongLongInt = Vector<unsigned long long int>;
    using VectorChar = Vector<char>;
    using VectorCharP = Vector<char*>;
    using VectorUnsignedChar = Vector<unsigned char>;
    using VectorDouble = Vector<double>;
    using VectorDoubleP = Vector<double*>;
    using VectorComplex = Vector<Complex>;
    using VectorBool = Vector<bool>;

// matrix types
    using MatrixInt = Matrix<int>;
    using MatrixUnsignedInt = Matrix<unsigned int>;
    using MatrixLongLongInt = Matrix<long long int>;
    using MatrixUnsignedLongLongInt = Matrix<unsigned long long int>;
    using MatrixChar = Matrix<char>;
    using MatrixUnsignedChar = Matrix<unsigned char>;
    using MatrixDouble = Matrix<double>;
    using MatrixBool = Matrix<bool>;

// 3D matrix types
    using Matrix3DDouble = Matrix3D<double>;

}

#endif /* DNR_H */

