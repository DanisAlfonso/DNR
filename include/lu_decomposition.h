/**************************************************************************
 * Danny Ramírez
 * LU decomposition
 *
 * Example:
 * const int n = ...;
 * MatrixDouble a(n, n);
 * VectorDouble b(n), x(n);
 * ...
 * LU alu(a);
 * alu.solve(b, x);
 **************************************************************************/

#include "math_utilities.h"

namespace DNR {

    /**
     * Object for solving equations Ax = b using LU
     * decomposition
     */
    struct LU {
        int n;
        MatrixDouble lu;
        VectorInt index;
        double d;
        LU(const MatrixDouble& a);

        void luDecomposition(const MatrixDouble& b, MatrixDouble& x);
        void solve(const VectorDouble& b, VectorDouble& x);
        void solve(const MatrixDouble& b, MatrixDouble& x);
        void inverse(MatrixDouble& ainv);
        double det();
        void mprove(const VectorDouble& b, VectorDouble& x);
        const MatrixDouble& aref;
    };

    /**
     * Here is the implementation of the constructor, whose argument
     * is the input matrix that is to be LU decomposed.
     */
    LU::LU(const MatrixDouble& a) : n(a.nrows()), lu(a), aref(a), index(n) {
        const double TINY = 1.0e-40;
        VectorDouble vv(n);
        d = 1.0;
        double big, temp;

        for (int i = 0; i < n; i++) {
            big = 0.0;
            for (int j = 0; j < n; j++) {
                temp = std::abs(lu[i][j]);
                if (temp > big) big = temp;
            }
            if (big == 0.0) throw std::runtime_error("Singular matrix in LU decomposition");
            vv[i] = 1.0 / big;
        }

        for (int k = 0; k < n; k++) {
            big = 0.0;
            int i, j, imax = k;
            for (i = k; i < n; i++) {
                temp = vv[i] * std::abs(lu[i][k]);
                if (temp > big) {
                    big = temp;
                    imax = i;
                }
            }
            if (k != imax) {
                for (j = 0; j < n; j++) {
                    std::swap(lu[imax][j], lu[k][j]);
                }
                d = -d;
                vv[imax] = vv[k];
            }
            index[k]  = imax;
            if (lu[k][k] == 0.0) lu[k][k] = TINY;

            for (i = k + 1; i < n; i++) {
                lu[i][k] /= lu[k][k];
                temp = lu[i][k];
                for (j = k + 1; j < n; j++)
                    lu[i][j] -= temp * lu[k][j];
            }
        }

    }

    /**
     * Solves the set of n linear equations Ax = b using the stored LU
     * decomposition of A.
     */
    void LU::solve(const VectorDouble& b, VectorDouble& x)
    {
        if (b.size() != n || x.size() != n)
            throw std::runtime_error("LU::solve bad sizes");

        for (int i = 0; i < n; i++) x[i] = b[i];

        // Forward substitution
        double sum = 0;
        int ii = 0;
        for (int i = 0; i < n; i++) {
            int ip = index[i];
            sum = x[ip];
            x[ip] = x[i];

            if (ii != 0) {
                for (int j = ii - 1; j < i; j++)  sum -= lu[i][j] * x[j];
            }
            else if (sum != 0.0)
                ii = i + 1;
            x[i] = sum;

        }

        // Now we do the backsubstitution
        for (int i = n - 1; i >= 0; i--) {
            sum = x[i];
            for (int j = i + 1; j < n; j++) sum -= lu[i][j] * x[j];
            x[i] = sum / lu[i][i];
        }
    }

    /**
     * Solves m sets of n linear equations Ax = B using the stored LU
     * decomposition of A. The matrix b[0..n-1][0..m-1] inputs the
     * right-hand sides, while x[0..n-1][0..m-1] returns the solution
     * the inverse of A times B.
     */
    void LU::solve(const MatrixDouble& b, MatrixDouble& x)
    {
        if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
            throw std::runtime_error("LU::solve bad sizes");

        VectorDouble xx(n);
        for (int j = 0; j < b.ncols(); j++) {
            for (int i = 0; i < n; i++) xx[i] = b[i][j];
            solve(xx, xx);
            for (int i = 0; i < n; i++) x[i][j] = xx[i];
        }
    }

    /**
     * Inverse of a Matrix.
     * LU has a number function that gives the inverse of the matrix A.
     * Simply, it creates an identity matrix and then invokes the
     * appropriate solve method.
     */
    void LU::inverse(MatrixDouble& ainv)
    {
        ainv.resize(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) ainv[i][j] = 0.;
            ainv[i][i] = 1.;
        }
        solve(ainv, ainv);
    }

    /**
     * Using the stored LU decomposition, return the determinant
     * of the matrix A.
     */
    double LU::det()
    {
        double dd = d;
        for (int i = 0; i < n; i++) dd *= lu[i][i];
        return dd;
    }

    /**
     * The stored LU decomposition.
     */
    void LU::luDecomposition(const MatrixDouble& b, MatrixDouble& x)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                x[i][j] = lu[i][j];
        }
    }
}

