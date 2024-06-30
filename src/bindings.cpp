#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "lu_decomposition.h"

namespace py = pybind11;
using namespace DNR;

PYBIND11_MODULE(DNR, m) {
    m.doc() = "Python bindings for DNRR";

    py::class_<Matrix<double>>(m, "MatrixDouble")
        .def(py::init<int, int>())
        .def("nrows", &Matrix<double>::nrows)
        .def("ncols", &Matrix<double>::ncols)
        .def("__getitem__", [](const Matrix<double> &m, std::pair<size_t, size_t> i) {
            return m[i.first][i.second];
        })
        .def("__setitem__", [](Matrix<double> &m, std::pair<size_t, size_t> i, double v) {
            m[i.first][i.second] = v;
        });

    py::class_<Vector<double>>(m, "VectorDouble")
        .def(py::init<int>())
        .def("size", &Vector<double>::size)
        .def("__getitem__", [](const Vector<double> &v, size_t i) {
            return v[i];
        })
        .def("__setitem__", [](Vector<double> &v, size_t i, double val) {
            v[i] = val;
        });

    py::class_<LU>(m, "LU")
        .def(py::init<const Matrix<double>&>())
        .def("solve", py::overload_cast<const Vector<double>&, Vector<double>&>(&LU::solve))
        .def("solve_matrix", py::overload_cast<const Matrix<double>&, Matrix<double>&>(&LU::solve))
        .def("det", &LU::det)
        .def("inverse", &LU::inverse);
}

