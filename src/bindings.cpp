#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "math_utilities.h"
#include "lu_decomposition.h"

namespace py = pybind11;

using namespace DNR;

PYBIND11_MODULE(DNR, m) {
    m.doc() = "Python bindings for MyProject";

    m.def("test_function", []() {
        return 42.0;
    }, "A simple test function");

    py::class_<MatrixDouble>(m, "MatrixDouble")
        .def(py::init<int, int>())
        .def("nrows", &MatrixDouble::nrows)
        .def("ncols", &MatrixDouble::ncols)
        .def("__getitem__", [](const MatrixDouble &m, std::pair<int, int> i) {
            return m[i.first][i.second];
        })
        .def("__setitem__", [](MatrixDouble &m, std::pair<int, int> i, double v) {
            m[i.first][i.second] = v;
        });

    py::class_<VectorDouble>(m, "VectorDouble")
        .def(py::init<int>())
        .def("size", &VectorDouble::size)
        .def("__getitem__", [](const VectorDouble &v, int i) {
            return v[i];
        })
        .def("__setitem__", [](VectorDouble &v, int i, double val) {
            v[i] = val;
        });

    py::class_<LU>(m, "LU")
        .def(py::init<const MatrixDouble&>())
        .def("solve", (void (LU::*)(const VectorDouble&, VectorDouble&)) &LU::solve)
        .def("solve_matrix", (void (LU::*)(const MatrixDouble&, MatrixDouble&)) &LU::solve)
        .def("inverse", &LU::inverse)
        .def("det", &LU::det)
        .def("luDecomposition", &LU::luDecomposition);
}

