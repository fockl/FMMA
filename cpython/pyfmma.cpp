#include<cstdio>
#include "../include/fmma/fmma.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

namespace py = pybind11;

template<typename TYPE, std::size_t DIM>
inline void make_python(py::module &m, const std::string& type_str){
  using Class = fmma::FMMA<TYPE, DIM>;
  std::string pyname = "fmma" + type_str + std::to_string(DIM);

  py::class_<Class>(m, pyname.c_str())
    .def(py::init<>())
    .def("solve", static_cast<void (Class::*)(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans)>(&Class::solve))
    .def("solve", static_cast<void (Class::*)(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans)>(&Class::solve))
    .def("check_blas", &Class::check_blas)
    .def("set_fn", static_cast<void (Class::*)(const std::function<TYPE(const std::array<TYPE, DIM>& target_source)>& fn)>(&Class::set_fn))
    .def("set_fn", static_cast<void (Class::*)(const std::function<TYPE(const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source)>& fn)>(&Class::set_fn))
    .def("set_solver_type", &Class::set_solver_type)
    .def("set_nrn_N", &Class::set_nrn_N)
    .def("set_poly_ord", &Class::set_poly_ord)
    .def("set_Depth", &Class::set_Depth)
    .def("set_trans_sym_flag", &Class::set_trans_sym_flag);

  return;
}

PYBIND11_MODULE(pyfmma, m){
  make_python<double, 1>(m, "d");
  make_python<float, 1>(m, "f");
  make_python<double, 2>(m, "d");
  make_python<float, 2>(m, "f");
  make_python<double, 3>(m, "d");
  make_python<float, 3>(m, "f");
}
