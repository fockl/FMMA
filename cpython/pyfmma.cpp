#include<cstdio>
#include "../include/fmma/fmma.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

namespace py = pybind11;

namespace fmma{

template<typename TYPE, std::size_t DIM>
class pyFMMA : public FMMA<TYPE, DIM>{
  public:

  void pysolve(const py::array_t<TYPE>& pytarget, const py::array_t<TYPE>& pysource_weight, const py::array_t<TYPE>& pysource, py::array_t<TYPE>& pyans){
    py::buffer_info buf_target = pytarget.request();
    py::buffer_info buf_source_weight = pysource_weight.request();
    py::buffer_info buf_source = pysource.request();
    py::buffer_info buf_ans = pyans.request();

    if(buf_target.ndim != 2){
      std::string str = "target shape must be 2 but " + std::to_string(buf_target.ndim);
      throw std::runtime_error(str.c_str());
    }
    if(buf_source_weight.ndim != 1){
      std::string str = "source_weight shape must be 2 but " + std::to_string(buf_source_weight.ndim);
      throw std::runtime_error(str.c_str());
    }
    if(buf_source.ndim != 2){
      std::string str = "source shape must be 2 but " + std::to_string(buf_source.ndim);
      throw std::runtime_error(str.c_str());
    }
    if(buf_ans.ndim != 1){
      std::string str = "ans shape must be 1 but " + std::to_string(buf_ans.ndim);
      throw std::runtime_error(str.c_str());
    }
    if(buf_target.shape[1] != DIM){
      std::string str = "target dimension must be " + std::to_string(DIM) + " but " + std::to_string(buf_target.shape[1]);
      throw std::runtime_error(str.c_str());
    }
    if(buf_source.shape[1] != DIM){
      std::string str = "source dimension must be " + std::to_string(DIM) + " but " + std::to_string(buf_source.shape[1]);
      throw std::runtime_error(str.c_str());
    }
    if(buf_ans.shape[0] != buf_target.shape[0]){
      std::string str = "ans num " + std::to_string(buf_ans.shape[0]) + " is different from target num " + std::to_string(buf_target.shape[0]);
      throw std::runtime_error(str.c_str());
    }
    if(buf_source.shape[0] != buf_source_weight.shape[0]){
      std::string str = "source num " + std::to_string(buf_source.shape[0]) + " is different from source_weight num " + std::to_string(buf_source_weight.shape[0]);
      throw std::runtime_error(str.c_str());
    }

    std::vector<std::array<TYPE, DIM>> target(buf_target.shape[0]), source(buf_source.shape[0]);
    std::vector<TYPE> source_weight(buf_source_weight.shape[0]), ans(buf_target.shape[0]);

    TYPE *pytarget_ptr = static_cast<TYPE*>(buf_target.ptr);
    TYPE *pysource_ptr = static_cast<TYPE*>(buf_source.ptr);
    TYPE *pysource_weight_ptr = static_cast<TYPE*>(buf_source_weight.ptr);
    for(std::size_t i=0; i<buf_target.shape[0]; ++i){
      for(std::size_t j=0; j<DIM; ++j){
        target[i][j] = pytarget_ptr[i*DIM+j];
      }
    }
    for(std::size_t i=0; i<buf_source.shape[0]; ++i){
      for(std::size_t j=0; j<DIM; ++j){
        source[i][j] = pysource_ptr[i*DIM+j];
      }
      source_weight[i] = pysource_weight_ptr[i];
    }

    this->solve(target, source_weight, source, ans);

    TYPE *pyans_ptr = static_cast<TYPE*>(buf_ans.ptr);
    for(std::size_t i=0; i<ans.size(); ++i){
      pyans_ptr[i] = ans[i];
    }
    return;
  }

  void pysolve(const py::array_t<TYPE>& pytarget, const py::array_t<TYPE>& pysource, py::array_t<TYPE>& pyans){
    py::buffer_info buf_target = pytarget.request();
    py::buffer_info buf_source = pysource.request();
    py::buffer_info buf_ans = pyans.request();

    if(buf_target.ndim != 2){
      std::string str = "target shape must be 2 but " + std::to_string(buf_target.ndim);
      throw std::runtime_error(str.c_str());
    }
    if(buf_source.ndim != 2){
      std::string str = "source shape must be 2 but " + std::to_string(buf_source.ndim);
      throw std::runtime_error(str.c_str());
    }
    if(buf_ans.ndim != 1){
      std::string str = "ans shape must be 1 but " + std::to_string(buf_ans.ndim);
      throw std::runtime_error(str.c_str());
    }
    if(buf_target.shape[1] != DIM){
      std::string str = "target dimension must be " + std::to_string(DIM) + " but " + std::to_string(buf_target.shape[1]);
      throw std::runtime_error(str.c_str());
    }
    if(buf_source.shape[1] != DIM){
      std::string str = "source dimension must be " + std::to_string(DIM) + " but " + std::to_string(buf_source.shape[1]);
      throw std::runtime_error(str.c_str());
    }
    if(buf_ans.shape[0] != buf_target.shape[0]){
      std::string str = "ans num " + std::to_string(buf_ans.shape[0]) + " is different from target num " + std::to_string(buf_target.shape[0]);
      throw std::runtime_error(str.c_str());
    }

    std::vector<std::array<TYPE, DIM>> target(buf_target.shape[0]), source(buf_source.shape[0]);
    std::vector<TYPE> ans(buf_target.shape[0]);

    TYPE *pytarget_ptr = static_cast<TYPE*>(buf_target.ptr);
    TYPE *pysource_ptr = static_cast<TYPE*>(buf_source.ptr);
    for(std::size_t i=0; i<buf_target.shape[0]; ++i){
      for(std::size_t j=0; j<DIM; ++j){
        target[i][j] = pytarget_ptr[i*DIM+j];
      }
    }
    for(std::size_t i=0; i<buf_source.shape[0]; ++i){
      for(std::size_t j=0; j<DIM; ++j){
        source[i][j] = pysource_ptr[i*DIM+j];
      }
    }

    this->solve(target, source, ans);

    TYPE *pyans_ptr = static_cast<TYPE*>(buf_ans.ptr);
    for(std::size_t i=0; i<ans.size(); ++i){
      pyans_ptr[i] = ans[i];
    }
    return;
  }
};

} // end namespace fmma

template<typename TYPE, std::size_t DIM>
inline void make_python(py::module &m, const std::string& type_str){
  using Class = fmma::pyFMMA<TYPE, DIM>;
  std::string pyname = "fmma" + type_str + std::to_string(DIM);

  py::class_<Class>(m, pyname.c_str())
    .def(py::init<>())
    .def("solve", static_cast<void (Class::*)(const py::array_t<TYPE>& pytarget, const py::array_t<TYPE>& pysource_weight, const py::array_t<TYPE>& pysource, py::array_t<TYPE>& pyans)>(&Class::pysolve))
    .def("solve", static_cast<void (Class::*)(const py::array_t<TYPE>& pytarget, const py::array_t<TYPE>& pysource, py::array_t<TYPE>& pyans)>(&Class::pysolve))
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
