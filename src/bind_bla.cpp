#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "vector.hpp"

using namespace nanoblas;
namespace py = pybind11;


// Python module nanoblas_impl, imported by nanoblas/__init__.py.
// Each .def maps a Python method (e.g. __getitem__ for v[i]) to a C++ lambda.
PYBIND11_MODULE(nanoblas_impl, m) {
    m.doc() = "Basic linear algebra module"; // optional module docstring
    
    py::class_<Vector<double>> (m, "Vector")
      .def(py::init<size_t>(),
           py::arg("size"), "create vector of given size")

      .def("__len__", &Vector<double>::size,
           "return size of vector")
      
      .def("__setitem__", [](Vector<double> & self, int i, double v) {
        if (i < 0) i += self.size();
        if (i < 0 || i >= self.size()) throw py::index_error("vector index out of range");
        self(i) = v;
      })
      .def("__getitem__", [](Vector<double> & self, int i) { return self(i); })
      
      .def("__setitem__", [](Vector<double> & self, py::slice inds, double val)
      {
        size_t start, stop, step, n;
        inds.compute(self.size(), &start, &stop, &step, &n);
        self.range(start, stop).slice(0,step) = val;
      })
      
      .def("__add__", [](Vector<double> & self, Vector<double> & other)
      { 
        if (self.size() != other.size())
          throw std::runtime_error("Vector sizes do not match: " + std::to_string(self.size()) + " and " + std::to_string(other.size()));
  
        return Vector<double> (self+other); 
      })

      .def("__rmul__", [](Vector<double> & self, double scal)
      { return Vector<double> (scal*self); })
      
      .def("__str__", [](const Vector<double> & self)
      {
        std::stringstream str;
        str << self;
        return str.str();
      })

    ;
}
