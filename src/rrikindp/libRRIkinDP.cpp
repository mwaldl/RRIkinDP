/**
 * The libRRIkinDP module exposing C++ RRIkinDP functionality to Python
 *
 * This file specifies the Python interface to the low level C++ library.
 * It is compiled into the Python module rrikindp.librrikindp.
 * (C) Maria Waldl 2024
 */


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "RRIkinDP.cpp"
#include <IntaRNA/RnaSequence.h>
#include <IntaRNA/Interaction.h>


namespace py = pybind11;


PYBIND11_MODULE(libRRIkinDP, m) {
    // Bind the EM (energy matrix) class
    py::class_< EM >(m, "EM")
        //.def(py::init<IntaRNA::Interaction, std::string, std::string, std::string, std::string, bool, bool, double>())
        .def(py::init<const IntaRNA::Interaction, const std::string&, const std::string&, const std::string&,const std::string&, bool, bool, double>())
        .def("get_e", &EM::get_e)
        .def("get_hybride_e", &EM::get_hybride_e)
        .def("get_accessibility", &EM::get_accessibility)
        .def("get_ED1", &EM::get_ED1)
        .def("get_ED2", &EM::get_ED2)
        .def("get_minBarrier", &EM::get_minBarrier);
    // Bind IntaRNA RnaSequence
    py::class_<IntaRNA::RnaSequence>(m, "RnaSequence") // IntaRNA::RnaSequence s2(seq_b_id, seq_b);
        .def(py::init<const std::string, const std::string>());
    // Bind the Intarna Interaction class  // IntaRNA::Interaction interaction(s1, s2); interaction.basePairs = bps_list_b0;
    py::class_<IntaRNA::Interaction>(m, "Interaction")
        .def(py::init<const IntaRNA::RnaSequence, const IntaRNA::RnaSequence>())
        .def_readwrite("basePairs", &IntaRNA::Interaction::basePairs)
        .def("is_Valid", &IntaRNA::Interaction::isValid)
        .def("is_Empty", &IntaRNA::Interaction::isEmpty);
}


