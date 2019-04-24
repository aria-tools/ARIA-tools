//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Author: Ravi Lanka
// Copyright 2019, by the California Institute of Technology. ALL RIGHTS
// RESERVED. United States Government Sponsorship acknowledged.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#include <Python.h>
#include "relaxIVdriver.h"
#include "unwcompmodule.h"

static const char * const __doc__ =
    "unwcomp module for 2-stage unwrapping ";

PyModuleDef moduledef = {
    // header
    PyModuleDef_HEAD_INIT,
    // name of the module
    "unwcomp",
    // module documentation string
    __doc__,
    // size of the per-interpreter state of the module;
    // -1 if this state is global
    -1,
    unwcomp_methods,
};

// initialization function for the module
// *must* be called PyInit_unwcomp
PyMODINIT_FUNC
PyInit_unwcomp()
{
    // create the module using moduledef struct defined above
    PyObject * module = PyModule_Create(&moduledef);
    // check whether module creation succeeded and raise an exception if not
    if (!module) {
        return module;
    }
    // otherwise, we have an initialized module
    // and return the newly created module
    return module;
}

PyObject* relaxIVwrapper_C(PyObject* self, PyObject* args)
{
  char *fileName;

  if (!PyArg_ParseTuple(args, "s", &fileName))
  {
    return NULL;
  }

  // Call MCF with the File Name
  std::vector<int> mcfRet(driver(fileName));

  // Wrap it using binder for returning to python
  PyObject* retVal = PyList_New(0);
  for (int i = 0; i < mcfRet.size() ;i++){
    PyList_Append(retVal, Py_BuildValue("i",mcfRet[i]));
  }

  return retVal;
}
