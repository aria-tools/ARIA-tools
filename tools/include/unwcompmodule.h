//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Author: Ravi Lanka
// Copyright 2019, by the California Institute of Technology. ALL RIGHTS
// RESERVED. United States Government Sponsorship acknowledged.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef snaphumodule_h
#define snaphumodule_h

#include <Python.h>
extern "C"
{
    PyObject *relaxIVwrapper_C(PyObject *self,PyObject *args);
}

static PyMethodDef unwcomp_methods[] = {
    {"relaxIVwrapper_Py",relaxIVwrapper_C,METH_VARARGS," "},
    {NULL,NULL,0,NULL}
};
#endif

