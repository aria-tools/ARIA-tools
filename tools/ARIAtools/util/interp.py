# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import scipy.interpolate


class InterpCube(object):
    """Class to interpolate intersection of cube with DEM."""

    def __init__(self, inobj, hgtobj, latobj, lonobj):
        """Init with h5py dataset."""
        self.data = inobj[:]
        self.hgts = hgtobj[:]
        self.offset = None
        self.interp = []
        self.latobj = latobj[:]
        self.lonobj = lonobj[:]

        self.createInterp()

    def createInterp(self):
        """Create interpolators."""
        self.offset = np.mean(self.data)
        for i in range(len(self.hgts)):
            self.interp.append(scipy.interpolate.RectBivariateSpline(
                self.latobj, self.lonobj, self.data[i] - self.offset))

    def __call__(self, line, pix, h):
        """Interpolate at a single point."""
        vals = np.array([x(line, pix)[0, 0] for x in self.interp])
        est = scipy.interpolate.interp1d(self.hgts, vals, kind='cubic')
        return est(h) + self.offset
