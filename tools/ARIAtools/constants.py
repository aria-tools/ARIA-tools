# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author(s): Simran Sangha
# Copyright 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Record global variables for use in other scripts
"""

# Create lists of all supported models and all aria layers
ARIA_EXTERNAL_CORRECTIONS = ['troposphereHydrostatic',
                             'troposphereWet',
                             'troposphereTotal',
                             'solidEarthTide',
                             'gacos_corrections']

ARIA_INTERNAL_CORRECTIONS = ['ionosphere']

ARIA_TROPO_INTERNAL = ['ERA5',
                       'GMAO',
                       'HRES',
                       'HRRR']

ARIA_TROPO_MODELS = ARIA_TROPO_INTERNAL + ['GACOS']

ARIA_LAYERS = ['unwrappedPhase',
               'coherence',
               'connectedComponents',
               'amplitude']
ARIA_LAYERS += ARIA_EXTERNAL_CORRECTIONS
ARIA_LAYERS += ARIA_INTERNAL_CORRECTIONS
ARIA_LAYERS += ARIA_TROPO_MODELS

ARIA_STANDARD_LAYERS = [
    'unwrappedPhase', 'coherence', 'incidenceAngle', 'azimuthAngle']

ARIA_STACK_DEFAULTS = ['unwrappedPhase',
                       'coherence',
                       'connectedComponents',
                       'troposphereTotal',
                       'ionosphere',
                       'solidEarthTide']

ARIA_STACK_OUTFILES = {
    'unwrappedPhase': 'unwrapStack',
    'gacos_corrections': 'gacosStack',
    'coherence': 'cohStack',
    'connectedComponents': 'connCompStack',
    'bParallel': 'bParStack',
    'amplitude': 'ampStack',
    'troposphereHydrostatic': 'tropoHydrostaticStack',
    'troposphereWet': 'tropoWetStack',
    'troposphereTotal': 'tropoStack',
    'ionosphere': 'ionoStack',
    'solidEarthTide': 'setStack'
}
ARIA_STACK_OUTFILES.update({i: i + 'Stack' for i in ARIA_TROPO_MODELS})
