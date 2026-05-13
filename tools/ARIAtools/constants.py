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

# List of supported pixel sizes (in arcseconds)
ARIA_PX_SIZES = [0.000277778, 0.000833334]

# Create lists of all supported models and all aria layers
ARIA_EXTERNAL_CORRECTIONS = [
    "troposphereHydrostatic",
    "troposphereWet",
    "troposphereTotal",
    "solidEarthTide",
]

ARIA_INTERNAL_CORRECTIONS = ["ionosphere"]

ARIA_TROPO_INTERNAL = [
    "ERA5",
    #    'HRES', # HRES model currently not available
    "HRRR",
]

ARIA_TROPO_MODELS = ARIA_TROPO_INTERNAL

ARIA_LAYERS = ["unwrappedPhase", "coherence", "connectedComponents", "amplitude"]
ARIA_LAYERS += ARIA_EXTERNAL_CORRECTIONS
ARIA_LAYERS += ARIA_INTERNAL_CORRECTIONS

ARIA_STANDARD_INTF_LAYERS = ["unwrappedPhase", "coherence"]
ARIA_STANDARD_GEOM_LAYERS = ["incidenceAngle", "azimuthAngle"]
ARIA_STANDARD_LAYERS = ARIA_STANDARD_INTF_LAYERS + ARIA_STANDARD_GEOM_LAYERS


ARIA_STACK_DEFAULTS = [
    "unwrappedPhase",
    "coherence",
    "connectedComponents",
    "troposphereTotal",
    "ionosphere",
    "solidEarthTide",
]

ARIA_STACK_OUTFILES = {
    "unwrappedPhase": "unwrapStack",
    "coherence": "cohStack",
    "connectedComponents": "connCompStack",
    "bParallel": "bParStack",
    "amplitude": "ampStack",
    "troposphereHydrostatic": "tropoHydrostaticStack",
    "troposphereWet": "tropoWetStack",
    "troposphereTotal": "tropoStack",
    "ionosphere": "ionoStack",
    "solidEarthTide": "setStack",
}
ARIA_STACK_OUTFILES.update({i: i + "Stack" for i in ARIA_TROPO_MODELS})
