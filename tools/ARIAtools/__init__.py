import warnings
import importlib.metadata


try:
    __version__ = importlib.metadata.version('ARIAtools')
except importlib.metadata.PackageNotFoundError:
    __version__ = None
    warnings.warn(
        "ARIA-tools is not installed! Install in editable/develop mode via "
        "(from the top of this repo): python -m pip install -e .",
        RuntimeWarning)
