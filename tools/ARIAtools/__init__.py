import warnings

with warnings.catch_warnings():
    from pkg_resources import get_distribution
    __version__ = get_distribution('ARIAtools').version
