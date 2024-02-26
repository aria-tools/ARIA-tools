import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pkg_resources import get_distribution
    try:
        print('ARIA-tools Version:', get_distribution('ARIAtools').version)
    except BaseException:
        pass

