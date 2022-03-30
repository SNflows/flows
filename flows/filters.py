import logging

FILTERS = {
    'up': 'u_mag',
    'gp': 'g_mag',
    'rp': 'r_mag',
    'ip': 'i_mag',
    'zp': 'z_mag',
    'B': 'B_mag',
    'V': 'V_mag',
    'J': 'J_mag',
    'H': 'H_mag',
    'K': 'K_mag',
}

FALLBACK_FILTER = 'gp'
logger = logging.getLogger(__name__)


def get_reference_filter(photfilter: str) -> str:
    """
    Translate photometric filter into table column.

    Parameters:
        photfilter (str): photometric filter corresponding to key of FILTERS
    """

    _ref_filter = FILTERS.get(photfilter, None)
    if _ref_filter is None:
        logger.warning(f"Could not find filter {photfilter} in catalogs. "
                       f"Using default {FALLBACK_FILTER} filter.")
        _ref_filter = FILTERS[FALLBACK_FILTER]
    return _ref_filter
